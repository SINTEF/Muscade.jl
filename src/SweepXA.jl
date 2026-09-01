
### add and decrement

function Newmarkβdecrement!{OX}(Xₐ::NTuple{NDX,𝕣2},ΔXₐ,Xgr,Agr,c) where{OX,NDX}
    # xₐ    -=    Δxₐ
    # xₐ′   -= a₁*Δxₐ
    # xₐ″   -= b₁*Δxₐ
    f  = (1.,c.a₁,c.b₁)  
    nA = size(Xₐ[1],2)
    for ider = 1:OX+1
        for igA ∈ eachindex(Agr.iA)
            iA,jA,sA = Agr.iA[igA], Agr.jA[igA], 1 / Agr.scaleA[igA] # inverse scaleA, because Xₐ is dX / dA 
            for igX ∈ eachindex(Xgr.iX)
                iX,jX,sX = Xgr.iX[igX], Xgr.jX[igX], Xgr.scaleX[igX]
                Xₐ[ider][iX,iA] -= f[ider] .* ΔXₐ[jX,jA] * sX * sA
            end
        end 
    end
end

### The solver

"""
	SweepXA{OX}

A solver that iterates over the sequence    
1.  non-linear, time domain solution, time-step by time-step of the the `X`-problem,
    leaving  `U`-dofs and `A`-dofs unchanged.
2.  adjust `A`-dofs to minimize the total cost of the model.
until convergence.

This solver is not suitable to get a numerical solution to shadow measurement data
(see [`DirectXUA](@ref) or [`FreqXU`](@ref)).  It is intended as a design optimisation
tool.

- `SweepXA{0}` uses Newton-Raphson. 
- `SweepXA{1}` uses implicit Euler. 
- `SweepXA{2}` uses Newmark-β, with Newton-Raphson iterations.

IMPORTANT NOTE: Muscade does not allow elements to have state variables, for example, plastic strain,
or shear-free position for dry friction.  Where the element implements such physics, this 
is implemented by introducing the state as a degree of freedom of the element, and solving
for its evolution, *even in a quasi-static problem*, requires the use of `ORDER≥1`.

An analysis is carried out by a call with the following syntax:

```
initialstate    = initialize!(model)
setdof!(initialstate,1.;class=:X,field=:tx1,order=1) 
states           = solve(SweepXA{2};initialstate=initialstate,time=0:10)
```
# Named arguments to `solve`:
- `dbg=(;)`           a named tuple to trace the call tree (for debugging)
- `verbose=true`      set to false to suppress printed output (for testing)
- `silenterror=false` set to true to suppress print out of error (for testing) 
- `initialstate`      a `State`, obtain from `initialize!` or `SweepX`, that
                      serves as an initial condition (`X`-dofs), seed (`A`-dofs)
                      and immutable value (`U`-dofs).
- `time`              an `AbstractVector` of strictly increasing times at which to
                      compute a step. `time[1]>initialstate.time` is required.  
- `β=1/4`,`γ=1/2`     parameters to the Newmark-β algorithm. Dummy if `OX<2`
- `maxXiter=50`       maximum number of equilibrium iterations at each step.
- `maxΔx=1e-5`        convergence criteria: norm of `X`. 
- `maxLλ=∞`           convergence criteria: norm of the residual. 
- `maxAiter=50`       maximum number of A-iterations
- `maxΔa=1e-5`        convergence criteria: norm of `A`. 
- `maxLa=∞`           convergence criteria: norm of the residual. 

# Output

A vector of length equal to that of the named input argument `time` containing the states at the time steps.                       

See also: [`solve`](@ref), [`initialize!`](@ref), [`findlastassigned`](@ref), [`study_singular`](@ref), [`DirectXUA`](@ref), [`FreqXU`](@ref)  
"""
struct                SweepXA{OX} <: AbstractSolver end
function solve(::Type{SweepXA{OX}},pstate,verbose,dbg;
               time::AbstractVector{𝕣},
               initialstate::State,
               β::𝕣=1/4,γ::𝕣=1/2,
               maxXiter::ℤ=50,maxΔx::ℝ=1e-5,maxLλ::ℝ=∞,
               maxAiter::ℤ=50,maxΔa::ℝ=1e-5,maxLa::ℝ=∞) where{OX}

    model,dis        = initialstate.model,initialstate.dis
    outX,asmX,        Xdofgr                = prepare(AssemblySweepX{OX},model,dis                              ) # assembler for std forward analysis
    outA,asmA,(Λdofgr,Xdofgr,Udofgr,Adofgr) = prepare(AssemblyDirect    ,model,dis,Wanted{1,OX+1,1,1}(:all,:all)) # assembler for A-update 
    nXdof            = getndof(Xdofgr)
    nAdof            = getndof(Adofgr)
    buffer           = NTuple{6}(𝕣1(undef,nXdof) for i=1:6)  # TODO 6?
    cΔX²,cLλ²        = maxΔx^2,maxLλ^2
    cΔA²,cLa²        = maxΔa^2,maxLa^2
    cXiter           = 0
    cAiter           = 0

    state            = State{1,OX+1,1}(copy(initialstate))
    state.Λ[1]      .= 0.
    states           = allocate(pstate,Vector{typeof(state)}(undef,length(time))) 

    Qa               = 𝕣1(undef,nAdof      )
    Qaa              = 𝕣2(undef,nAdof,nAdof)
    Xₐ               = NTuple{OX+1}(𝕣2(undef,nXdof,nAdof) for idx=1:OX+1)
    local ◺Lλx, ◺K♯, ΔA, cAiter   # declare variables to be local to the present scope, i.e. scope the function.  
                           # This way we can initialise within a loop, but have a variable that exists at next iter/step or even after the loop

    # main part
    for iAiter ∈ 1:maxAiter
        verbose && @printf "    A-iteration %3d\n" iAiter
        firstAiter = iAiter==1
        # time-independent Acost
        assembleA!{:matrices}(outA,asmA,dis,model,state,(dbg...,solver=:SweepXA,phase=:fixedAcost,iAiter=iAiter))
        Qa               .= outA.L1[ind.A][1]   
        Qaa              .= outA.L2[ind.A,ind.A][1,1] 
        zero!(Xₐ)
        for (istep,t)    ∈ enumerate(time)
            #verbose && @printf "        step %3d\n" istep
            oldt         = state.time
            state.time   = t
            Δt           = t-oldt
            if Δt≤0 
                if istep==1 muscadeerror(@sprintf("Time step length not strictly positive at istep=1 (time[1]=%g, initialstate.time=%g, Δt=%g)",time[1],initialstate.time,Δt))
                else        muscadeerror(@sprintf("Time step length not strictly positive at istep=%3d (time[istep]=%g, time[istep-1]=%g, Δt=%g)",istep,time[istep],time[istep-1],Δt))
                end
            end
            outX.c = c = Newmarkβcoefficients{OX}(Δt,β,γ)

            # step and iterations

            #verbose && @printf "            Equilibrium iterations\n"
            for iXiter   ∈ 1:maxXiter
                cXiter  += 1
                firstXiter = iXiter==1 
                if   firstXiter assemble!{:step}(outX,asmX,dis,model,state,Δt,(dbg...,solver=:SweepXA,mission=:step,iAiter=iAiter,istep=istep,iXiter=iXiter))
                else            assemble!{:iter}(outX,asmX,dis,model,state,Δt,(dbg...,solver=:SweepXA,mission=:iter,iAiter=iAiter,istep=istep,iXiter=iXiter))
                end
                try if istep==1 && firstXiter && firstAiter 
                        ◺Lλx = lu(       outX.Lλx) 
                else           lu!(◺Lλx, outX.Lλx) 
                end catch;    muscadeerror(@sprintf("matrix factorization failed at Newmark-β, iAiter=%i, step=%i, iXiter=%i",iAiter,istep,iXiter)) end
                ΔX       = ◺Lλx\outX.Lλ
                Newmarkβdecrement!{OX}(state,ΔX ,Xdofgr,outX.c,firstXiter,buffer...)

                ΔX²,Lλ²  = sum(ΔX.^2),sum(outX.Lλ.^2)
                if ΔX²≤cΔX² && Lλ²≤cLλ² 
                    verbose && @printf "        step %3d converged in %3d iterations. |Δx|=%7.1e |Lλ|=%7.1e\n" istep iXiter √(ΔX²) √(Lλ²)
                    states[istep] = State(state.time,state.Λ,deepcopy(state.X),state.U,state.A,state.SP,model,dis)
                    break#out of the iXiter loop
                end
                iXiter==maxXiter && muscadeerror(@sprintf("no convergence of step %3d for iAiter %3d after %3d iterations |Δx|=%g / %g, |Lλ|=%g / %g",istep,iAiter,iXiter,√(ΔX²),maxΔx,√(Lλ²)^2,maxLλ))
            end

            #verbose && @printf "            Assembling matrices for optimisation phase\n" 
 
            assemble!{:matrices}(outA,asmA,dis,model,state,Δt,(dbg...,solver=:SweepXA,phase=:Asensitivity,iAiter=iAiter,istep=istep))
            Lx,La,Lxx,Lax,Laa = outA.L1[ind.X], outA.L1[ind.A], outA.L2[ind.X,ind.X], outA.L2[ind.A,ind.X], outA.L2[ind.A,ind.A]
            Lλx,Lλa           = outA.L2[ind.Λ,ind.X],outA.L2[ind.Λ,ind.A]

            # sensitivity of X to A

            #verbose && @printf "            Evaluate sensitivity of X to A\n" 
            if istep>1
                a₁,a₂,a₃,b₁,b₂,b₃ = c.a₁,c.a₂,c.a₃,c.b₁,c.b₂,c.b₃
                if OX≥2 b₂♯,b₃♯   = b₂/(1-a₂), a₃/(1-a₂)+b₃       end
                if OX≥1 Xₐ[2]    .= (1-a₂)  .* Xₐ[2] .- a₃      .* Xₐ[3]   end #         xₐ′-= aₐ
                if OX≥2 Xₐ[3]    .= -b₂♯    .* Xₐ[2] .+ (1-b₃♯) .* Xₐ[3]   end # same as xₐ″-= bₐ but in place

                Rₐ                 = Matrix(Lλa[1,1]) 
                Rₐ               .+= Lλx[1,1]*Xₐ[1]
                if OX≥1 Rₐ       .+= Lλx[1,2]*Xₐ[2] end 
                if OX≥2 Rₐ       .+= Lλx[1,3]*Xₐ[3] end
                K♯                 = Lλx[1,1]
                if OX≥1 K♯.nzval .+= Lλx[1,2].nzval*c.a₁ end # operate on nzval directly, lest Sparse resize it
                if OX≥2 K♯.nzval .+= Lλx[1,3].nzval*c.b₁ end 

                try if istep==2 && firstAiter
                        ◺K♯ = lu(      K♯) 
                else          lu!(◺K♯, K♯) 
                end catch;    muscadeerror(@sprintf("matrix factorization failed at Xₐ evaluation, iAiter=%i, step=%i",iAiter,istep)) end
                ΔXₐ       = ◺K♯\Rₐ
                Newmarkβdecrement!{OX}(Xₐ,ΔXₐ ,Xdofgr,Adofgr,c) # increments without a or b

                # time-integrate total gradient and Hessian of L wrt A

                #verbose && @printf "            Time-integrate total gradient and Hessian of L wrt A\n" 
                Qa         .+= La[1]      # this is gradient of cost over time Δt
                Qaa        .+= Laa[1,1]
                for idx ∈ 1:OX+1
                    Qa     .+=             Lx[     idx] ∘₁ Xₐ[idx]
                    Qaa    .+=             Lax[1  ,idx] ∘₁ Xₐ[idx]
                    Qaa    .+= Xₐ[idx]' ∘₁ Lax[1  ,idx]'  
                    for jdx ∈ 1:OX+1  
                        Qaa.+= Xₐ[idx]' ∘₁ Lxx[idx,jdx] ∘₁ Xₐ[jdx] 
                    end
                end
            end
        end # istep

        # update A
        # verbose && @printf "        Increment total cost given A\n" 
        try ΔA   = Qaa\Qa
        catch;     muscadeerror(@sprintf("matrix factorization failed at A-update, iAiter=%i",iAiter)) end
        decrement!(state,1,ΔA,Adofgr) 

        # Aiter convergence
        ΔA²,La²          = sum(ΔA.^2),sum(Qa.^2)
        verbose && @printf "        In A-iteration %3d, |ΔA|=%7.1e |Qa|=%7.1e\n" iAiter √(ΔA²) √(La²)
        if ΔA²≤cΔA² && La²≤cLa² 
            verbose && @printf "    SweepXA converged in %3d A-iterations. |ΔA|=%7.1e / %g |La|=%7.1e / %g\n" iAiter √(ΔA²) maxΔa √(La²) maxLa
            break#out of the iAiter loop
        end
        iAiter==maxAiter && muscadeerror(@sprintf("no convergence of SweepXA after %3d A-iterations |ΔA|=%g / %g, |La|=%g / %g",iAiter,√(ΔA²),maxΔa,√(La²)^2,maxLa))
        cAiter = iAiter

        # reset state to initial conditions for a new step-sweep
        state.time = initialstate.time
        for i=1:min(OX+1,length(initialstate.X))
            state.X[i]     .= initialstate.X[i]
        end
        for i= length(initialstate.X)+1:OX+1
            state.X[i]     .= 0.
        end
    end # iAiter
    verbose && @printf "\n    nel=%d, nXdof=%d, nAdof=%d, nstep=%d, ΣnXiter=%d, mean(nXiter)=%d\n" getnele(model) getndof(Xdofgr) getndof(Adofgr) length(time) cXiter cXiter/length(time)/cAiter

    return
end # solve