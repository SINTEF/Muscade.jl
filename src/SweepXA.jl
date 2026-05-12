struct   propagate!{OX} end
function propagate!{OX}(Xₐ,c) where{OX}
    a₁,a₂,a₃,b₁,b₂,b₃ = c.a₁,c.a₂,c.a₃,c.b₁,c.b₂,c.b₃
    if OX≥2 b₂♯,b₃♯   = b₂/(1-a₂), a₃/(1-a₂)+b₃       end
    if OX≥1 Xₐ[2]   .-= a₂  .* Xₐ[2] .+ a₃ .* Xₐ[3]   end #         xₐ′-= aₐ
    if OX≥2 Xₐ[3]   .-= b₂♯ .* Xₐ[2] .+ b₃♯ .*Xₐ[3]   end # same as xₐ″-= bₐ but in place
    return nothing
end

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
	SweepXA{OX,OSX1,OSX2}

A non-linear, time domain solver, that solves the problem time-step by time-step.
Only the `X`-dofs of the model are solved for, while `U`-dofs and `A`-dofs are unchangeδ

- `SweepXA{0,OSX1,OSX2}` is Newton-Raphson. 
- `SweepXA{1,OSX1,OSX2}` is implicit Euler. 
- `SweepXA{2,OSX1,OSX2}` is Newmark-β, with Newton-Raphson iterations.

`OSX1` and `OSX2` refer to the order of time derivatives of `X` actualy used in the evaluation of `X`-costs.
For example, a dynamic problem can have strain-measurement only, allowing to use `OXS1=OSX2=0`.
`Qa♯` is computed using `OSX1`, while `Qaa♯` uses `OSX2`, so `OSX1>OSX2` introduces a pseudo-Newton step
in the update of `A`. This accelerates each iteration, but makes convergence slover.  

IMPORTANT NOTE: Muscade does not allow elements to have state variables, for example, plastic strain,
or shear-free position for dry friction.  Where the element implements such physics, this 
is implemented by introducing the state as a degree of freedom of the element, and solving
for its evolution, *even in a quasi-static problem*, requires the use of `ORDER≥1`.

An analysis is carried out by a call with the following syntax:

```
initialstate    = initialize!(model)
setdof!(initialstate,1.;class=:U,field=:λcsr)
states           = solve(SweepXA{2};initialstate=initialstate,time=0:10)                                     TODO
```
# Named arguments to `solve`:
- `dbg=(;)`           a named tuple to trace the call tree (for debugging)
- `verbose=true`      set to false to suppress printed output (for testing)
- `silenterror=false` set to true to suppress print out of error (for testing) 
- `initialstate`      a `State`, obtain from `ìnitialize!` or `SweepXA`.
- `time`              maximum number of Newton-Raphson iterations 
- `β=1/4`,`γ=1/2`     parameters to the Newmark-β algorithm. Dummy if `OX<2`
- `maxXiter=50`        maximum number of equilibrium iterations at each step.
- `maxΔx=1e-5`        convergence criteria: norm of `X`. 
- `maxLλ=∞`           convergence criteria: norm of the residual. 
- `maxLineIter=50`    Maximum number of iteration in the feasibility line search.
                      set to 0 to skip the line search (not recommended for models
                      with inequality constraints).
- `sfac=0.5`          Parameter in the line search for a feasible point. If a 
                      tentative result is not feasible, backtrack by a factor `sfac`.
                      If still not feasible, backtrack what is left by a factor `sfac`,
                      and so forth, up to `maxLineIter` times.
- `γfac=0.5`          Parameter for feasibility. For an inequality constraint `g(X)`
                      with reaction force `λ`, require `g(X)*λ==γ`, and multiply
                      `γ *= γfac` at each iteration.                            

# Output

A vector of length equal to that of the named input argument `time` containing the states at the time steps.                       

See also: [`solve`](@ref), [`initialize!`](@ref), [`findlastassigned`](@ref), [`study_singular`](@ref), [`DirectXUA`](@ref), [`FreqXU`](@ref)  
"""
struct        SweepXA{OX,OSX1,OSX2} <: AbstractSolver end
function solve(SX::Type{SweepXA{OX,OSX1,OSX2}},pstate,verbose,dbg;
                    time::AbstractVector{𝕣},
                    initialstate::State,
                    β::𝕣=1/4,γ::𝕣=1/2,
                    maxXiter::ℤ=50,maxΔx::ℝ=1e-5,maxLλ::ℝ=∞,
                    maxAiter::ℤ=50,maxΔa::ℝ=1e-5,maxLa::ℝ=∞) where{OX,OSX1,OSX2}

    model,dis        = initialstate.model,initialstate.dis
    vectors          = ((ind.A,      1  ), 
                        NTuple{OSX1+1     }((ind.X,idx         ) for idx=1:OSX1+1             )...)
    matrices         = ((ind.A,ind.A,1,1), 
                        NTuple{OSX1+1    }((ind.A,ind.X,1  ,idx) for idx=1:OSX1+1             )...,
                        NTuple{OSX1+1    }((ind.X,ind.A,idx,1  ) for idx=1:OSX1+1             )...,  
                        NTuple{(OSX2+1)^2}((ind.X,ind.X,idx,jdx) for idx=1:OSX2+1,jdx=1:OSX2+1)...) 
    outX,asmX,        Xdofgr                = prepare(AssemblySweepX{OX},model,dis                              ) # assembler for std forward analysis
    outA,asmA,(Λdofgr,Xdofgr,Udofgr,Adofgr) = prepare(AssemblyDirect    ,model,dis,Wanted{1,OX+1,1,1}(:all,:all)) # assembler for A-update 
    nXdof            = getndof(Xdofgr)
    nAdof            = getndof(Adofgr)
    buffer           = NTuple{6}(𝕣1(undef,nXdof) for i=1:6)  # TODO 6?
    cΔX²,cLλ²        = maxΔx^2,maxLλ^2
    cΔA²,cLa²        = maxΔa^2,maxLa^2
    cXiter           = 0

    state            = State{1,OX+1,1}(copy(initialstate)) 
    state.Λ[1]      .= 0.
    states           = allocate(pstate,Vector{typeof(state)}(undef,length(time))) 

    La♯              = 𝕣1(undef,nAdof      )
    Laa♯             = 𝕣2(undef,nAdof,nAdof)
    Xₐ               = NTuple{OX+1}(𝕣2(undef,nXdof,nAdof) for idx=1:OX+1)
    local ◺Lλx, ◺K♯, ΔA   # declare variables to be local to the present scope, i.e. scope the function.  
                           # This way we can initialise within a loop, but have a variable that exists at next iter/step or even after the loop

    Ra      = 𝕣1(undef,length(time))  ### dbg
    deltaXa = 𝕣1(undef,length(time))  ### dbg
    Xa      = 𝕣1(undef,length(time))  ### dbg
    Va      = 𝕣1(undef,length(time))  ### dbg
    Aa      = 𝕣1(undef,length(time))  ### dbg
    ΔXs     = 𝕣1(undef,length(time))  ### dbg
    LA      = 𝕣1(undef,length(time))  ### dbg
    LAA     = 𝕣1(undef,length(time))  ### dbg

            fig      = Figure(size = (1000,800))
            axeX      = Axis(fig[1,1])



    # main part
    for iAiter ∈ 1:maxAiter
        firstAiter = iAiter==1
        # time-independent Acost
        assembleA!{:matrices}(outA,asmA,dis,model,state,(dbg...,solver=:SweepXA,phase=:fixedAcost,iAiter=iAiter))
        La♯              .= outA.L1[ind.A][1]   
        Laa♯             .= outA.L2[ind.A,ind.A][1,1] 

        zero!(Xₐ)
        zero!(LA)
        zero!(LAA)

        for (istep,t)    ∈ enumerate(time)
            oldt         = state.time
            state.time   = t
            Δt           = t-oldt
            Δt ≤ 0 && OX>0 && muscadeerror(@sprintf("Time step length not strictly positive at step=%3d",istep))
            outX.c = c = Newmarkβcoefficients{OX}(Δt,β,γ)

            # step and iterations

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

                if iXiter == 1
                    ΔXs[istep] = ΔX[1] ### dbg
                end

                ΔX²,Lλ²  = sum(ΔX.^2),sum(outX.Lλ.^2)
                if ΔX²≤cΔX² && Lλ²≤cLλ² 
                    #verbose && @printf "        step %3d converged in %3d iterations. |Δx|=%7.1e |Lλ|=%7.1e\n" istep iXiter √(ΔX²) √(Lλ²)
                    states[istep] = State(state.time,state.Λ,deepcopy(state.X),state.U,state.A,state.SP,model,dis)
                    break#out of the iXiter loop
                end
                iXiter==maxXiter && muscadeerror(@sprintf("no convergence of step %3d for iAiter %3d after %3d iterations |Δx|=%g / %g, |Lλ|=%g / %g",istep,iAiter,iXiter,√(ΔX²),maxΔx,√(Lλ²)^2,maxLλ))
            end

            assemble!{:matrices}(outA,asmA,dis,model,state,Δt,(dbg...,solver=:SweepXA,phase=:Asensitivity,iAiter=iAiter,istep=istep))
            Lx,La,Lxx,Lax,Laa = outA.L1[ind.X], outA.L1[ind.A], outA.L2[ind.X,ind.X], outA.L2[ind.A,ind.X], outA.L2[ind.A,ind.A]
            Lλx,Lλa           = outA.L2[ind.Λ,ind.X],outA.L2[ind.Λ,ind.A]

            # sensitivity of X to A

            propagate!{OX}(Xₐ,c)
            Rₐ           = Matrix(Lλa[1,1]) 
            Rₐ         .+= Lλx[1,1]*Xₐ[1]
            OX≥1 && Rₐ .+= Lλx[1,2]*Xₐ[2]
            OX≥2 && Rₐ .+= Lλx[1,3]*Xₐ[3]
            K♯           = Lλx[1,1]
            OX≥1 && K♯ .+= Lλx[1,2]*c.a₁ 
            OX≥2 && K♯ .+= Lλx[1,3]*c.b₁ 

            try if istep==1 && firstAiter
                    ◺K♯ = lu(      K♯) 
            else          lu!(◺K♯, K♯) 
            end catch;    muscadeerror(@sprintf("matrix factorization failed at Xₐ evaluation, iAiter=%i, step=%i",iAiter,istep)) end
            ΔXₐ       = ◺K♯\Rₐ
            Newmarkβdecrement!{OX}(Xₐ,ΔXₐ ,Xdofgr,Adofgr,c) # increments without a or b

            # time-integrate total gradient and Hessian of L wrt A

            La♯         .+= La[1]      # this is gradient of cost over time Δt
            Laa♯        .+= Laa[1,1]
            for idx ∈ 1:OX+1
                La♯     .+=             Lx[     idx] ∘₁ Xₐ[idx]
                Laa♯    .+=             Lax[1  ,idx] ∘₁ Xₐ[idx]
                Laa♯    .+= Xₐ[idx]' ∘₁ Lax[1  ,idx]'  
                for jdx ∈ 1:OX+1  
                    Laa♯.+= Xₐ[idx]' ∘₁ Lxx[idx,jdx] ∘₁ Xₐ[jdx] 
                end
            end

            LA[istep]         += La[1][1]      # this is gradient of cost over time Δt
            LAA[istep]        += Laa[1,1][1,1]
            for idx ∈ 1:OX+1
                LA[istep]     +=             Lx[     idx][1] * Xₐ[idx][1,1]
                LAA[istep]    +=             Lax[1  ,idx][1,1] * Xₐ[idx][1,1]
                LAA[istep]    += Xₐ[idx][1,1] * Lax[1  ,idx][1,1]  
                for jdx ∈ 1:OX+1  
                    LAA[istep]+= Xₐ[idx][1,1] * Lxx[idx,jdx][1,1] * Xₐ[jdx][1,1] 
                end
            end
            
            deltaXa[istep] =-ΔXₐ[1, 1]  ### dbg   
            Xa[istep]      =  Xₐ[1][1]  ### dbg   
            Va[istep]      =  Xₐ[2][1]  ### dbg   
            Aa[istep]      =  Xₐ[3][1]  ### dbg   

        end # istep

        Xs = [s.X[1][1] for s∈states] ### dbg   
        Vs = [s.X[2][1] for s∈states] ### dbg   

        if true#iAiter==1 ### dbg   
            # fig      = Figure(size = (1000,800))
            # axeX      = Axis(fig[1,1])
            lines!(axeX,time,Xs,color=:black)
        #   lines!(axeX,time,Va,color=:grey)
           #lines!(axeX,time,Ra,color=:blue)
           #lines!(axeX,time,120*deltaXa,color=:green)
           lines!(axeX,time,Xs+ .1*Xa,color=:lightgrey)
           #lines!(axeX,time,deltaXa,color=:blue)
           #lines!(axeX,time,Va*.1,color=:lightblue)
           lines!(axeX,time,LA./1000,color=:green)
           #lines!(axeX,time,LAA./1000,color=:lightgreen)
           #lines!(axeX,time,-10*Aa,color=:cyan)
           #lines!(axeX,time,ΔXs*-10,color=:orange)
           display(fig)
        end

        # update A
 
        try ΔA   = Laa♯\La♯
        catch;     muscadeerror(@sprintf("matrix factorization failed at A-update, iAiter=%i",iAiter)) end
        

        decrement!(state,1,ΔA,Adofgr) 

        @show Laa♯[1,1],La♯[1],ΔA[1],state.A[1]  ### dbg

        # Aiter convergence
        ΔA²,La²          = sum(ΔA.^2),sum(La♯.^2)
       ##### verbose && @printf "    In A-iteration %3d, |ΔA|=%7.1e |La♯|=%7.1e\n" iAiter √(ΔA²) √(La²)
        if ΔA²≤cΔA² && La²≤cLa² 
            verbose && @printf "    SweepXA converged in %3d A-iterations. |ΔA|=%7.1e / %g |La|=%7.1e / %g\n" iAiter √(ΔA²) maxΔa √(La²) maxLa
            break#out of the iAiter loop
        end
        iAiter==maxAiter && muscadeerror(@sprintf("no convergence of SweepXA after %3d A-iterations |ΔA|=%g / %g, |La|=%g / %g",iAiter,√(ΔA²),maxΔa,√(La²)^2,maxLa))

        # reset state to initial conditions for a new step-sweep
        state.time = initialstate.time
        for i=1:min(OX+1,length(initialstate.X))
            state.X[i]     .= initialstate.X[i]
        end
        for i= length(initialstate.X)+1:OX+1
            state.X[i]     .= 0.
        end
    end 
    verbose && @printf "\n    nel=%d, nXdof=%d, nstep=%d, ΣnXiter=%d, mean(nXiter)=%d\n" getnele(model) getndof(Xdofgr) length(time) cXiter cXiter/length(time)

    return
end