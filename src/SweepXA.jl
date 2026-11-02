
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


function reset!(state::State{1,NDX,1},initialstate::State) where{NDX}
    state.time = initialstate.time
    for i=1:min(NDX,length(initialstate.X))
        state.X[i]     .= initialstate.X[i]
    end
    for i= length(initialstate.X)+1:NDX
        state.X[i]     .= 0.
    end
end

struct Route{OX,T}
    data::T  # TODO explicit content
end
Route{OX}(data::T) where{OX,T} = Route{OX,T}(data)

function Route{OX}(initialstate,time,params) where{OX}
    model,dis        = initialstate.model,initialstate.dis
    outX,asmX,        Xdofgr                = prepare(AssemblySweepX{OX},model,dis                              ) # assembler for std forward analysis
    outA,asmA,(Λdofgr,Xdofgr,Udofgr,Adofgr) = prepare(AssemblyDirect    ,model,dis,Wanted{1,OX+1,1,1}(:all,:all)) # assembler for A-update 
    state            = State{1,OX+1,1}(copy(initialstate))
    state.Λ[1]      .= 0.
    states           = Vector{typeof(state)}(undef,length(time)) 
#    states           = allocate(pstate,Vector{typeof(state)}(undef,length(time))) 
    nXdof            = getndof(Xdofgr)
    nAdof            = getndof(Adofgr)
    buffer           = NTuple{6}(𝕣1(undef,nXdof) for i=1:6)  # TODO 6?
    cΔX²,cLλ²        = params.maxΔx^2,params.maxLλ^2
    cXiter           = Ref(0)
    Xₐ               = NTuple{OX+1}(𝕣2(undef,nXdof,nAdof) for idx=1:OX+1)

    return Route{OX}((;model,time,dis,outX,asmX,outA,asmA,initialstate,state,states,buffer,Adofgr,Xdofgr,cΔX²,cLλ²,
                      β=params.β,γ=params.γ,maxXiter=params.maxXiter,cXiter,Xₐ))
end

function walk!(Qa,Qaa,route::Route{OX},verbose,iAiter,dbg) where{OX}
    local ◺Lλx, ◺K♯
    # sweep and A-sensitivity
    # TODO later: allocate and zero Qa, Qaa
    #firstAiter = iAiter==1
    o = route.data
    zero!(o.Xₐ)

    for (istep,t)    ∈ enumerate(o.time)
        #verbose && @printf "        step %3d\n" istep
        oldt         = o.state.time
        o.state.time = t
        Δt           = t-oldt
        if Δt≤0 
            if istep==1 muscadeerror(@sprintf("Time step length not strictly positive at istep=1 (time[1]=%g, initialstate.time=%g, Δt=%g)",o.time[1],o.initialstate.time,Δt))
            else        muscadeerror(@sprintf("Time step length not strictly positive at istep=%3d (time[istep]=%g, time[istep-1]=%g, Δt=%g)",istep,o.time[istep],o.time[istep-1],Δt))
            end
        end
        o.outX.c = c = Newmarkβcoefficients{OX}(Δt,o.β,o.γ)

        # step and iterations

        #verbose && @printf "            Equilibrium iterations\n"
        for iXiter   ∈ 1:o.maxXiter
            o.cXiter[]  += 1
            firstXiter = iXiter==1 
            if   firstXiter assemble!{:step}(o.outX,o.asmX,o.dis,o.model,o.state,Δt,(dbg...,solver=:SweepXA,mission=:step,istep=istep,iXiter=iXiter))
            else            assemble!{:iter}(o.outX,o.asmX,o.dis,o.model,o.state,Δt,(dbg...,solver=:SweepXA,mission=:iter,istep=istep,iXiter=iXiter))
            end
            try if istep==1 && firstXiter #&& firstAiter 
                    ◺Lλx = lu(       o.outX.Lλx) 
            else           lu!(◺Lλx, o.outX.Lλx) 
            end catch;    muscadeerror(@sprintf("matrix factorization failed at Newmark-β, iAiter=%i, step=%i, iXiter=%i",iAiter,istep,iXiter)) end
            ΔX       = ◺Lλx\o.outX.Lλ
            Newmarkβdecrement!{OX}(o.state,ΔX ,o.Xdofgr,o.outX.c,firstXiter,o.buffer...)

            ΔX²,Lλ²  = sum(ΔX.^2),sum(o.outX.Lλ.^2)
            if ΔX²≤o.cΔX² && Lλ²≤o.cLλ² 
                verbose && @printf "        step %3d converged in %3d iterations. |Δx|=%7.1e |Lλ|=%7.1e\n" istep iXiter √(ΔX²) √(Lλ²)
                o.states[istep] = State(o.state.time,o.state.Λ,deepcopy(o.state.X),o.state.U,o.state.A,o.state.SP,o.model,o.dis)
                break#out of the iXiter loop
            end
            iXiter==o.maxXiter && muscadeerror(@sprintf("no convergence of step %3d for iAiter %3d after %3d iterations |Δx|=%g / %g, |Lλ|=%g / %g",istep,iAiter,iXiter,√(ΔX²),maxΔx,√(Lλ²)^2,maxLλ))
        end

        #verbose && @printf "            Assembling matrices for optimisation phase\n" 

        assemble!{:matrices}(o.outA,o.asmA,o.dis,o.model,o.state,Δt,(dbg...,solver=:SweepXA,phase=:Asensitivity,iAiter=iAiter,istep=istep))
        Lx,La,Lxx,Lax,Laa = o.outA.L1[ind.X], o.outA.L1[ind.A], o.outA.L2[ind.X,ind.X], o.outA.L2[ind.A,ind.X], o.outA.L2[ind.A,ind.A]
        Lλx,Lλa           = o.outA.L2[ind.Λ,ind.X],o.outA.L2[ind.Λ,ind.A]

        # sensitivity of X to A

        #verbose && @printf "            Evaluate sensitivity of X to A\n" 
        if istep>1
            a₁,a₂,a₃,b₁,b₂,b₃ = c.a₁,c.a₂,c.a₃,c.b₁,c.b₂,c.b₃
            if OX≥2 b₂♯,b₃♯   = b₂/(1-a₂), a₃/(1-a₂)+b₃       end
            if OX≥1 o.Xₐ[2]    .= (1-a₂)  .* o.Xₐ[2] .- a₃      .* o.Xₐ[3]   end #         xₐ′-= aₐ
            if OX≥2 o.Xₐ[3]    .= -b₂♯    .* o.Xₐ[2] .+ (1-b₃♯) .* o.Xₐ[3]   end # same as xₐ″-= bₐ but in place

            Rₐ                 = Matrix(Lλa[1,1]) 
            Rₐ               .+= Lλx[1,1]*o.Xₐ[1]
            if OX≥1 Rₐ       .+= Lλx[1,2]*o.Xₐ[2] end 
            if OX≥2 Rₐ       .+= Lλx[1,3]*o.Xₐ[3] end
            K♯                 = Lλx[1,1]
            if OX≥1 K♯.nzval .+= Lλx[1,2].nzval*c.a₁ end # operate on nzval directly, lest Sparse resize it
            if OX≥2 K♯.nzval .+= Lλx[1,3].nzval*c.b₁ end 

            try if istep==2 #&& firstAiter
                    ◺K♯ = lu(      K♯) 
            else          lu!(◺K♯, K♯) 
            end catch;    muscadeerror(@sprintf("matrix factorization failed at Xₐ evaluation, iAiter=%i, step=%i",iAiter,istep)) end
            ΔXₐ       = ◺K♯\Rₐ
            Newmarkβdecrement!{OX}(o.Xₐ,ΔXₐ ,o.Xdofgr,o.Adofgr,c) # increments without a or b

            # time-integrate total gradient and Hessian of L wrt A

            #verbose && @printf "            Time-integrate total gradient and Hessian of L wrt A\n" 
            Qa         .+= La[1]             # this is gradient of cost over time Δt
            Qaa        .+= Laa[1,1]
            for idx ∈ 1:OX+1
                Qa     .+=             Lx[     idx] ∘₁ o.Xₐ[idx]
                Qaa    .+=             Lax[1  ,idx] ∘₁ o.Xₐ[idx]
                Qaa    .+= o.Xₐ[idx]' ∘₁ Lax[1  ,idx]'  
                for jdx ∈ 1:OX+1  
                    Qaa.+= o.Xₐ[idx]' ∘₁ Lxx[idx,jdx] ∘₁ o.Xₐ[jdx] 
                end
            end
        end
    end # istep
    # TODO put! Qa,Qaa in channel

end # sweep and A-sensitivity


### The solver

"""
	SweepXA{OX}

A solver that iterates over the sequence    
1.  non-linear, time domain solution, time-step by time-step of the the `X`-problem,
    leaving  `U`-dofs and `A`-dofs unchanged.
2.  adjust `A`-dofs to minimize the total cost of the model.
until convergence.

This solver is not suitable to get a numerical solution to shadow measurement data
(see [`DirectXUA`](@ref) or [`FreqXU`](@ref)).  It is intended as a design optimisation
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

    cΔA²,cLa²        = maxΔa^2,maxLa^2
    route            = Route{OX}(initialstate,time,(;β,γ,maxXiter,maxΔx,maxLλ))
    o                = route.data
    pstate[]         = o.states
    nAdof            = getndof(o.Adofgr)
    Qa               = 𝕣1(undef,nAdof      )
    Qaa              = 𝕣2(undef,nAdof,nAdof)
    local  ΔA, cAiter   # declare variables to be local to the present scope, i.e. scope the function.  
                        # This way we can initialise within a loop, but have a variable that exists at next iter/step or even after the loop

    # main part
    for iAiter ∈ 1:maxAiter
        verbose && @printf "    A-iteration %3d\n" iAiter
        
        # time-independent Acost
        assembleA!{:matrices}(o.outA,o.asmA,o.dis,o.model,o.state,(dbg...,solver=:SweepXA,phase=:fixedAcost,iAiter=iAiter))
        Qa               .= o.outA.L1[ind.A][1]   
        Qaa              .= o.outA.L2[ind.A,ind.A][1,1] 

        walk!(Qa,Qaa,route,verbose,iAiter,(dbg...,iAiter=iAiter))

        # TODO take! from Channel and += to Qa, Qaa     
        # update A
        # verbose && @printf "        Increment total cost given A\n" 
        try ΔA   = Qaa\Qa
        catch;     muscadeerror(@sprintf("matrix factorization failed at A-update, iAiter=%i",iAiter)) end
        decrement!(o.state,1,ΔA,o.Adofgr) 

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
        reset!(o.state,o.initialstate)
    end # iAiter
    verbose && @printf "\n    nel=%d, nXdof=%d, nAdof=%d, nstep=%d, ΣnXiter=%d, mean(nXiter)=%d\n" getnele(o.model) getndof(o.Xdofgr) getndof(o.Adofgr) length(o.time) o.cXiter[] o.cXiter[]/length(o.time)/cAiter

    return
end # solve