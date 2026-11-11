
### add and decrement

#struct Newmarkβdecrement!{OX} end
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

# reset existing `state` to `initialstate`
function reset!(state::State{1,NDX,1},initialstate::State) where{NDX}
    state.time = initialstate.time
    for i=1:min(NDX,length(initialstate.X))
        state.X[i]     .= initialstate.X[i]
    end
    for i= length(initialstate.X)+1:NDX
        state.X[i]     .= 0.
    end
end

"""
    sweeps = [Sweep{OX}(time=0:.1:1),...]

NB: `Sweep` is not a solver, just a datastructure to
specify input for [`SweepXA`](@ref).    

# Type parameter
- `OX` 0, 1, 2  order of time derivatives on `X` (static, "viscous", dynamic)    

# Named arguments to `Sweep`:
- `time`             an `AbstractVector` of the times at which to step
- `dependson`        Describes which state to use to start the sweep
                     - A positive number `j` if the sweep depends on the last step
                     of sweeps[j].
                     - A negative number `-j` if the sweeps depends on `initialstate[-j]`
                     (See [`SweepXA`](@ref))
- `active= @functor(with(),active(t)=true)` a `Bool`-valued [`Functor`](@ref) of time::`ℝ`. Must be `true`
                     for steps that contribute to the A-cost-knowing-equilibrium.
- `maxXiter=50`      set to false to suppress printed output (for testing)
- `β=1/4`            Newmark-β parameter, dummy for `OX<2`
- `γ=1/2`            Newmark's γ, parameter dummy for `OX<2`
- `γ=1/2`            Newmark's γ, parameter dummy for `OX<2`
- `maxΔx=1e-5`       Convergence requirement: norm of the increment
- `maxLλ=∞`          Convergence requirement: norm of the residual

See also: [`SweepXA`](@ref)  

"""

struct Sweep{OX,N,C}
    isweep    :: 𝕫
    time      :: 𝕣1
    dependson :: 𝕫 # >0: sweep number <0: initialstate number
    active    :: Functor{N,C}
    maxXiter  :: 𝕫
    β         :: 𝕣
    γ         :: 𝕣
    cΔX²      :: 𝕣
    cLλ²      :: 𝕣
end
@functor with() evertrue(t)=true
Sweep{OX}(;isweep::𝕫,time::ℝ1,dependson::𝕫,active::Functor{N,C}=evertrue,maxXiter::𝕫=50,β::𝕣=1/4,γ::𝕣=1/2,maxΔx::𝕣=1e-5,maxLλ::𝕣=∞) where{OX,N,C} = 
       Sweep{OX,N,C}(isweep,collect(time),dependson,active,maxXiter,β,γ,maxΔx^2,maxLλ^2)
order(::Sweep{OX}) where{OX} = OX

struct SweepData{OX,T}
    data::T  # TODO explicit content
end
SweepData{OX}(data::T) where{OX,T} = SweepData{OX,T}(data) 

function makeSweepData(initialstate,sweep) 
    scheduler             = SingleDependencyScheduler([max(0,sweepᵢ.dependson) for sweepᵢ∈sweep])   
    nsweep                = length(sweep)
    sweepdata             = Vector{SweepData                        }(undef,nsweep) # deliberate abstract eltype
    states                = Vector{Vector{State{1,NDX,1}} where{NDX}}(undef,nsweep) # deliberate abstract eltype
    model,dis             = initialstate[1].model,initialstate[1].dis
    while true # traverse the sweeps in scheduled order
        isweep,_,alldone  = next!(scheduler)
        alldone && break
        done!(scheduler,isweep)
        sweepᵢ            = sweep[isweep]
        j                 = sweepᵢ.dependson
        init              = j>0 ? states[+j][end] : initialstate[-j] # thanks to scheduler, states[j][end] is def'ed
        OX                = order(sweepᵢ)
        outX,asmX,   Xdofgr           = prepare(AssemblySweepX{OX},model,dis                              ) # assembler for std forward analysis
        outA,asmA,(_,Xdofgr,_,Adofgr) = prepare(AssemblyDirect    ,model,dis,Wanted{1,OX+1,1,1}(:all,:all)) # assembler for A-update 
        nXdof,nAdof       = getndof(Xdofgr),getndof(Adofgr)
        buffer            = NTuple{6}(𝕣1(undef,nXdof) for i=1:6)  
        Xₐ                = NTuple{OX+1}(𝕣2(undef,nXdof,nAdof) for idx=1:OX+1)
        Qa,Qaa            = 𝕣1(undef,nAdof), 𝕣2(undef,nAdof,nAdof)
        states[isweep]    = Vector{State{1,OX+1,1}}(undef,length(sweepᵢ.time)) 
        sweepdata[isweep] = SweepData{OX}((;isweep,model,dis,sweep=sweepᵢ,outX,asmX,outA,asmA,initialstate=init,states=states[isweep],
                                    Adofgr,Xdofgr,cXiter=Ref(0),buffer,Xₐ,Qa,Qaa))
    end

    return sweepdata,states,scheduler
end

function addQa!(accum,sweepdata,dbg)
    accum.Qa  .+= sweepdata.data.Qa
    accum.Qaa .+= sweepdata.data.Qaa
end

function sweep!(sweepdata::SweepData{OX},error_ch,dbg;iAiter,verbose,kwargs...) where{OX}  
    local ◺Lλx, ◺K♯
       
    o = sweepdata.data
    zero!(o.Xₐ )
    zero!(o.Qa )
    zero!(o.Qaa)

    state            = State{1,OX+1,1}(copy(o.initialstate))
    state.Λ[1]      .= 0. 

#    @dbg verbose o.sweep.time
    for (istep,t)    ∈ enumerate(o.sweep.time)
        if isready(error_ch)
            verbose && println("Sweep interrupted")
            return
        end

       # verbose && @printf "        step %3d\n" istep
        oldt         = state.time
        state.time   = t
        Δt           = t-oldt
        if Δt≤0 
            if istep==1 muscadeerror(@sprintf("Time step length not strictly positive at istep=1 (time[1]=%g, initialstate.time=%g, Δt=%g)",o.time[1],o.initialstate.time,Δt))
            else        muscadeerror(@sprintf("Time step length not strictly positive at istep=%3d (time[istep]=%g, time[istep-1]=%g, Δt=%g)",istep,o.time[istep],o.time[istep-1],Δt))
            end
        end
        o.outX.c = c = Newmarkβcoefficients{OX}(Δt,o.sweep.β,o.sweep.γ)

        # step and iterations

     #   verbose && @printf "            Equilibrium iterations\n"
        for iXiter   ∈ 1:o.sweep.maxXiter
            o.cXiter[]  += 1
            firstXiter = iXiter==1 
            if   firstXiter assemble!{:step}(o.outX,o.asmX,o.dis,o.model,state,Δt,(;dbg...,solver=:SweepXA,mission=:step,istep,iXiter))
            else            assemble!{:iter}(o.outX,o.asmX,o.dis,o.model,state,Δt,(;dbg...,solver=:SweepXA,mission=:iter,istep,iXiter))
            end
            try if istep==1 && firstXiter 
                    ◺Lλx = lu(       o.outX.Lλx) 
            else           lu!(◺Lλx, o.outX.Lλx) 
            end catch;    muscadeerror(@sprintf("matrix factorization failed at Newmark-β, iAiter=%i, step=%i, iXiter=%i",iAiter,istep,iXiter)) end
            ΔX       = ◺Lλx\o.outX.Lλ
            Newmarkβdecrement!{OX}(state,ΔX ,o.Xdofgr,o.outX.c,firstXiter,o.buffer...)

            ΔX²,Lλ²  = sum(ΔX.^2),sum(o.outX.Lλ.^2)
            if ΔX²≤o.sweep.cΔX² && Lλ²≤o.sweep.cLλ² 
                verbose && @printf "        sweep %3d, step %4d converged in %3d iterations. |Δx|=%7.1e |Lλ|=%7.1e\n" o.isweep istep iXiter √(ΔX²) √(Lλ²)
                o.states[istep] = State(state.time,state.Λ,deepcopy(state.X),state.U,state.A,state.SP,o.model,o.dis) 
                break#out of the iXiter loop
            end
            iXiter==o.sweep.maxXiter && muscadeerror(@sprintf("no convergence of step %3d for iAiter %3d after %3d iterations |Δx|=%g / %g, |Lλ|=%g / %g",istep,iAiter,iXiter,√(ΔX²),√(o.sweep.cΔX²),√(Lλ²)^2,√(o.sweep.cΔλ²)))
        end

        if o.sweep.active(t) && istep>1
       #     verbose && @printf "            Assembling matrices for optimisation phase\n" 

            assemble!{:matrices}(o.outA,o.asmA,o.dis,o.model,state,Δt,(;dbg...,solver=:SweepXA,phase=:Asensitivity,iAiter,istep))
            Lx,La,Lxx,Lax,Laa = o.outA.L1[ind.X], o.outA.L1[ind.A], o.outA.L2[ind.X,ind.X], o.outA.L2[ind.A,ind.X], o.outA.L2[ind.A,ind.A]
            Lλx,Lλa           = o.outA.L2[ind.Λ,ind.X],o.outA.L2[ind.Λ,ind.A]

            # sensitivity of X to A

      #      verbose && @printf "            Evaluate sensitivity of X to A\n" 
            if istep>1
                a₁,a₂,a₃,b₁,b₂,b₃  = c.a₁,c.a₂,c.a₃,c.b₁,c.b₂,c.b₃
                if OX≥2 b₂♯,b₃♯    = b₂/(1-a₂), a₃/(1-a₂)+b₃       end
                if OX≥1 o.Xₐ[2]   .= (1-a₂)  .* o.Xₐ[2] .- a₃      .* o.Xₐ[3]   end #         xₐ′-= aₐ
                if OX≥2 o.Xₐ[3]   .= -b₂♯    .* o.Xₐ[2] .+ (1-b₃♯) .* o.Xₐ[3]   end # same as xₐ″-= bₐ but in place

                Rₐ                 = Matrix(Lλa[1,1]) 
                Rₐ               .+= Lλx[1,1]*o.Xₐ[1]
                if OX≥1 Rₐ       .+= Lλx[1,2]*o.Xₐ[2] end 
                if OX≥2 Rₐ       .+= Lλx[1,3]*o.Xₐ[3] end
                K♯                 = Lλx[1,1]
                if OX≥1 K♯.nzval .+= Lλx[1,2].nzval*c.a₁ end # operate on nzval directly, lest Sparse resize it
                if OX≥2 K♯.nzval .+= Lλx[1,3].nzval*c.b₁ end 

                try if istep==2 
                        ◺K♯ = lu(      K♯) 
                else          lu!(◺K♯, K♯) 
                end catch;    muscadeerror(@sprintf("matrix factorization failed at Xₐ evaluation, iAiter=%i, step=%i",iAiter,istep)) end
                ΔXₐ       = ◺K♯\Rₐ
                Newmarkβdecrement!{OX}(o.Xₐ,ΔXₐ ,o.Xdofgr,o.Adofgr,c) # increments without a or b

                # time-integrate total gradient and Hessian of L wrt A

         #       verbose && @printf "            Time-integrate total gradient and Hessian of L wrt A\n" 
                o.Qa         .+= La[1]             # this is gradient of cost over time Δt
                o.Qaa        .+= Laa[1,1]
                for idx ∈ 1:OX+1
                    o.Qa     .+=               Lx[     idx] ∘₁ o.Xₐ[idx]
                    o.Qaa    .+=               Lax[1  ,idx] ∘₁ o.Xₐ[idx]
                    o.Qaa    .+= o.Xₐ[idx]' ∘₁ Lax[1  ,idx]'  
                    for jdx ∈ 1:OX+1  
                        o.Qaa.+= o.Xₐ[idx]' ∘₁ Lxx[idx,jdx] ∘₁ o.Xₐ[jdx] 
                    end
                end
            end
        end
    end # istep
    for istep = 1:2
        println("218: states[isweep=$(o.isweep)][istep=$istep]")
        state_memory(o.states[istep])
    end

end # sweep and A-sensitivity


### The solver

#=
REPRISE
1) audit memory layout of `state`
2) am I incrementing A at all? Check what line 317 does
=#



"""
	SweepXA

A solver that iterates over the sequence    
1.  non-linear, time domain solution, time-step by time-step of the the `X`-problem,
    leaving  `U`-dofs and `A`-dofs unchanged (see [`SweepX`](@ref)).
2.  adjust `A`-dofs to minimize the total cost of the model.
until convergence.

This solver is not suitable to get a numerical solution to shadow measurement data
(see [`DirectXUA`](@ref) or [`FreqXU`](@ref)).  It is intended as a design optimisation
tool.

The solver allows to handle multiple load cases. A single `initialstate` is provided.
The input allows to describe a series of sweeps (essentialy: [`SweepX`](@ref) analyses).
Each sweep has its own `OX` (see [`SweepX`](@ref)), and an initial state that can be
either `initialstate`, or the state at the last step of another sweep. These dependencies
must form a tree structure, cyclic dependencies are an error.

The following is a simple example:

```
initialstates    = [initialize!(model)]
sweeps          = [Sweep{0}(time=[0,1,2]   ,dependson=-1),   # 1 depends on initialstate[1]
                   Sweep{2}(time=10:.1:12  ,dependson= 1),   # 2 depends on sweeps[1]   
                   Sweep{0}(time=[20,21,22],dependson=-1),   # 3 depends on initialstate[1]
                   Sweep{2}(time=30:.1:32) ,dependson= 3)]   # 4 depends on sweeps[3]
states          = solve(SweepXA;initialstates,sweeps)
```

# Named arguments to `solve`:
- `dbg=(;)`           a named tuple to trace the call tree (for debugging)
- `verbose=true`      set to false to suppress printed output (for testing)
- `silenterror=false` set to true to suppress print out of error (for testing) 
- `initialstates`     a `Vector` of `State`s, obtained from `initialize!` or `SweepX`, that
                      serves as an initial condition (`X`-dofs), seed (`A`-dofs)
                      and immutable value (`U`-dofs).
- `sweeps`            a `Vector` of [`Sweep`](@ref).
- `maxAiter=50`       maximum number of A-iterations
- `maxΔa=1e-5`        convergence criteria: norm of `A`. 
- `maxLa=∞`           convergence criteria: norm of the residual. 
- `multithread=true`  set to false to force sequential computations

# Output

A data structure accessed as `state[isweep][istep]`.                        

See also: [`solve`](@ref), [`initialize!`](@ref), [`findlastassigned`](@ref), [`SingleDependencyScheduler`](@ref), [`Sweep`](@ref)    
"""
struct                SweepXA <: AbstractSolver end
function solve(::Type{SweepXA},pstate,verbose,dbg;
               initialstates::AbstractVector{STATE},
               sweeps::AbstractVector{SWEEP},
               maxAiter::ℤ=50,maxΔa::ℝ=1e-5,maxLa::ℝ=∞,multithread=true) where{STATE<:State,SWEEP<:Sweep}

    cΔA²,cLa²                      = maxΔa^2,maxLa^2
    sweepdata,pstate[],scheduler   = makeSweepData(initialstates,sweeps)
    o1                             = sweepdata[1].data

    nAdof            = getndof(o1.Adofgr)
    Qa               = 𝕣1(undef,nAdof      )
    Qaa              = 𝕣2(undef,nAdof,nAdof)

    local  ΔA, cAiter   # declare variables to be local to the present scope, i.e. scope the function.  
                        # This way we can initialise within a loop, but have a variable that exists at next iter/step or even after the loop


    # main part
    for iAiter ∈ 1:maxAiter
        verbose && @printf "\n    A-iteration %3d\n" iAiter
        
        # time-independent Acost
        assembleA!{:matrices}(o1.outA,o1.asmA,o1.dis,o1.model,o1.initialstate,(;dbg...,solver=:SweepXA,phase=:fixedAcost,iAiter))
        Qa               .= o1.outA.L1[ind.A][1]   
        Qaa              .= o1.outA.L2[ind.A,ind.A][1,1] 

        mode = multithread==true ? :parallel : :sequential
        run_all_jobs!{mode}((;Qa,Qaa),sweep!,addQa!,sweepdata,scheduler,false,(;dbg...,solver=:SweepXA,phase=:fixedAcost,iAiter);iAiter,verbose) 

        try ΔA   = Qaa\Qa
        catch;     muscadeerror(@sprintf("matrix factorization failed at A-update, iAiter=%i",iAiter)) end

        decrement!(o1.states[1],1,ΔA,o1.Adofgr) # only once, as all A are ===


        # Aiter convergence
        ΔA²,La²          = sum(ΔA.^2),sum(Qa.^2)
        verbose && @printf "    A-iteration %3d:                                      |ΔA|=%7.1e |Qa|=%7.1e\n" iAiter √(ΔA²) √(La²)
        if ΔA²≤cΔA² && La²≤cLa² 
            verbose && @printf "    SweepXA converged in %3d A-iterations. |ΔA|=%7.1e / %g |La|=%7.1e / %g\n" iAiter √(ΔA²) maxΔa √(La²) maxLa
            break#out of the iAiter loop
        end
        iAiter==maxAiter && muscadeerror(@sprintf("no convergence of SweepXA after %3d A-iterations |ΔA|=%g / %g, |La|=%g / %g",iAiter,√(ΔA²),maxΔa,√(La²)^2,maxLa))
        cAiter = iAiter

        # # reset state to initial conditions for a new step-sweep
        # for sweepdataᵢ∈sweepdata
        #     reset!(sweepdataᵢ.data.states[1],sweepdataᵢ.data.initialstate)
        # end
    end # iAiter
    verbose && @printf "\n    nel=%d, nXdof=%d, nAdof=%d, nstep=%d, ΣnXiter=%d, mean(nXiter)=%d\n" getnele(o.model) getndof(o.Xdofgr) getndof(o.Adofgr) length(o.time) o.cXiter[] o.cXiter[]/length(o.time)/cAiter



    return
end # solve