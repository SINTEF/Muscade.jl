### Assembler

mutable struct AssemblySweepXA{OX,NDX} <: Assembly
    # up
    Lx        :: 𝕣1  
    La        :: 𝕣1  
    Lλa       :: 𝕣2 
    Lxx       :: Sparse𝕣2 
    Lax       :: 𝕣2 
    Laa       :: 𝕣2 
    # down
    c         :: Newmarkβcoefficients{OX}
end   

function prepare(::Type{AssemblySweepXA{OX}},model,dis) where{OX}
    Xdofgr             = allXdofs(model,dis) 
    Adofgr             = allAdofs(model,dis)
    nXdof  = nΛdof     = getndof(Xdofgr)
    nAdof              = getndof(Adofgr)
    narray,neletyp     = 6,getneletyp(model)
    asm                = Matrix{𝕫2}(undef,narray,neletyp)  # asm[iarray,ieletyp][ieledof,iele]
    Lx                 = asmvec!(view(asm, 1,:),Xdofgr,dis)
    La                 = asmvec!(view(asm, 2,:),Adofgr,dis)
    Lλa                = asmfullmat!(view(asm, 3,:),view(asm,1,:),view(asm,2,:),nΛdof,nAdof)  
    Lxx                = asmmat!(view(asm, 4,:),view(asm,1,:),view(asm,1,:),nXdof,nXdof)
    Lax                = asmfullmat!(view(asm, 5,:),view(asm,2,:),view(asm,1,:),nAdof,nXdof)  
    Laa                = asmfullmat!(view(asm, 6,:),view(asm,2,:),view(asm,2,:),nAdof,nAdof)

    out                = AssemblySweepXA{OX,OX+1}(Lx,La,Lλa,Lxx,Lax,Laa, Newmarkβcoefficients{OX}()) 
    return out,asm,Xdofgr,Adofgr
end
function zero!(out::AssemblySweepXA) # TODO
    zero!(out.Lx )
    zero!(out.La )
    zero!(out.Lλa)
    zero!(out.Lxx)
    zero!(out.Lax)
    zero!(out.Laa)
end

#=        TODO
solver
write specific adiff for ElementCost
SweepXA for order 0 and 1
Multi load cases        
=#


function addin!{:sensitivity}(out::AssemblySweepXA{2},asm,iele,scale,eleobj,Λ,X::NTuple{Nxder,<:SVector{Nx}},U,A::SVector{Na},t,Δt,SP,dbg) where{Nxder,Nx,Na}
    a₁,a₂,a₃,b₁,b₂,b₃ = out.c.a₁,out.c.a₂,out.c.a₃,out.c.b₁,out.c.b₂,out.c.b₃
    x,x′,x″,λ         = ∂0(X),∂1(X),∂2(X),∂0(Λ)
    δΛ,δX,δA          = reδ{2}((;Λ=λ,X=x,A),(;Λ=scale.Λ,X=scale.X,A=scale.A)) 
    iΛ,iX,iA,Nz       = revariate_indices(λ,x,A) 
    vx                = x  +    δX
    vx′               = x′ + a₁*δX   
    vx″               = x″ + b₁*δX  
    L,FB              = getlagrangian(eleobj,λ+δΛ,(vx,vx′,vx″),U,A+δA,t,SP,dbg) # TODO jump over elements with residual
    ∇L                = ∂{2,Nz}(L)
    add_value!(out.Lx , asm[1], iele, ∇L, iX    ;Δt)     
    add_value!(out.La , asm[2], iele, ∇L, iA    ;Δt)             
    add_∂!{1 }(out.Lλa, asm[3], iele, ∇L, iΛ, iA;Δt)    
    add_∂!{1 }(out.Lxx, asm[4], iele, ∇L, iX, iX;Δt)  
    add_∂!{1 }(out.Lax, asm[5], iele, ∇L, iA, iX;Δt)  
    add_∂!{1 }(out.Laa, asm[6], iele, ∇L, iA, iA;Δt)  
end

function addin!{:Acost}(out::AssemblySweepXA,asm,iele,scale,eleobj::Acost,A::SVector{Na},dbg) where{Na} 
    d      = revariate{2}((;A),(;A=scale.A)) # careful: revariate returns a NamedTuple
    ø      = nothing
    C,_    = lagrangian(eleobj,ø,ø,ø,d.A,ø,ø ,dbg)
    ∇ₐC    = ∂{2,Na}(C)
    add_value!(out.La ,asm[2],iele,∇ₐC)
    add_∂!{1 }(out.Laa,asm[6],iele,∇ₐC)
end



"""
	SweepXA{OX}

A non-linear, time domain solver, that solves the problem time-step by time-step.
Only the `X`-dofs of the model are solved for, while `U`-dofs and `A`-dofs are unchangeδ

- `SweepXA{0}` is Newton-Raphson, with feasibility line-search, to handle inequality constraints. 
- `SweepXA{1}` is implicit Euler, with feasibility line-search. 
- `SweepXA{2}` is Newmark-β, with Newton-Raphson iterations and feasibility line search

IMPORTANT NOTE: Muscade does not allow elements to have state variables, for example, plastic strain,
or shear-free position for dry friction.  Where the element implements such physics, this 
is implemented by introducing the state as a degree of freedom of the element, and solving
for its evolution, *even in a quasi-static problem*, requires the use of `ORDER≥1`.

An analysis is carried out by a call with the following syntax:

```
initialstate    = initialize!(model)
setdof!(initialstate,1.;class=:U,field=:λcsr)
states           = solve(SweepXA{2};initialstate=initialstate,time=0:10)
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
struct        SweepXA{OX} <: AbstractSolver end
function solve(SX::Type{SweepXA{OX}},pstate,verbose,dbg;
                    time::AbstractVector{𝕣},
                    initialstate::State,
                    β::𝕣=1/4,γ::𝕣=1/2,
                    maxXiter::ℤ=50,maxΔx::ℝ=1e-5,maxLλ::ℝ=∞,
                    maxAiter::ℤ=50,maxΔa::ℝ=1e-5,maxLa::ℝ=∞) where{OX}

    model,dis        = initialstate.model,initialstate.dis
    outX ,asmX ,Xdofgr          = prepare(AssemblySweepX{ OX},model,dis)  
    outXA,asmXA,Xdofgr,Adofgr   = prepare(AssemblySweepXA{OX},model,dis)  
    nXdof            = getndof(Xdofgr)
    nAdof            = getndof(Adofgr)
    buffer           = ntuple(i->𝕣1(undef,nXdof), 6)  
    cΔX²,cLλ²        = maxΔx^2,maxLλ^2
    cΔA²,cLa²        = maxΔa^2,maxLa^2
    cXiter           = 0

    state            = State{1,OX+1,1}(copy(initialstate)) 
    state.Λ[1]      .= 0.
    states           = allocate(pstate,Vector{typeof(state)}(undef,length(time))) 

    La♯              = 𝕣1(undef,nAdof      )
    Laa♯             = 𝕣2(undef,nAdof,nAdof)
    local Lλx,ΔX # declare Lλx to scope the function, without having to actualy initialize the variable

    # main part
    for iAiter = 1:maxAiter

        assembleA!{:Acost}(outXA,asmXA,dis,model,state,(dbg...,solver=:SweepXA,phase=:Acost,iAiter=iAiter))
        La♯              .= outXA.La   
        Laa♯             .= outXA.Laa  

        # forward sweep
        for (step,t)     ∈ enumerate(time)
            oldt         = state.time
            state.time   = t
            Δt           = t-oldt
            Δt ≤ 0 && OX>0 && muscadeerror(@sprintf("Time step length not strictly positive at step=%3d",step))
            outX.c        = Newmarkβcoefficients{OX}(Δt,β,γ)
            for iXiter   = 1:maxXiter
                cXiter  += 1
                firstiter = iXiter==1
                if   firstiter assemble!{:step}(outX,asmX,dis,model,state,Δt,(dbg...,solver=:SweepXA,step=step,iXiter=iXiter))
                else           assemble!{:iter}(outX,asmX,dis,model,state,Δt,(dbg...,solver=:SweepXA,step=step,iXiter=iXiter))
                end
                try if step==1  && firstiter  Lλx = lu(outX.Lλx) # here we do not write "local Lλx", so we refer to the variable defined outside the loops (we do not shadow Lλx)
                else                          lu!(Lλx, outX.Lλx) 
                end catch;    muscadeerror(@sprintf("matrix factorization failed at step=%i, iXiter=%i",step,iXiter)) end
                ΔX       = Lλx\outX.Lλ
                ΔX²,Lλ²  = sum(ΔX.^2),sum(outX.Lλ.^2)
                Newmarkβdecrement!{OX}(state,ΔX ,Xdofgr,outX.c,firstiter,buffer...)

                if ΔX²≤cΔX² && Lλ²≤cLλ² 
                    #verbose && @printf "        step %3d converged in %3d iterations. |Δx|=%7.1e |Lλ|=%7.1e\n" step iXiter √(ΔX²) √(Lλ²)
                    states[step] = State(state.time,state.Λ,deepcopy(state.X),state.U,state.A,state.SP,model,dis)
                    break#out of the iXiter loop
                end
                iXiter==maxXiter && muscadeerror(@sprintf("no convergence of step %3d after %3d iterations |Δx|=%g / %g, |Lλ|=%g / %g",step,iXiter,√(ΔX²),maxΔx,√(Lλ²)^2,maxLλ))
            end

            assemble!{:sensitivity}(outXA,asmXA,dis,model,state,Δt,(dbg...,solver=:SweepXA,step=step,iAiter=iAiter))
            ΔXₐ           = Lλx\outXA.Lλa 
            LaxΔXₐ        = outXA.Lax  ∘₁ ΔXₐ        # aa 
            ΔXₐ′Lxx       = ΔXₐ'       ∘₁ outXA.Lxx  # ax 
            La♯         .+= outXA.La  .+ ΔXₐ′Lxx ∘₁ ΔX  .- outXA.Lax ∘₁ ΔX  .- outXA.Lx ∘₁ ΔXₐ  
            Laa♯        .+= outXA.Laa .+ ΔXₐ′Lxx ∘₁ ΔXₐ .- LaxΔXₐ           .- LaxΔXₐ'             
        end # istep

        # update A
        ΔA               = Laa♯\La♯
        ΔA²,La²          = sum(ΔA.^2),sum(La♯.^2)
        decrement!(state,1,ΔA,Adofgr) 
        verbose && @printf "    In A-iteration %3d, |ΔA|=%7.1e |La♯|=%7.1e\n" iAiter √(ΔA²) √(La²)

        # Aiter convergence
        if ΔA²≤cΔA² && La²≤cLa² 
            verbose && @printf "    SweepXA converged in %3d A-iterations. |ΔA|=%7.1e / %g |La|=%7.1e / %g\n" iAiter √(ΔA²) maxΔa √(La²) maxLa
            break#out of the iAiter loop
        end
        iAiter==maxAiter && muscadeerror(@sprintf("no convergence of SweepXA after %3d A-iterations |ΔA|=%g / %g, |La|=%g / %g",iAiter,√(ΔA²),maxΔa,√(La²)^2,maxLa))

        # reset state to initial conditions
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