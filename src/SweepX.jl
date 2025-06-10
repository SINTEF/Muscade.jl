### Assembler

mutable struct AssemblySweepX{ORDER,Tλ,Tλx} <: Assembly
    # up
    Lλ        :: Tλ                
    Lλx       :: Tλx
    ming      :: 𝕣
    minλ      :: 𝕣
    Σλg       :: 𝕣
    npos      :: 𝕫
    # down
    c         :: @NamedTuple{a₁::𝕣, a₂::𝕣, a₃::𝕣, b₁::𝕣, b₂::𝕣, b₃::𝕣}
    firstiter :: 𝕓   
    line      :: 𝕓
end   
function prepare(::Type{AssemblySweepX{ORDER}},model,dis) where{ORDER}
    Xdofgr             = allXdofs(model,dis)  # dis: the model's disassembler
    ndof               = getndof(Xdofgr)
    narray,neletyp     = 2,getneletyp(model)
    asm                = Matrix{𝕫2}(undef,narray,neletyp)  # asm[iarray,ieletyp][ieledof,iele]
    Lλ                 = asmvec!(view(asm,1,:),Xdofgr,dis) 
    Lλx                = asmmat!(view(asm,2,:),view(asm,1,:),view(asm,1,:),ndof,ndof) 
    out                = AssemblySweepX{ORDER,typeof(Lλ),typeof(Lλx)}(Lλ,Lλx,∞,∞,0.,0,(a₁=0.,a₂=0.,a₃=0.,b₁=0.,b₂=0.,b₃=0.),false,false) 
    return out,asm,Xdofgr
end
function zero!(out::AssemblySweepX) 
    zero!(out.Lλ)
    zero!(out.Lλx)
    out.ming = ∞    
    out.minλ = ∞
    out.Σλg  = 0.
    out.npos = 0    
end
@inline function lineFB!(out,FB)
    if hasfield(typeof(FB),:mode) && FB.mode==:positive
        out.ming   = min(out.ming,VALUE(FB.g))
        out.minλ   = min(out.minλ,VALUE(FB.λ))
        out.Σλg   += VALUE(FB.g)*VALUE(FB.λ)
        out.npos  += 1
    end
end
function addin!(out::AssemblySweepX{ORDER},asm,iele,scale,eleobj::E,Λ,X::NTuple{Nxder,<:SVector{Nx}},U,A,t,SP,dbg) where{ORDER,E,Nxder,Nx}
    if Nx==0; return end   
    a₁,a₂,a₃,b₁,b₂,b₃ = out.c.a₁,out.c.a₂,out.c.a₃,out.c.b₁,out.c.b₂,out.c.b₃
    if ~out.line
        if ORDER==2 && out.firstiter
            i          = SVector{Nx}(1:Nx)
            δXr        = δ{1,Nx+1,𝕣}(SVector{Nx+1}(scale.X...,1.))      
            δX         = δXr[i]        
            δr         = δXr[Nx+1]     # Newmark-β special: we need C⋅a and M⋅b
            x,x′,x″    = ∂0(X),∂1(X),∂2(X)
            a          = a₂*x′ + a₃*x″
            b          = b₂*x′ + b₃*x″
            vx         = x  +    δX
            vx′        = x′ + a₁*δX - a*δr 
            vx″        = x″ + b₁*δX - b*δr
            Lλ,FB      = getresidual(eleobj,(vx,vx′,vx″),U,A,t,SP,dbg)
            Lλ         = Lλ .* scale.X
            add_value!(out.Lλ ,asm[1],iele,Lλ             )
            add_∂!{1}( out.Lλ ,asm[1],iele,Lλ,ia=1:Nx,ida=(Nx+1,))  # rhs = R - C⋅a - M⋅b 
            add_∂!{1}( out.Lλx,asm[2],iele,Lλ,ia=1:Nx,ida=1:Nx   )
        else
            δX         = δ{1,Nx,𝕣}(scale.X)
            if     ORDER==0  Lλ,FB = getresidual(eleobj,(∂0(X)+δX,                         ),U,A,t,SP,dbg)
            elseif ORDER==1  Lλ,FB = getresidual(eleobj,(∂0(X)+δX, ∂1(X)+a₁*δX             ),U,A,t,SP,dbg)
            elseif ORDER==2  Lλ,FB = getresidual(eleobj,(∂0(X)+δX, ∂1(X)+a₁*δX, ∂2(X)+b₁*δX),U,A,t,SP,dbg)
            end
            Lλ         = Lλ .* scale.X
            add_value!(out.Lλ ,asm[1],iele,Lλ)
            add_∂!{1}( out.Lλx,asm[2],iele,Lλ)
        end
    else # if out.line
        if ORDER==2 && out.firstiter
            δr         = δ{1}()              # Newmark-β special: we need C⋅a and M⋅b
            x,x′,x″    = ∂0(X),∂1(X),∂2(X)
            a          = a₂*x′ + a₃*x″
            b          = b₂*x′ + b₃*x″
            vx         = x 
            vx′        = x′ - a .*δr 
            vx″        = x″ - b .*δr 
            Lλ,FB      = getresidual(eleobj,promote(vx,vx′,vx″),U,A,t,SP,dbg)
            Lλ         = Lλ .* scale.X
            add_value!(out.Lλ ,asm[1],iele,Lλ)
            add_∂!{1}( out.Lλ ,asm[1],iele,Lλ)  # rhs = R - C⋅a - M⋅b 
            lineFB!(out,FB)
        else         
            Lλ,FB      = getresidual(eleobj,X,U,A,t,SP,dbg)
            Lλ         = Lλ .* scale.X
            add_value!(out.Lλ ,asm[1],iele,Lλ)
            lineFB!(out,FB)
        end
    end
end


"""
	SweepX{ORDER}

A non-linear, time domain solver, that solves the problem time-step by time-step.
Only the `X`-dofs of the model are solved for, while `U`-dofs and `A`-dofs are unchanged.

- `SweepX{0}` is Newton-Raphson, with feasibility line-search, to handle inequality constraints. 
- `SweepX{1}` is implicit Euler, with feasibility line-search. 
- `SweepX{2}` is Newmark-β, with Newton-Raphson iterations and feasibility line search

IMPORTANT NOTE: Muscade does not allow elements to have state variables, for example, plastic strain,
or shear-free position for dry friction.  Where the element implements such physics, this 
is implemented by introducing the state as a degree of freedom of the element, and solving
for its evolution, *even in a static problem*, requires the use of `ORDER≥1`

An analysis is carried out by a call with the following syntax:

```
initialstate    = initialize!(model)
setdof!(initialstate,1.;class=:U,field=:λcsr)
states           = solve(SweepX{2};initialstate=initialstate,time=0:10)
```
# Named arguments to `solve`:
- `dbg=(;)`           a named tuple to trace the call tree (for debugging)
- `verbose=true`      set to false to suppress printed output (for testing)
- `silenterror=false` set to true to suppress print out of error (for testing) 
- `initialstate`      a `State`, obtain from `ìnitialize!` or `SweepX`.
- `time`              maximum number of Newton-Raphson iterations 
- `β=1/4`,`γ=1/2`     parameters to the Newmark-β algorithm. Dummy if `ORDER<2`
- `maxiter=50`        maximum number of equilibrium iterations at each step.
- `maxΔx=1e-5`        convergence criteria: norm of `X`. D
- `maxLλ=∞`           convergence criteria: norm of the residual. 
- `saveiter=false`    set to true so that output `states` contains the state
                      at the iteration of the last step analysed.  Useful to study
                      a step that fails to converge. 
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

See also: [`solve`](@ref), [`initialize!`](@ref), [`findlastassigned`](@ref), [`studysingular`](@ref), [`DirectXUA`](@ref), [`FreqXU`](@ref)  
"""
struct        SweepX{ORDER} <: AbstractSolver end
function solve(SX::Type{SweepX{ORDER}},pstate,verbose,dbg;
                    time::AbstractVector{𝕣},
                    initialstate::State,
                    β::𝕣=1/4,γ::𝕣=1/2,
                    maxiter::ℤ=50,maxΔx::ℝ=1e-5,maxLλ::ℝ=∞,
                    saveiter::𝔹=false,
                    maxLineIter::ℤ=50,sfac::𝕣=.5,γfac::𝕣=.5) where{ORDER}
    model,dis        = initialstate.model,initialstate.dis
    out,asm,Xdofgr   = prepare(AssemblySweepX{ORDER},model,dis)  
    ndof             = getndof(Xdofgr)
    if ORDER≥1    x′ = 𝕣1(undef,ndof) end 
    if ORDER≥2    x″ = 𝕣1(undef,ndof) end 
    citer            = 0
    cΔx²,cLλ²        = maxΔx^2,maxLλ^2
    state            = State{1,ORDER+1,1}(copy(initialstate,SP=(γ=0.,))) 


    states           = allocate(pstate,Vector{typeof(state)}(undef,saveiter ? maxiter : length(time))) # states is not a return argument of this function.  Hence it is not lost in case of exception
    local facLλx 
    for (step,t)     ∈ enumerate(time)
        oldt         = state.time
        state.time   = t
        Δt           = t-oldt
        Δt ≤ 0 && ORDER>0 && muscadeerror(@sprintf("Time step length not strictly positive at step=%3d",step))
        if     ORDER==0 out.c= (a₁=0.      , a₂=0. , a₃=0.         , b₁=0.        , b₂=0.      , b₃=0.  )
        elseif ORDER==1 out.c= (a₁=1/Δt    , a₂=0  , a₃=0.         , b₁=0.        , b₂=0.      , b₃=0.  )
        elseif ORDER==2 out.c= (a₁=γ/(β*Δt), a₂=γ/β, a₃=(γ/2β-1)*Δt, b₁=1/(β*Δt^2), b₂=1/(β*Δt), b₃=1/2β) # γ, as in Newmark's β and γ
        end
        state.time   = t
        out.firstiter= true
        out.line     = true
        assemble!(out,asm,dis,model,state,(dbg...,solver=:SweepX,phase=:preliminary,step=step))
        out.ming ≤ 0 && muscadeerror(@sprintf("Initial point is not strictly primal-feasible at step=%3d",step)) # This is going to suck
        out.minλ ≤ 0 && muscadeerror(@sprintf("Initial point is not strictly dual-feasible at step=%3d"  ,step)) # This is going to suck
        state.SP     = (γ=out.Σλg/out.npos * γfac,)   # γ, is in interior point, g(X)*λ=γ
        for iiter    = 1:maxiter
            citer   += 1
            out.firstiter = firstiter = iiter==1
            out.line = false
            assemble!(out,asm,dis,model,state,(dbg...,solver=:SweepX,step=step,iiter=iiter))
            try if step==1 && firstiter  facLλx = lu(out.Lλx) 
            else                         lu!(facLλx, out.Lλx) 
            end catch; muscadeerror(@sprintf("matrix factorization failed at step=%i, iiter=%i",step,iiter)) end
            Δx       = facLλx\out.Lλ
            Δx²,Lλ²  = sum(Δx.^2),sum(out.Lλ.^2)
            if     ORDER==0  decr0!(state,Δx ,Xdofgr                      )
            elseif ORDER==1  decr1!(state,Δx ,Xdofgr,out.c                )
            elseif ORDER==2  decr2!(state,Δx ,Xdofgr,out.c,firstiter,x′,x″)
            end

            out.line = true    
            s = 1.    
            for iline = 1:maxLineIter
                assemble!(out,asm,dis,model,state,(dbg...,solver=:SweepX,phase=:linesearch,step=step,iiter=iiter,iline=iline))
                out.minλ > 0 && out.ming > 0 &&  break
                iline==maxLineIter && muscadeerror(@sprintf("Line search failed at step=%3d, iiter=%3d, iline=%3d, s=%7.1e",step,iiter,iline,s))
                Δs    = s*(sfac-1)
                s    += Δs
                if     ORDER==0  decr0!(state,Δs*Δx ,Xdofgr                      )
                elseif ORDER==1  decr1!(state,Δs*Δx ,Xdofgr,out.c                )
                elseif ORDER==2  decr2!(state,Δs*Δx ,Xdofgr,out.c,firstiter,x′,x″)
                end
            end

            verbose && saveiter && @printf("        iteration %3d, γ= %7.1e\n",iiter,γ)
            saveiter && (states[iiter]=State(state.time,state.Λ,deepcopy(state.X),state.U,state.A,state.SP,model,dis))
            if Δx²≤cΔx² && Lλ²≤cLλ² 
                verbose && @printf "    step %3d converged in %3d iterations. |Δx|=%7.1e |Lλ|=%7.1e\n" step iiter √(Δx²) √(Lλ²)
                ~saveiter && (states[step]=State(state.time,state.Λ,deepcopy(state.X),state.U,state.A,state.SP,model,dis))
                break#out of the iiter loop
            end
            iiter==maxiter && muscadeerror(@sprintf("no convergence of step %3d after %3d iterations |Δx|=%g / %g, |Lλ|=%g / %g",step,iiter,√(Δx²),maxΔx,√(Lλ²)^2,maxLλ))
            state.SP     = (γ=state.SP.γ*γfac,)
        end
    end
    verbose && @printf "\n    nel=%d, ndof=%d, nstep=%d, niter=%d, niter/nstep=%5.2f\n" getnele(model) getndof(Xdofgr) length(time) citer citer/length(time)
    return
end
function decr2!(state,Δx ,Xdofgr,c,firstiter,x′,x″)
    a₁,a₂,a₃,b₁,b₂,b₃ = c.a₁,c.a₂,c.a₃,c.b₁,c.b₂,c.b₃
    if firstiter
        getdof!(state,1,x′,Xdofgr) 
        getdof!(state,2,x″,Xdofgr) 
        a        = a₂*x′+a₃*x″
        b        = b₂*x′+b₃*x″
        Δx′      = a₁*Δx + a
        Δx″      = b₁*Δx + b
    else
        Δx′      = a₁*Δx 
        Δx″      = b₁*Δx 
    end
    decrement!(state,1,Δx ,Xdofgr)
    decrement!(state,2,Δx′,Xdofgr)
    decrement!(state,3,Δx″,Xdofgr)
end
function decr1!(state,Δx ,Xdofgr,c)
    Δx′      = c.a₁*Δx 
    decrement!(state,1,Δx ,Xdofgr)
    decrement!(state,2,Δx′,Xdofgr)
end
function decr0!(state,Δx ,Xdofgr)
    decrement!(state,1,Δx ,Xdofgr)
end
