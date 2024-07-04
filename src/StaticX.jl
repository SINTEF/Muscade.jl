
###--------------------- ASMstaticX: for good old static FEM

mutable struct AssemblyStaticX{Tλ,Tλx} <:Assembly
    Lλ    :: Tλ
    Lλx   :: Tλx 
end   
function prepare(::Type{AssemblyStaticX},model,dis) 
    dofgr              = allXdofs(model,dis)
    ndof               = getndof(dofgr)
    narray,neletyp     = 2,getneletyp(model)
    asm                = Matrix{𝕫2}(undef,narray,neletyp)  
    Lλ                 = asmvec!(view(asm,1,:),dofgr,dis) 
    Lλx                = asmmat!(view(asm,2,:),view(asm,1,:),view(asm,1,:),ndof,ndof) 
#    out                = one_for_each_thread(AssemblyStaticX(Lλ,Lλx,∞)) # KEEP - parallel
    out                = AssemblyStaticX(Lλ,Lλx) # sequential
    return out,asm,dofgr
end
function zero!(out::AssemblyStaticX)
    zero!(out.Lλ)
    zero!(out.Lλx)
end
# function add!(out1::AssemblyStaticX,out2::AssemblyStaticX) 
#     add!(out1.Lλ,out2.Lλ)
#     add!(out1.Lλx,out2.Lλx)
# end
function addin!(out::AssemblyStaticX,asm,iele,scale,eleobj::E,Λ,X::NTuple{Nxder,<:SVector{Nx}},U,A,t,SP,dbg) where{E,Nxder,Nx}
    if Nx==0; return end # don't waste time on Acost elements...  
    ΔX         = δ{1,Nx,𝕣}(scale.X)                 # NB: precedence==1, input must not be Adiff 
    Lλ,FB    = getresidual(eleobj,(∂0(X)+ΔX,),U,A,t,SP,dbg) #  no feedback FB
    Lλ         = Lλ .* scale.X
    add_value!(out.Lλ ,asm[1],iele,Lλ)
    add_∂!{1}( out.Lλx,asm[2],iele,Lλ)
end
###--------------------- ASMstaticXline: for line search

mutable struct AssemblyStaticXline{Tλ} <:Assembly
    Lλ    :: Tλ
    ming  :: 𝕣
    minλ  :: 𝕣
    Σλg   :: 𝕣
    npos  :: 𝕫
end   
function prepare(::Type{AssemblyStaticXline},model,dis) 
    dofgr              = allXdofs(model,dis)
    narray,neletyp     = 1,getneletyp(model)
    asm                = Matrix{𝕫2}(undef,narray,neletyp)  
    Lλ                 = asmvec!(view(asm,1,:),dofgr,dis) 
    out                = AssemblyStaticXline(Lλ,∞,∞,0.,0) 
    return out,asm
end
function zero!(out::AssemblyStaticXline)
    zero!(out.Lλ)
    out.ming = ∞    
    out.minλ = ∞
    out.Σλg  = 0.
    out.npos = 0    
end
# function add!(out1::AssemblyStaticXline,out2::AssemblyStaticXline) 
#     add!(out1.Lλ,out2.Lλ)
#     out1.ming = min(out1.ming,out2.ming)
#     out1.minλ = min(out1.minλ,out2.minλ)
#     out1.Σλg += out2.Σλg
#     out1.npos+= out2.npos
# end
function addin!(out::AssemblyStaticXline,asm,iele,scale,eleobj::E,Λ,X::NTuple{Nxder,<:SVector{Nx}},U,A,t,SP,dbg) where{E,Nxder,Nx}
    if Nx==0; return end # don't waste time on Acost elements...  
    Lλ,FB = getresidual(eleobj,X,U,A,t,SP,dbg)
    Lλ      = Lλ .* scale.X
    add_value!(out.Lλ ,asm[1],iele,Lλ) 
    if hasfield(typeof(FB),:mode) && FB.mode==:positive
        out.ming   = min(out.ming,VALUE(FB.g))
        out.minλ   = min(out.minλ,VALUE(FB.λ))
        out.Σλg   += VALUE(FB.g)*VALUE(FB.λ)
        out.npos  += 1
    end
end

"""
    StaticX

A non-linear static solver for forward (not inverse, optimisation) FEM.
The current implementation does not handle element memory. 

An analysis is carried out by a call with the following syntax:

```
initialstate    = initialize!(model)
state           = solve(StaticX;initialstate=initialstate,time=[0.,1.])
```

# Named arguments
- `dbg=(;)`           a named tuple to trace the call tree (for debugging)
- `verbose=true`      set to false to suppress printed output (for testing)
- `silenterror=false` set to true to suppress print out of error (for testing) 
- `initialstate`      a single `state` - obtained from a call to `initialize!`, or 
                      from a previous analysis
- `time`              an `AbstractVector` vector of the times at which to compute 
                      equilibrium.  While this solver does not account for dynamic
                      effect, the model will typicaly describe some loads as time 
                      dependent. 
- `maxiter=50`        maximum number of Newton-Raphson iteration at any given step 
- `maxΔx=1e-5`        convergence criteria: a norm on the scaled `X` increment 
- `maxincrement=∞`    convergence criteria: a norm on the scaled residual
- `saveiter=false`    set to true so that the output `state` is a vector describing 
                      the states of the model at the last iteration (for debugging 
                      non-convergence) 
- `γfac=0.5`          at each iteration, the barrier parameter γ is multiplied 

# Output
A vector of length equal to that of `time` containing the state of the model at each of these steps                       

See also: [`solve`](@ref), [`StaticXUA`](@ref), [`initialize!`](@ref)
"""
struct StaticX <: AbstractSolver end 
function solve(::Type{StaticX},pstate,verbose,dbg;
                    time::AbstractVector{𝕣},
                    initialstate::State,
                    maxiter::ℤ=50,maxΔx::ℝ=1e-5,maxresidual::ℝ=∞,
                    saveiter::𝔹=false,
                    maxLineIter::ℤ=50,α::𝕣=.1,β::𝕣=.5,γfac::𝕣=.5)
    model,dis        = initialstate.model,initialstate.dis
    out1,asm1,Xdofgr = prepare(AssemblyStaticX    ,model,dis)
    out2,asm2        = prepare(AssemblyStaticXline,model,dis)
    citer            = 0
    cΔx²,cLλ²        = maxΔx^2,maxresidual^2
    state            = State{1,1,1}(initialstate,(γ=0.,))
    states           = allocate(pstate,Vector{typeof(state)}(undef,saveiter ? maxiter : length(time))) # state is not a return argument of this function.  Hence it is not lost in case of exception
    local facLλx 
    for (step,t)     ∈ enumerate(time)
        state.time   = t
        assemble!(out2,asm2,dis,model,state,(dbg...,solver=:StaticX,phase=:preliminary,step=step))
        out2.ming ≤ 0 && muscadeerror(@sprintf("Initial point is not strictly primal-feasible at step=%3d",step))
        out2.minλ ≤ 0 && muscadeerror(@sprintf("Initial point is not strictly dual-feasible at step=%3d"  ,step))
        state.SP     = (γ=out2.Σλg/out2.npos * γfac,)
        for iiter    = 1:maxiter
            citer   += 1
            assemble!(out1,asm1,dis,model,state,(dbg...,solver=:StaticX,phase=:direction,step=step,iiter=iiter))
            try if step==1 && iiter==1
                facLλx = lu(firstelement(out1).Lλx) 
            else
                lu!(facLλx,firstelement(out1).Lλx) 
            end catch; muscadeerror(@sprintf("matrix factorization failed at step=%3d, iiter=%3d",step,iiter)) end
            Δx       = facLλx\firstelement(out1).Lλ
            Δx²,Lλ²  = sum(Δx.^2),sum(firstelement(out1).Lλ.^2)
            decrement!(state,0,Δx,Xdofgr)

            s = 1.    
            for iline = 1:maxLineIter
                assemble!(out2,asm2,dis,model,state,(dbg...,solver=:StaticX,phase=:linesearch,step=step,iiter=iiter,iline=iline))
                out2.minλ > 0 && out2.ming > 0 && sum(firstelement(out2).Lλ.^2) ≤ Lλ²*(1-α*s)^2 && break
                iline==maxLineIter && muscadeerror(@sprintf("Line search failed at step=%3d, iiter=%3d, iline=%3d, s=%7.1e",step,iiter,iline,s))
                Δs = s*(β-1)
                s += Δs
                decrement!(state,0,Δs*Δx,Xdofgr)
            end

            verbose && saveiter && @printf("        iteration %3d, γ= %7.1e\n",iiter,γ)
            saveiter && (states[iiter]=State(state.time,state.Λ,deepcopy(state.X),state.U,state.A,state.SP,model,dis))
            if Δx²*s^2≤cΔx² && Lλ²≤cLλ² 
                verbose && @printf "    step %3d converged in %3d iterations. |Δx|=%7.1e |Lλ|=%7.1e\n" step iiter √(Δx²) √(Lλ²)
                ~saveiter && (states[step]=State(state.time,state.Λ,deepcopy(state.X),state.U,state.A,state.SP,model,dis))
                break#out of the iiter loop
            end
            iiter==maxiter && muscadeerror(@sprintf("no convergence at step=%3d, iiter=%3d, |Δx|=%7.1e / %7.1e, |Lλ|=%7.1e / %7.1e",step,iiter,√(Δx²),maxΔx,√(Lλ²)^2,maxresidual))
            state.SP     = (γ=state.SP.γ*γfac,)
        end
    end
    verbose && @printf "\n    nel=%d, ndof=%d, nstep=%d, niter=%d, niter/nstep=%5.2f\n" getnele(model) getndof(Xdofgr) length(time) citer citer/length(time)
    return
end
