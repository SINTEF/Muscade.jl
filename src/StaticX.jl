
###--------------------- ASMstaticX: for good old static FEM

mutable struct AssemblyStaticX{Tλ,Tλx} <:Assembly
    Lλ    :: Tλ
    Lλx   :: Tλx 
    α     :: 𝕣
end   
function prepare(::Type{AssemblyStaticX},model,dis) 
    dofgr              = allXdofs(model,dis)
    ndof               = getndof(dofgr)
    narray,neletyp     = 2,getneletyp(model)
    asm                = Matrix{𝕫2}(undef,narray,neletyp)  
    Lλ                 = asmvec!(view(asm,1,:),dofgr,dis) 
    Lλx                = asmmat!(view(asm,2,:),view(asm,1,:),view(asm,1,:),ndof,ndof) 
#    out                = one_for_each_thread(AssemblyStaticX(Lλ,Lλx,∞)) # KEEP - parallel
    out                = AssemblyStaticX(Lλ,Lλx,∞) # sequential
    return out,asm,dofgr
end
function zero!(out::AssemblyStaticX)
    zero!(out.Lλ)
    zero!(out.Lλx)
    out.α = ∞    
end
function add!(out1::AssemblyStaticX,out2::AssemblyStaticX) 
    add!(out1.Lλ,out2.Lλ)
    add!(out1.Lλx,out2.Lλx)
    out1.α = min(out1.α,out2.α)
end
function addin!(out::AssemblyStaticX,asm,iele,scale,eleobj::E,Λ,X::NTuple{Nxder,<:SVector{Nx}},U,A,t,SP,dbg) where{E,Nxder,Nx}
    if Nx==0; return end # don't waste time on Acost elements...  
    ΔX         = δ{1,Nx,𝕣}(scale.X)                 # NB: precedence==1, input must not be Adiff 
    Lλ,χ,FB    = getresidual(implemented(eleobj)...,eleobj,(∂0(X)+ΔX,),U,A,t,nothing,nothing,SP,dbg)
    Lλ         = Lλ .* scale.X
    add_value!(out.Lλ ,asm[1],iele,Lλ)
    add_∂!{1}( out.Lλx,asm[2],iele,Lλ)
    out.α      = min(out.α,default{:α}(FB,∞))
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
- `γ0=1.`             an initial value of the barrier coefficient for the handling of contact
                      using an interior point method
- `γfac1=0.5`         at each iteration, the barrier parameter γ is multiplied 
- `γfac2=100.`        by γfac1*exp(-min(αᵢ)/γfac2)^2), where αᵢ is computed by the i-th
                      interior point savvy element as αᵢ=abs(λ-g)/γ                                               

# Output
A vector of length equal to that of `time` containing the state of the model at each of these steps                       

See also: [`solve`](@ref), [`StaticXUA`](@ref), [`initialize!`](@ref)
"""
struct StaticX <: AbstractSolver end 
function solve(::Type{StaticX},pstate,verbose,dbg;
                    time::AbstractVector{𝕣},
                    initialstate::State,
                    maxiter::ℤ=50,maxΔx::ℝ=1e-5,maxresidual::ℝ=∞,
                    saveiter::𝔹=false,γ0::𝕣=1.,γfac1::𝕣=.5,γfac2::𝕣=100.)
    # important: this code assumes that there is no χ in state.
    model,dis        = initialstate.model,initialstate.dis
    out,asm,dofgr    = prepare(AssemblyStaticX,model,dis)
    citer            = 0
    cΔx²,cLλ²        = maxΔx^2,maxresidual^2
    s                = State{1,1}(initialstate,(γ=0.,))
    state            = allocate(pstate,Vector{typeof(s)}(undef,saveiter ? maxiter : length(time))) # state is not a return argument of this function.  Hence it is not lost in case of exception
    local facLλx 
    for (step,t)     ∈ enumerate(time)
        s.time       = t
        s.SP         = (γ=γ0,)
        for iiter    = 1:maxiter
            citer   += 1
            assemble!(out,asm,dis,model,s,(dbg...,solver=:StaticX,step=step,iiter=iiter))

            try if step==1 && iiter==1
                facLλx = lu(firstelement(out).Lλx) 
            else
                lu!(facLλx,firstelement(out).Lλx) 
            end catch; muscadeerror(@sprintf("matrix factorization failed at step=%i, iiter=%i",step,iiter)) end
            Δx       = facLλx\firstelement(out).Lλ
            Δx²,Lλ²  = sum(Δx.^2),sum(firstelement(out).Lλ.^2)
            decrement!(s,0,Δx,dofgr)
            verbose && saveiter && @printf("        iteration %3d, γ= %7.1e\n",iiter,γ)
            saveiter && (state[iiter]=State(s.Λ,deepcopy(s.X),s.U,s.A,s.time,s.SP,model,dis))
            if Δx²≤cΔx² && Lλ²≤cLλ² 
                verbose && @printf "    step %3d converged in %3d iterations. |Δx|=%7.1e |Lλ|=%7.1e\n" step iiter √(Δx²) √(Lλ²)
                ~saveiter && (state[step]=State(s.Λ,deepcopy(s.X),s.U,s.A,s.time,s.SP,model,dis))
                break#out of the iiter loop
            end
            iiter==maxiter && muscadeerror(@sprintf("no convergence in step %3d after %3d iterations |Δx|=%g / %g, |Lλ|=%g / %g",step,iiter,√(Δx²),maxΔx,√(Lλ²)^2,maxresidual))
            Δγ       = γfac1*exp(-(firstelement(out).α/γfac2)^2)
            s.SP     = (γ=s.SP.γ*Δγ,)
        end
    end
    verbose && @printf "\n    nel=%d, ndof=%d, nstep=%d, niter=%d, niter/nstep=%5.2f\n" getnele(model) getndof(dofgr) length(time) citer citer/length(time)
    return
end
