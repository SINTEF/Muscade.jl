# standard Newmark-β solver

mutable struct AssemblyDynamicX{Tλ,Tλx} <:Assembly
    # from model to solver
    Lλ      :: Tλ                # Newmark-β rhs
    Lλx     :: Tλx               # Newmark-β incremental matrix  
    α       :: 𝕣                 # feedback to solver for interior point
    # from solver to assembler
    a₁      :: 𝕣                 # coefficients for linear combinations in Newmark-β
    a₂      :: 𝕣                 # ... 
    a₃      :: 𝕣
    b₁      :: 𝕣
    b₂      :: 𝕣
    b₃      :: 𝕣
end   
function prepare(::Type{AssemblyDynamicX},model,dis,β,γ) 
    dofgr              = allXdofs(model,dis)  # dis: the model's disassembler
    ndof               = getndof(dofgr)
    narray,neletyp     = 2,getneletyp(model)
    asm                = Matrix{𝕫2}(undef,narray,neletyp)  # asm[iarray,ieletyp][ieledof,iele]
    Lλ                 = asmvec!(view(asm,1,:),dofgr,dis) 
    Lλx                = asmmat!(view(asm,2,:),view(asm,1,:),view(asm,1,:),ndof,ndof) 
#    out                = one_for_each_thread(AssemblyDynamicX(Lλ,Lλx,∞)) # KEEP - parallel 
    out                = AssemblyDynamicX(Lλ,Lλx,∞,0.,0.,0.,0.,0.,0.) # sequential
    return out,asm,dofgr
end
function zero!(out::AssemblyDynamicX)
    zero!(out.Lλ)
    zero!(out.Lλx)
    out.α = ∞    
end
function add!(out1::AssemblyDynamicX,out2::AssemblyDynamicX) 
    add!(out1.Lλ,out2.Lλ)
    add!(out1.Lλx,out2.Lλx)
    out1.α = min(out1.α,out2.α)
end
function addin!(out::AssemblyDynamicX,asm,iele,scale,eleobj::E,Λ,X::NTuple{Nxder,<:SVector{Nx}},U,A,t,SP,dbg) where{E,Nxder,Nx}
    # asm[iarray][ientry,iel]
    if Nx==0; return end # don't waste time on Acost elements...  
    i          = SVector{Nx}(1:Nx)
    δkr        = δ{1,Nx+1,𝕣}(SVector{Nx+1}(scale.X...,1.))      
    δk         = δkr[i]
    δr         = δkr[Nx+1]
    x,x′,x″    = ∂0(X),∂1(X),∂2(X)
    a          = out.a₂*x′ + out.a₃*x″
    b          = out.b₂*x′ + out.b₃*x″
    vx         = x  +        δk
    vx′        = x′ + out.a₁*δk - a*δr 
    vx″        = x″ + out.b₁*δk - b*δr                                      # χo     ,χcv
#    Lλ,χ,FB    = getresidual(implemented(eleobj)...,eleobj,(vx,vx′,vx″),U,A,t,nothing,nothing,SP,dbg)
    Lλ,FB      = getresidual(eleobj,(vx,vx′,vx″),U,A,t,SP,dbg)
    add_value!(out.Lλ ,asm[1],iele,Lλ             )
    add_∂!{1}( out.Lλ ,asm[1],iele,Lλ,1:Nx,(Nx+1,))
    add_∂!{1}( out.Lλx,asm[2],iele,Lλ,1:Nx,1:Nx   )
    out.α      = min(out.α,default{:α}(FB,∞))
end

"""
	DynamicX

A non-linear dynamic time domain solver.
The algorutm is Newmark-β
The current implementation does not handle element memory. 

An analysis is carried out by a call with the following syntax:

```
initialstate    = initialize!(model)
setdof!(initialstate,1.;class=:U,field=:λcsr)
state           = solve(DynamicX  ;initialstate=initialstate,time=0:10)
```

# Named arguments to `solve`:
- `dbg=(;)`           a named tuple to trace the call tree (for debugging)
- `verbose=true`      set to false to suppress printed output (for testing)
- `silenterror=false` set to true to suppress print out of error (for testing) 
- `initialstate`      a `State`, obtain from `ìnitialize!` or `StaticX`.
- `time`              maximum number of Newton-Raphson iterations 
- `β=1/5`,`γ=1/2`     parameters to the Newmark-β algorithm.
- `maxiter=50`        maximum number of equilibrium iterations at each step.
- `maxΔx=1e-5`        convergence criteria: norm of `X`. Default.
- `maxresidual=∞`     convergence criteria: norm of the residual. 

# Output

A vector of length equal to that of the named input argument `time` containing the states at the time steps.                       

See also: [`solve`](@ref), [`StaticX`](@ref), [`setdof!`](@ref) 
"""
struct DynamicX <: AbstractSolver end
function solve(::Type{DynamicX},pstate,verbose,dbg;
                    time::AbstractVector{𝕣},
                    initialstate::State,
                    β::ℝ=1/4,γ::ℝ=1/2,
                    maxiter::ℤ=50,maxΔx::ℝ=1e-5,maxresidual::ℝ=∞,
                    saveiter::𝔹=false,γ0::𝕣=1.,γfac1::𝕣=.5,γfac2::𝕣=100.)
    # important: this code assumes that there is no χ in state.
    model,dis        = initialstate.model,initialstate.dis
    out,asm,dofgr    = prepare(AssemblyDynamicX,model,dis,β,γ)
    citer            = 0
    cΔx²,cLλ²        = maxΔx^2,maxresidual^2
    s                = State{1,3,1}(initialstate,(γ=0.,))
    state            = allocate(pstate,Vector{typeof(s)}(undef,saveiter ? maxiter : length(time))) # state is not a return argument of this function.  Hence it is not lost in case of exception
    local facLλx 
    for (step,t)     ∈ enumerate(time)
        oldt         = s.time
        s.time       = t
        Δt           = t-oldt
        s.SP         = (γ=γ0,)
        for iiter    = 1:maxiter
            if iiter == 1
                out.a₁,out.a₂,out.a₃ = γ/(β*Δt),   γ/β,      (γ/2β-1)*Δt
                out.b₁,out.b₂,out.b₃ = 1/(β*Δt^2), 1/(β*Δt), 1/2β
            else
                out.a₂,out.a₃        = 0., 0.
                out.b₂,out.b₃        = 0., 0.
            end
            citer   += 1
            assemble!(out,asm,dis,model,s,(dbg...,solver=:DynamicX,step=step,iiter=iiter))
            try if step==1 && iiter==1
                facLλx = lu(firstelement(out).Lλx) 
            else
                lu!(facLλx,firstelement(out).Lλx) 
            end catch; muscadeerror(@sprintf("matrix factorization failed at step=%i, iiter=%i",step,iiter)) end
            Δx       = facLλx\firstelement(out).Lλ
            Δx²,Lλ²  = sum(Δx.^2),sum(firstelement(out).Lλ.^2)
            x′ ,x″   = Vector{𝕣}(undef,length(Δx)), Vector{𝕣}(undef,length(Δx))
            getdof!(s,1,x′,dofgr) 
            getdof!(s,2,x″,dofgr) 
            Δx′      = out.a₁*Δx+out.a₂*x′+out.a₃*x″ 
            Δx″      = out.b₁*Δx+out.b₂*x′+out.b₃*x″
            decrement!(s,0,Δx ,dofgr)
            decrement!(s,1,Δx′,dofgr)
            decrement!(s,2,Δx″,dofgr)
            verbose && saveiter && @printf("        iteration %3d, γ= %7.1e\n",iiter,γ)
            saveiter && (state[iiter]=State(s.time,s.Λ,deepcopy(s.X),s.U,s.A,s.SP,model,dis))
            if Δx²≤cΔx² && Lλ²≤cLλ² 
                verbose && @printf "    step %3d converged in %3d iterations. |Δx|=%7.1e |Lλ|=%7.1e\n" step iiter √(Δx²) √(Lλ²)
                ~saveiter && (state[step]=State(s.time,s.Λ,deepcopy(s.X),s.U,s.A,s.SP,model,dis))
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
