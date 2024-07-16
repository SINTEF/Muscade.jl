# standard Newmark-β solver

mutable struct AssemblyNewmarkX{Tλ,Tλx} <:Assembly
    # from model to solver
    Lλ      :: Tλ                # Newmark-β rhs
    Lλx     :: Tλx               # Newmark-β incremental matrix  
#    α       :: 𝕣                 # feedback to solver for interior point
    # from solver to assembler
    a₁      :: 𝕣                 # coefficients for linear combinations in Newmark-β
    a₂      :: 𝕣                 # ... 
    a₃      :: 𝕣
    b₁      :: 𝕣
    b₂      :: 𝕣
    b₃      :: 𝕣
end   
function prepare(::Type{AssemblyNewmarkX},model,dis) 
    Xdofgr             = allXdofs(model,dis)  # dis: the model's disassembler
    ndof               = getndof(Xdofgr)
    narray,neletyp     = 2,getneletyp(model)
    asm                = Matrix{𝕫2}(undef,narray,neletyp)  # asm[iarray,ieletyp][ieledof,iele]
    Lλ                 = asmvec!(view(asm,1,:),Xdofgr,dis) 
    Lλx                = asmmat!(view(asm,2,:),view(asm,1,:),view(asm,1,:),ndof,ndof) 
#    out                = AssemblyNewmarkX(Lλ,Lλx,∞,0.,0.,0.,0.,0.,0.) # sequential
    out                = AssemblyNewmarkX(Lλ,Lλx,0.,0.,0.,0.,0.,0.) # sequential
    return out,asm,Xdofgr
end
function zero!(out::AssemblyNewmarkX)
    zero!(out.Lλ)
    zero!(out.Lλx)
#    out.α = ∞    
end
function addin!(out::AssemblyNewmarkX,asm,iele,scale,eleobj::E,Λ,X::NTuple{Nxder,<:SVector{Nx}},U,A,t,SP,dbg) where{E,Nxder,Nx}
    # asm[iarray][ientry,iel]
    if Nx==0; return end # don't waste time on Acost elements...  
    i          = SVector{Nx}(1:Nx)
    δkr        = δ{1,Nx+1,𝕣}(SVector{Nx+1}(scale.X...,1.))      
    δk         = δkr[i]        # classic
    δr         = δkr[Nx+1]     # Newmark-β special: we need Lλx′⋅a and Lλx″⋅b
    x,x′,x″    = ∂0(X),∂1(X),∂2(X)
    a          = out.a₂*x′ + out.a₃*x″
    b          = out.b₂*x′ + out.b₃*x″
    vx         = x  +        δk
    vx′        = x′ + out.a₁*δk - a*δr 
    vx″        = x″ + out.b₁*δk - b*δr 
    Lλ,FB      = getresidual(eleobj,(vx,vx′,vx″),U,A,t,SP,dbg)
    Lλ         = Lλ .* scale.X
    add_value!(out.Lλ ,asm[1],iele,Lλ             )
    add_∂!{1}( out.Lλ ,asm[1],iele,Lλ,1:Nx,(Nx+1,))  # rhs = Lλ - Lλx′⋅a - Lλx″⋅b 
    add_∂!{1}( out.Lλx,asm[2],iele,Lλ,1:Nx,1:Nx   )
end

###

mutable struct AssemblyNewmarkXline{Tλ} <:Assembly
    Lλ    :: Tλ
    ming  :: 𝕣
    minλ  :: 𝕣
    Σλg   :: 𝕣
    npos  :: 𝕫
    # from solver to assembler
    a₁      :: 𝕣                 # coefficients for linear combinations in Newmark-β
    a₂      :: 𝕣                 # ... 
    a₃      :: 𝕣
    b₁      :: 𝕣
    b₂      :: 𝕣
    b₃      :: 𝕣
end   
function prepare(::Type{AssemblyNewmarkXline},model,dis) 
    dofgr              = allXdofs(model,dis)
    narray,neletyp     = 1,getneletyp(model)
    asm                = Matrix{𝕫2}(undef,narray,neletyp)  
    Lλ                 = asmvec!(view(asm,1,:),dofgr,dis) 
    out                = AssemblyNewmarkXline(Lλ,∞,∞,0.,0, 0.,0.,0.,0.,0.,0.) 
    return out,asm
end
function zero!(out::AssemblyNewmarkXline)
    zero!(out.Lλ)
    out.ming = ∞    
    out.minλ = ∞
    out.Σλg  = 0.
    out.npos = 0    
end
function addin!(out::AssemblyNewmarkXline,asm,iele,scale,eleobj::E,Λ,X::NTuple{Nxder,<:SVector{Nx}},U,A,t,SP,dbg) where{E,Nxder,Nx}
    if Nx==0; return end # don't waste time on Acost elements...  
    δr         = δ{1}()              # Newmark-β special: we need Lλx′⋅a and Lλx″⋅b
    x,x′,x″    = ∂0(X),∂1(X),∂2(X)
    a          = out.a₂*x′ + out.a₃*x″
    b          = out.b₂*x′ + out.b₃*x″
    vx         = x  .+ 0*δr
    vx′        = x′ - a.*δr 
    vx″        = x″ - b.*δr 
    Lλ,FB      = getresidual(eleobj,(vx,vx′,vx″),U,A,t,SP,dbg)
    Lλ         = Lλ .* scale.X
    add_value!(out.Lλ ,asm[1],iele,Lλ)
    add_∂!{1}( out.Lλ ,asm[1],iele,Lλ)  # rhs = Lλ - Lλx′⋅a - Lλx″⋅b 
end


"""
	NewmarkX

A non-linear dynamic time domain solver.
The algorutm is Newmark-β
The current implementation does not handle element memory. 

An analysis is carried out by a call with the following syntax:

```
initialstate    = initialize!(model)
setdof!(initialstate,1.;class=:U,field=:λcsr)
states           = solve(NewmarkX  ;initialstate=initialstate,time=0:10)
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
struct NewmarkX <: AbstractSolver end
function solve(::Type{NewmarkX},pstate,verbose,dbg;
                    time::AbstractVector{𝕣},
                    initialstate::State,
                    β::ℝ=1/4,γ::ℝ=1/2,
                    maxiter::ℤ=50,maxΔx::ℝ=1e-5,maxresidual::ℝ=∞,
                    saveiter::𝔹=false,
                    maxLineIter::ℤ=50,line1::𝕣=.1,line2::𝕣=.5,γfac::𝕣=.5)
    model,dis        = initialstate.model,initialstate.dis
    out1,asm1,Xdofgr = prepare(AssemblyNewmarkX    ,model,dis)
    out2,asm2        = prepare(AssemblyNewmarkXline,model,dis)
    ndof             = getndof(Xdofgr)
    x′ ,x″           = 𝕣1(undef,ndof), 𝕣1(undef,ndof) 
    citer            = 0
    cΔx²,cLλ²        = maxΔx^2,maxresidual^2
    state            = State{1,3,1}(initialstate,(γ=0.,))
    states           = allocate(pstate,Vector{typeof(state)}(undef,saveiter ? maxiter : length(time))) # states is not a return argument of this function.  Hence it is not lost in case of exception
    local facLλx 
    for (step,t)     ∈ enumerate(time)
        oldt         = state.time
        state.time   = t
        Δt           = t-oldt
        for iiter    = 1:maxiter
            if iiter == 1
                out1.a₁,out1.a₂,out1.a₃ = γ/(β*Δt),   γ/β,      (γ/2β-1)*Δt
                out1.b₁,out1.b₂,out1.b₃ = 1/(β*Δt^2), 1/(β*Δt), 1/2β
                out2.a₁,out2.a₂,out2.a₃ = γ/(β*Δt),   γ/β,      (γ/2β-1)*Δt
                out2.b₁,out2.b₂,out2.b₃ = 1/(β*Δt^2), 1/(β*Δt), 1/2β
            else
                out1.a₂,out1.a₃        = 0., 0.
                out1.b₂,out1.b₃        = 0., 0.
                out2.a₂,out2.a₃        = 0., 0.
                out2.b₂,out2.b₃        = 0., 0.
            end
            citer   += 1
            assemble!(out1,asm1,dis,model,state,(dbg...,solver=:NewmarkX,step=step,iiter=iiter))
            try if step==1 && iiter==1
                facLλx = lu(out1.Lλx) 
            else
                lu!(facLλx,out1.Lλx) 
            end catch; muscadeerror(@sprintf("matrix factorization failed at step=%i, iiter=%i",step,iiter)) end
            Δx       = facLλx\out1.Lλ
            Δx²,Lλ²  = sum(Δx.^2),sum(out1.Lλ.^2)
            getdof!(state,1,x′,Xdofgr) 
            getdof!(state,2,x″,Xdofgr) 
            a        = out1.a₂*x′+out1.a₃*x″
            b        = out1.b₂*x′+out1.b₃*x″
            Δx′      = out1.a₁*Δx + a
            Δx″      = out1.b₁*Δx + b
            decrement!(state,0,Δx ,Xdofgr)
            decrement!(state,1,Δx′,Xdofgr)
            decrement!(state,2,Δx″,Xdofgr)

            s = 1.    
            for iline = 1:maxLineIter
                assemble!(out2,asm2,dis,model,state,(dbg...,solver=:NewmarkX,phase=:linesearch,step=step,iiter=iiter,iline=iline))
                # TODO
                #
                # The requirement       sum(out2.Lλ.^2) ≤ Lλ²*(1-line1*s)^2      leads to failure, also in the absence of any contraint. 
                # 1) Does this imply a bug in update "decrement(...δx...)
                # 2) Is that because requiring Lλ² to always decrease is a bad idea for non-convex problems?
                #
                # @show step,iiter,iline, s, sum(out2.Lλ.^2), Lλ²*(1-line1*s)^2
                # out2.minλ > 0 && out2.ming > 0 && sum(out2.Lλ.^2) ≤ Lλ²*(1-line1*s)^2 && break
                out2.minλ > 0 && out2.ming > 0 &&  break
                iline==maxLineIter && muscadeerror(@sprintf("Line search failed at step=%3d, iiter=%3d, iline=%3d, s=%7.1e",step,iiter,iline,s))
                Δs    = s*(line2-1)
                s    += Δs
                δx    = Δs*Δx
                δx′   = out1.a₁*δx + a
                δx″   = out1.b₁*δx + b
                decrement!(state,0,δx ,Xdofgr)
                decrement!(state,1,δx′,Xdofgr)  
                decrement!(state,2,δx″,Xdofgr)
            end

            verbose && saveiter && @printf("        iteration %3d, γ= %7.1e\n",iiter,γ)
            saveiter && (states[iiter]=State(state.time,state.Λ,deepcopy(state.X),state.U,state.A,state.SP,model,dis))
            if Δx²≤cΔx² && Lλ²≤cLλ² 
                verbose && @printf "    step %3d converged in %3d iterations. |Δx|=%7.1e |Lλ|=%7.1e\n" step iiter √(Δx²) √(Lλ²)
                ~saveiter && (states[step]=State(state.time,state.Λ,deepcopy(state.X),state.U,state.A,state.SP,model,dis))
                break#out of the iiter loop
            end
            iiter==maxiter && muscadeerror(@sprintf("no convergence in step %3d after %3d iterations |Δx|=%g / %g, |Lλ|=%g / %g",step,iiter,√(Δx²),maxΔx,√(Lλ²)^2,maxresidual))
            state.SP     = (γ=state.SP.γ*γfac,)
        end
    end
    verbose && @printf "\n    nel=%d, ndof=%d, nstep=%d, niter=%d, niter/nstep=%5.2f\n" getnele(model) getndof(Xdofgr) length(time) citer citer/length(time)
    return
end
