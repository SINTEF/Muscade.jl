### Assembler

mutable struct AssemblySweepXA{ORDER} <: Assembly
    # up
    Lλ         :: 𝕣1  
    Lx         :: 𝕣1  
    Lr         :: Base.RefValue{𝕣}  
    La         :: 𝕣1  
    Lλx        :: Sparse𝕣2 
    Lλa        :: Sparse𝕣2 
    Lxx        :: Sparse𝕣2 
    Lxr        :: 𝕣1 
    Lrr        :: Base.RefValue{𝕣} 
    Lax        :: Sparse𝕣2 
    Lar        :: 𝕣1 
    Laa        :: Sparse𝕣2 

    ming      :: 𝕣
    minλ      :: 𝕣
    Σλg       :: 𝕣
    npos      :: 𝕫
    # down
    c         :: @NamedTuple{a₁::𝕣, a₂::𝕣, a₃::𝕣, b₁::𝕣, b₂::𝕣, b₃::𝕣}
    firstiter :: 𝕓   
    line      :: 𝕓
end   

function prepare(::Type{AssemblySweepXA{ORDER}},model,dis) where{ORDER}
    Λdofgr             = allΛdofs(model,dis)
    Xdofgr             = allXdofs(model,dis) 
    Adofgr             = allAdofs(model,dis)
    nΛdof              = getndof(Λdofgr)
    nXdof              = getndof(Xdofgr)
    nAdof              = getndof(Adofgr)
    narray,neletyp     = 12,getneletyp(model)
    asm                = Matrix{𝕫2}(undef,narray,neletyp)  # asm[iarray,ieletyp][ieledof,iele]
    Lλ                 = asmvec!(view(asm, 1,:),Λdofgr,dis)
    Lx                 = asmvec!(view(asm, 2,:),Xdofgr,dis)
    Lr                 = Ref(0.)
    La                 = asmvec!(view(asm, 3,:),Adofgr,dis)
    Lλx                = asmmat!(view(asm, 4,:),view(asm,1,:),view(asm,2,:),nXdof,nXdof)
    Lλa                = asmmat!(view(asm, 5,:),view(asm,1,:),view(asm,4,:),nXdof,nAdof)
    Lxx                = asmmat!(view(asm, 6,:),view(asm,2,:),view(asm,2,:),nXdof,nXdof)
    Lxr                = asmvec!(view(asm, 7,:),Xdofgr,dis) 
    Lrr                = Ref(0.)
    Lax                = asmmat!(view(asm, 8,:),view(asm,4,:),view(asm,2,:),nAdof,nXdof)
    Lar                = asmvec!(view(asm, 9,:),Adofgr,dis)  
    Laa                = asmmat!(view(asm,10,:),view(asm,4,:),view(asm,4,:),nAdof,nAdof)

    out                = AssemblySweepXA{ORDER,typeof(Lλ),typeof(Lλx)}(Lλ,Lx,Lr,La,Lλx,Lλa,Lxx,Lxr,Lrr,Lax,Lar,Laa,
                                                                       ∞,∞,0.,0,(a₁=0.,a₂=0.,a₃=0.,b₁=0.,b₂=0.,b₃=0.),false,false) 
    return out,asm,Xdofgr
end
function zero!(out::AssemblySweepXA) # TODO
    zero!(out.Lλ )
    zero!(out.Lx )
    zero!(out.Lr )
    zero!(out.La )
    zero!(out.Lλx)
    zero!(out.Lλa)
    zero!(out.Lxx)
    zero!(out.Lxr)
    zero!(out.Lrr)
    zero!(out.Lax)
    zero!(out.Lar)
    zero!(out.Laa)
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

#=
REPRISE
1) test the whole assembly process
2) solver
3) this way of adiffing is more readable than DirectXUA/addin!.  
   Is there a performance penalty to SweepXA's style? Make DirectXUA (and other solvers' addin!) more readable?
4) It seems that DirectXUA/addin! uses adiff to the second order also when only vectors are required.  This would be very significant for FreqXU.
=#

function addin!{:newmark}(out::AssemblySweepXA,asm,iele,scale,eleobj,Λ,X::NTuple{Nxder,<:SVector{Nx}},U,A::SVector{Na},t,SP,dbg) where{Nxder,Nx,Na}
    a₁,a₂,a₃,b₁,b₂,b₃ = out.c.a₁,out.c.a₂,out.c.a₃,out.c.b₁,out.c.b₂,out.c.b₃
    Nz                = 2Nx+Na
    iΛ                = SVector{Nx  ,𝕫}(    1: Nx  )
    iX                = SVector{Nx  ,𝕫}( Nx+1:2Nx  )
    iA                = SVector{Na  ,𝕫}(2Nx+1: Nz  )
    ir                = SVector{1   ,𝕫}(       Nz+1)
    s                 = SVector{Nz+1,𝕣}(scale.Λ...,scale.X...,1.,scale.A...)
    δZ                = δ{1,Nz+1,𝕣}(s) + δ{2,Nz+1,𝕣}(s)      
    δΛ                = δZ[iΛ]        
    δX                = δZ[iX]        
    δA                = δZ[iA]        
    δr                = δZ[ir]     # Newmark-β special: we need C⋅a and M⋅b
    x,x′,x″           = ∂0(X),∂1(X),∂2(X)
    a                 = a₂*x′ + a₃*x″
    b                 = b₂*x′ + b₃*x″
    vx                = x     +    δX
    vx′               = x′    + a₁*δX + a*δr 
    vx″               = x″    + b₁*δX + b*δr
    L,FB              = getlagrangian(eleobj,Λ+δΛ,(vx,vx′,vx″),U,A+δA,t,SP,dbg)
    ∇L                = ∂{2,Nz+1}(L)
    add_value!(      out.Lλ , asm[ 1], iele, ∇L, ia=iΛ        )  # Lλ  = R    
    add_∂!{1,:minus}(out.Lλ , asm[ 1], iele, ∇L, ia=iΛ, ida=ir)  # Lλ -=   C⋅a + M⋅b   
    add_value!(      out.Lx , asm[ 2], iele, ∇L, ia=iX        )  # Lx    
    add_value!(      out.Lr ,                ∇L, ia=ir        )     
    add_value!(      out.La , asm[ 3], iele, ∇L, ia=iA        )             
    add_∂!{1       }(out.Lλx, asm[ 4], iele, ∇L, ia=iΛ, ida=iX)  # Lλx = K + a₁C + b₁M - there is no Lλr
    add_∂!{1       }(out.Lλa, asm[ 5], iele, ∇L, ia=iΛ, ida=iA)    
    add_∂!{1       }(out.Lxx, asm[ 6], iele, ∇L, ia=iX, ida=iX)  
    add_∂!{1       }(out.Lxr, asm[ 7], iele, ∇L, ia=iX, ida=ir)  
    add_∂!{1       }(out.Lrr,                ∇L, ia=ir, ida=ir)   
    add_∂!{1       }(out.Lax, asm[ 8], iele, ∇L, ia=iA, ida=iX)  
    add_∂!{1       }(out.Lar, asm[ 9], iele, ∇L, ia=iA, ida=ir)  
    add_∂!{1       }(out.Laa, asm[10], iele, ∇L, ia=iA, ida=iA)  
end
function addin!{:iter}(out::AssemblySweepXA{ORDER},asm,iele,scale,eleobj,Λ,X::NTuple{Nxder,<:SVector{Nx}},U,A::SVector{Na},t,SP,dbg) where{ORDER,Nxder,Nx,Na}
    a₁,b₁             = out.c.a₁,out.c.b₁₃
    Nz                = 2Nx+Na
    iΛ                = SVector{Nx ,𝕫}(    1: Nx  )
    iX                = SVector{Nx ,𝕫}( Nx+1:2Nx  )
    iA                = SVector{Na ,𝕫}(2Nx+1: Nz  )
    s                 = SVector{Nzr,𝕣}(scale.Λ...,scale.X...,scale.A...)
    δZ                = δ{1,Nz,𝕣}(s) + δ{2,Nz,𝕣}(s)      
    δΛ                = δZ[iΛ]        
    δX                = δZ[iX]        
    δA                = δZ[iA]        
    if     ORDER==0  L,FB = getlagrangian(eleobj,Λ+δΛ²,(∂0(X)+δX,                         ),U,A+δA,t,SP,dbg)
    elseif ORDER==1  L,FB = getlagrangian(eleobj,Λ+δΛ²,(∂0(X)+δX, ∂1(X)+a₁*δX             ),U,A+δA,t,SP,dbg)
    elseif ORDER==2  L,FB = getlagrangian(eleobj,Λ+δΛ²,(∂0(X)+δX, ∂1(X)+a₁*δX, ∂2(X)+b₁*δX),U,A+δA,t,SP,dbg)
    end
    ∇L²              = ∂{2,Nz}(L)
    add_value!(out.Lλ , asm[ 1], iele, ∇L², ia=iΛ        )  # Lλ  = R    
    add_value!(out.Lx , asm[ 2], iele, ∇L², ia=iX        )  # Lx         
    add_value!(out.La , asm[ 3], iele, ∇L², ia=iA        )             
    add_∂!{1 }(out.Lλx, asm[ 4], iele, ∇L², ia=iΛ ,ida=iX)  # Lλx = K + a₁C + b₁M - there is no Lλr
    add_∂!{1 }(out.Lλa, asm[ 5], iele, ∇L², ia=iΛ ,ida=iA)    
    add_∂!{1 }(out.Lxx, asm[ 6], iele, ∇L², ia=iX ,ida=iX)  
    add_∂!{1 }(out.Lax, asm[ 8], iele, ∇L², ia=iA ,ida=iX)  
    add_∂!{1 }(out.Laa, asm[10], iele, ∇L², ia=iA ,ida=iA)  
end
function addin!{:linesearch}(out::AssemblySweepXA,asm,iele,scale,eleobj,Λ,X,U,A,t,SP,dbg) 
    _,FB             = getlagrangian(eleobj,Λ,X,U,A,t,SP,dbg)
    lineFB!(out,FB)
end
function addin!{mission}(out::AssemblySweepXA,asm,iele,scale,eleobj::Acost,A::SVector{Na},dbg) where{Na,mission} # addin Atarget element
    A∂  = SVector{Na,∂ℝ{2,Na,∂ℝ{1,Na,𝕣}}}(∂²ℝ{1,Na}(A[idof],idof, scale.A[idof])   for idof=1:Na)
    ø   = nothing
    C,_ = lagrangian(eleobj,ø,ø,ø,A∂,ø,ø ,dbg)
    ∇ₐC = ∂{2,Na}(C)
    add_value!(out.La,asm[arrnum(ind.A)],iele,∇ₐC)
    if mission==:matrices
        add_∂!{1}(out.Laa,asm[arrnum(ind.A,ind.A)],iele,∇ₐC)
    end
end
addin!{:linesearch}(args...) = nothing


"""
	SweepXA{ORDER}

A non-linear, time domain solver, that solves the problem time-step by time-step.
Only the `X`-dofs of the model are solved for, while `U`-dofs and `A`-dofs are unchanged.

- `SweepXA{0}` is Newton-Raphson, with feasibility line-search, to handle inequality constraints. 
- `SweepXA{1}` is implicit Euler, with feasibility line-search. 
- `SweepXA{2}` is Newmark-β, with Newton-Raphson iterations and feasibility line search

IMPORTANT NOTE: Muscade does not allow elements to have state variables, for example, plastic strain,
or shear-free position for dry friction.  Where the element implements such physics, this 
is implemented by introducing the state as a degree of freedom of the element, and solving
for its evolution, *even in a static problem*, requires the use of `ORDER≥1`

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
- `β=1/4`,`γ=1/2`     parameters to the Newmark-β algorithm. Dummy if `ORDER<2`
- `maxiter=50`        maximum number of equilibrium iterations at each step.
- `maxΔx=1e-5`        convergence criteria: norm of `X`. 
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

See also: [`solve`](@ref), [`initialize!`](@ref), [`findlastassigned`](@ref), [`study_singular`](@ref), [`DirectXUA`](@ref), [`FreqXU`](@ref)  
"""
struct        SweepXA{ORDER} <: AbstractSolver end
function solve(SX::Type{SweepXA{ORDER}},pstate,verbose,dbg;
                    time::AbstractVector{𝕣},
                    initialstate::State,
                    β::𝕣=1/4,γ::𝕣=1/2,
                    maxiter::ℤ=50,maxΔx::ℝ=1e-5,maxLλ::ℝ=∞,
                    saveiter::𝔹=false,
                    maxLineIter::ℤ=50,sfac::𝕣=.5,γfac::𝕣=.5) where{ORDER}
    model,dis        = initialstate.model,initialstate.dis
    out,asm,Xdofgr   = prepare(AssemblySweepXA{ORDER},model,dis)  
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
        assemble!(out,asm,dis,model,state,(dbg...,solver=:SweepXA,phase=:preliminary,step=step))
        out.ming ≤ 0 && muscadeerror(@sprintf("Initial point is not strictly primal-feasible at step=%3d",step)) # This is going to suck
        out.minλ ≤ 0 && muscadeerror(@sprintf("Initial point is not strictly dual-feasible at step=%3d"  ,step)) # This is going to suck
        state.SP     = (γ=out.Σλg/out.npos * γfac,)   # γ, is in interior point, g(X)*λ=γ
        for iiter    = 1:maxiter
            citer   += 1
            out.firstiter = firstiter = iiter==1
            out.line = false
            assemble!(out,asm,dis,model,state,(dbg...,solver=:SweepXA,step=step,iiter=iiter))
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
                assemble!(out,asm,dis,model,state,(dbg...,solver=:SweepXA,phase=:linesearch,step=step,iiter=iiter,iline=iline))
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
