### Assembler

mutable struct AssemblySweepXA{ORDER} <: Assembly
    # up
    Lλ         :: 𝕣1  
    Lx         :: 𝕣1  
    Lr         :: 𝕣0   
    La         :: 𝕣1  
    Lλx        :: Sparse𝕣2 
    Lλa        :: Sparse𝕣2 
    Lxx        :: Sparse𝕣2 
    Lxr        :: 𝕣1 
    Lrr        :: 𝕣0 
    Lax        :: Sparse𝕣2 
    Lar        :: 𝕣1 
    Laa        :: Sparse𝕣2 

    ming      :: 𝕣
    minλ      :: 𝕣
    Σλg       :: 𝕣
    npos      :: 𝕫
    # down
    c         :: @NamedTuple{a₁::𝕣, a₂::𝕣, a₃::𝕣, b₁::𝕣, b₂::𝕣, b₃::𝕣}
end   

function prepare(::Type{AssemblySweepXA{ORDER}},model,dis) where{ORDER}
    Λdofgr             = allΛdofs(model,dis)
    Xdofgr             = allXdofs(model,dis) 
    Adofgr             = allAdofs(model,dis)
    nΛdof              = getndof(Λdofgr)
    nXdof              = getndof(Xdofgr)
    nAdof              = getndof(Adofgr)
    narray,neletyp     = 10,getneletyp(model)
    asm                = Matrix{𝕫2}(undef,narray,neletyp)  # asm[iarray,ieletyp][ieledof,iele]
    Lλ                 = asmvec!(view(asm, 1,:),Λdofgr,dis)
    Lx                 = asmvec!(view(asm, 2,:),Xdofgr,dis)
    Lr                 = 𝕣0()
    La                 = asmvec!(view(asm, 3,:),Adofgr,dis)
    Lλx                = asmmat!(view(asm, 4,:),view(asm,1,:),view(asm,2,:),nXdof,nXdof)
    Lλa                = asmmat!(view(asm, 5,:),view(asm,1,:),view(asm,3,:),nXdof,nAdof)
    Lxx                = asmmat!(view(asm, 6,:),view(asm,2,:),view(asm,2,:),nXdof,nXdof)
    Lxr                = asmvec!(view(asm, 7,:),Xdofgr,dis) 
    Lrr                = 𝕣0()
    Lax                = asmmat!(view(asm, 8,:),view(asm,3,:),view(asm,2,:),nAdof,nXdof)
    Lar                = asmvec!(view(asm, 9,:),Adofgr,dis)  
    Laa                = asmmat!(view(asm,10,:),view(asm,3,:),view(asm,3,:),nAdof,nAdof)

    out                = AssemblySweepXA{ORDER}(Lλ,Lx,Lr,La,Lλx,Lλa,Lxx,Lxr,Lrr,Lax,Lar,Laa, ∞,∞,0.,0, (a₁=0.,a₂=0.,a₃=0.,b₁=0.,b₂=0.,b₃=0.)) 
    return out,asm,Λdofgr,Xdofgr,Adofgr
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

#=
REPRISE
2) solver
3) use revariate, and write specific addiff for ElementCost++
=#

function addin!{:newmark}(out::AssemblySweepXA,asm,iele,scale,eleobj,Λ,X::NTuple{Nxder,<:SVector{Nx}},U,A::SVector{Na},t,SP,dbg) where{Nxder,Nx,Na}
    a₁,a₂,a₃,b₁,b₂,b₃ = out.c.a₁,out.c.a₂,out.c.a₃,out.c.b₁,out.c.b₂,out.c.b₃
    x,x′,x″,λ         = ∂0(X),∂1(X),∂2(X),∂0(Λ)
    δΛ,δX,δA,δr       = reδ{2}((;Λ=λ,X=x,A,r=0.),(;Λ=scale.Λ,X=scale.X,A=scale.A,r=1.)) 
    iΛ,iX,iA,ir,Nz    = revariate_indices(λ,x,A,0.) 
    a                 = a₂*x′ + a₃*x″
    b                 = b₂*x′ + b₃*x″
    vx                = x     +    δX
    vx′               = x′    + a₁*δX + a*δr  
    vx″               = x″    + b₁*δX + b*δr 
    L,FB              = getlagrangian(eleobj,λ+δΛ,(vx,vx′,vx″),U,A+δA,t,SP,dbg)
    ∇L                = ∂{2,Nz}(L)
    add_value!(      out.Lλ , asm[ 1], iele, ∇L, iΛ    )  # Lλ  = R    
    add_∂!{1,:minus}(out.Lλ , asm[ 1], iele, ∇L, iΛ, ir)  # Lλ -=   C⋅a + M⋅b   
    add_value!(      out.Lx , asm[ 2], iele, ∇L, iX    )  # Lx    
    add_value!(      out.Lr ,                ∇L, ir    )     
    add_value!(      out.La , asm[ 3], iele, ∇L, iA    )             
    add_∂!{1       }(out.Lλx, asm[ 4], iele, ∇L, iΛ, iX)  # Lλx = K + a₁C + b₁M - there is no Lλr
    add_∂!{1       }(out.Lλa, asm[ 5], iele, ∇L, iΛ, iA)    
    add_∂!{1       }(out.Lxx, asm[ 6], iele, ∇L, iX, iX)  
    add_∂!{1       }(out.Lxr, asm[ 7], iele, ∇L, iX, ir) 
    add_∂!{1       }(out.Lrr,                ∇L, ir, ir)   
    add_∂!{1       }(out.Lax, asm[ 8], iele, ∇L, iA, iX)  
    add_∂!{1       }(out.Lar, asm[ 9], iele, ∇L, iA, ir)  
    add_∂!{1       }(out.Laa, asm[10], iele, ∇L, iA, iA)  
end
function addin!{:iter}(out::AssemblySweepXA{ORDER},asm,iele,scale,eleobj,Λ,X::NTuple{Nxder,<:SVector{Nx}},U,A::SVector{Na},t,SP,dbg) where{ORDER,Nxder,Nx,Na}
    a₁,b₁             = out.c.a₁,out.c.b₁₃
    x,x′,x″,λ         = ∂0(X),∂1(X),∂2(X),∂0(Λ)
    δΛ,δX,δA          = reδ{2}((;Λ=λ,X=x,A),(;Λ=scale.Λ,X=scale.X,A=scale.A)) 
    if     ORDER==0  L,FB = getlagrangian(eleobj,λ+δΛ,(x+δX,                   ),U,A+δA,t,SP,dbg)
    elseif ORDER==1  L,FB = getlagrangian(eleobj,λ+δΛ,(x+δX, x′+a₁*δX          ),U,A+δA,t,SP,dbg)
    elseif ORDER==2  L,FB = getlagrangian(eleobj,λ+δΛ,(x+δX, x′+a₁*δX, x″+b₁*δX),U,A+δA,t,SP,dbg)
    end
    ∇L               = ∂{2,Nz}(L)
    add_value!(out.Lλ , asm[ 1], iele, ∇L, iΛ    )  # Lλ  = R    
    add_value!(out.Lx , asm[ 2], iele, ∇L, iX    )  # Lx         
    add_value!(out.La , asm[ 3], iele, ∇L, iA    )             
    add_∂!{1 }(out.Lλx, asm[ 4], iele, ∇L, iΛ ,iX)  # Lλx = K + a₁C + b₁M - there is no Lλr
    add_∂!{1 }(out.Lλa, asm[ 5], iele, ∇L, iΛ ,iA)    
    add_∂!{1 }(out.Lxx, asm[ 6], iele, ∇L, iX ,iX)  
    add_∂!{1 }(out.Lax, asm[ 8], iele, ∇L, iA ,iX)  
    add_∂!{1 }(out.Laa, asm[10], iele, ∇L, iA ,iA)  
end
function addin!{mission}(out::AssemblySweepXA,asm,iele,scale,eleobj::Acost,A::SVector{Na},dbg) where{Na,mission} # addin Atarget element
    ∂A     = revariate{2}((;A),(;A=scale.A)) 
    ø      = nothing
    C,_    = lagrangian(eleobj,ø,ø,ø,∂A,ø,ø ,dbg)
    ∇ₐC    = ∂{2,Na}(C)
    add_value!(out.La ,asm[ 3],iele,∇ₐC)
    add_∂!{1 }(out.Laa,asm[10],iele,∇ₐC)
end


"""
	SweepXA{ORDER}

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
- `β=1/4`,`γ=1/2`     parameters to the Newmark-β algorithm. Dummy if `ORDER<2`
- `maxiter=50`        maximum number of equilibrium iterations at each step.
- `maxΔx=1e-5`        convergence criteria: norm of `X`. 
- `maxLλ=∞`           convergence criteria: norm of the residual. 
- `saveiter=false`    set to true so that output `states` contains the state
                      at the iteration of the last step analyseδ  Useful to study
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
                    maxXiter::ℤ=50,maxΔx::ℝ=1e-5,maxLλ::ℝ=∞,
                    maxAiter::ℤ=50,maxΔa::ℝ=1e-5,maxLa::ℝ=∞) where{ORDER}

    model,dis        = initialstate.model,initialstate.dis
    outX ,asmX ,       Xdofgr          = prepare(AssemblySweepX{ ORDER},model,dis)  
    outXA,asmXA,Λdofgr,Xdofgr,Adofgr   = prepare(AssemblySweepXA{ORDER},model,dis)  
    nXdof            = getndof(Xdofgr)
    nAdof            = getndof(Adofgr)
    nstep            = length(time)
    cΔx²,cLλ²        = maxΔx^2,maxLλ^2
    cΔa²,cLa²        = maxΔa^2,maxLa^2
    cXiter           = 0
    cAiter           = 0

    state            = OffsetVector{State{1,ORDER+1,1}}(0,nstep)
    pstate[]         = state 
    state[0] = s     = State{1,ORDER+1,1}(copy(initialstate)) 
    for istep        = 1:nstep # share A (and U, which won't be touched)
        state[istep] = State{1,ORDER+1,1}(time[istep],copy(s.Λ),copy(s.X),s.U,s.A,s.SP,s.model,s.dis)
    end 
    
    buffer           = ntuple(i->𝕣1(undef,nXdof), 6)  
    if ORDER≥1    x′ = 𝕣1(undef,nXdof       ) end 
    if ORDER≥2    x″ = 𝕣1(undef,nXdof       ) end 
    Lλx              = Vector{LU𝕣     }(undef,nstep)
    Lxx              = Vector{Sparse𝕣2}(undef,nstep)
    Lx               = [𝕣1(undef,nXdof      ) for istep=1:nstep]    
    Lax              = [𝕣2(undef,nAdof,nXdof) for istep=1:nstep]   
    Δxₐ              = [𝕣2(undef,nXdof,nAdof) for istep=1:nstep]     
    Δx               = [𝕣1(undef,nXdof      ) for istep=1:nstep]    
    LxxΔx            = [𝕣1(undef,nXdof      ) for istep=1:nstep] 
    La♯              =  𝕣1(undef,nAdof      )
    Laa♯             =  𝕣2(undef,nAdof,nAdof)
    Lx♯              =  𝕣1(undef,nXdof      ) 
    δx               =  𝕣1(undef,nXdof      )
    ΔΛ               =  𝕣1(undef,nXdof      )
    LxxΔxₐ           =  𝕣2(undef,nXdof,nAdof)
    LxΔxₐ            =  𝕣1(undef,nAdof      )
    LaxΔxₐ           =  𝕣2(undef,nAdof,nAdof)
    LaxΔx            =  𝕣1(undef,nAdof      )
    ΔxₐLxxΔx         =  𝕣1(undef,nAdof      )
    ΔxₐLxxΔxₐ        =  𝕣2(undef,nAdof,nAdof)

    for iAiter = 1:maxAiter
        assembleA!{:ok}(outXA,asmXA,dis,model,state,(dbg...,solver=:SweepXA,phase=:Acost,iAiter=iAiter))
        La♯ .= outXA.La   # Lₐ*  in the theory
        Laa♯.= outXA.Laa  # Lₐₐ* in the theory

        # forward sweep
        for istep        = 1:nstep
            t            = state[istep  ].time
            oldt         = state[istep-1].time
            Δt           = t-oldt
            Δt≤0 && ORDER>0 && muscadeerror(@sprintf("Time step length not strictly positive at istep=%3d",istep))
            outX.c        = Newmarkβcoefficients(ORDER,Δt,β,γ)

            state[istep].X .= state[istep-1].X   
            Δx[istep]   .= 0.

            # std Newmark-β
            for iXiter   = 1:maxXiter
                firstiter = iXiter==1
                if ORDER==2 && firstiter  assemble!{:newmark}(outX,asmX,dis,model,state,(dbg...,solver=:SweepX,step=step,iiter=iiter))
                else                      assemble!{:iter   }(outX,asmX,dis,model,state,(dbg...,solver=:SweepX,step=step,iiter=iiter))
                end
                try if iAiter==1 && firstiter  Lλx[istep] = lu(outX.Lλx) 
                else                           lu!(Lλx[istep], outX.Lλx) 
                end catch;                     muscadeerror(@sprintf("Lλx matrix factorization failed at Aiter=%3d, istep=%i, iiter=%i",iAiter,istep,iiter)) end
                δx         .= Lλx[istep]\outX.Lλ
                Δx[istep] .+= δx
                Newmarkβdecrement!{ORDER}(state[istep],δx ,Xdofgr,outX.c,outX.firstiter,buffer...) 
                δx²,Lλ²     = sum(δx.^2),sum(outX.Lλ.^2)
                cXiter     += 1
                if δx²≤cΔx² && Lλ²≤cLλ² 
                    verbose && @printf "    In Aiter %3d, step %3d converged in %3d X-iterations. |Δx|=%7.1e |Lλ|=%7.1e\n" iAiter istep iXiter √(δx²) √(Lλ²)
                    break#out of the iXiter loop
                end
                iXiter==maxXiter && muscadeerror(@sprintf("no X-convergence at iAiter=%3d, istep=%3d after %3d X-iterations |Δx|=%g / %g, |Lλ|=%g / %g",iAiter,istep,iXiter,√(δx²),maxΔx,√(Lλ²)^2,maxLλ))
            end

            # sensitivity
            if ORDER==2 && firstiter  assemble!{:newmark}(outXA,asmXA,dis,model,state,(dbg...,solver=:SweepXA,iAiter=iAiter,istep=istep))
            else                      assemble!{:iter   }(outXA,asmXA,dis,model,state,(dbg...,solver=:SweepXA,iAiter=iAiter,istep=istep))
            end
            assemble!{TODO}(outXA,asmXA,dis,model,state,(dbg...,solver=:SweepXA,iAiter=iAiter,istep=istep,iXiter=iXiter))
            Lx[ istep] .= outXA.Lx        # TODO instead of copying, why not let 'outXA' point to the correct memory?
            Lxx[istep]  = copy(outXA.Lxx) # because Lxx is sparse
            Lax[istep] .= outXA.Lax

            try lu!(Lλx[istep], outXA.Lλx) catch; muscadeerror(@sprintf("Lλx matrix factorization failed at Aiter=%3d, istep=%i, sensitivity",iAiter,istep)) end
            Δxₐ[istep] .= Lλx[istep]\outXA.Lλa 
            δx         .= Lλx[istep]\outXA.Lλ  
            Δx[ istep].+= δx
            Newmarkβdecrement!{ORDER}(state[istep],δx ,Xdofgr,outXA.c,firstiter,buffer...) 

            # TODO causing allocations here?
            LxxΔx[istep] .= Lxx[istep]  ∘₁ Δx[   istep]                              .+ outXA.Lr   # x
            LxxΔxₐ       .= Lxx[istep]  ∘₁ Δxₐ[  istep]                              .+ outXA.Lxr  # xa
            LxΔxₐ        .= Lx[ istep]  ∘₁ Δxₐ[  istep]                              .+ outXA.Lr   # a
            LaxΔxₐ       .= Lax[istep]  ∘₁ Δxₐ[  istep]                              .+ outXA.Lar  # aa 
            LaxΔx        .= Lax[istep]  ∘₁ Δx[   istep]                               + outXA.Lar  # a
            ΔxₐLxxΔx     .= Δxₐ[istep]' ∘₁ LxxΔx[istep]  .+ outXA.Lxr' ∘₁ Δx[ istep] .+ outXA.Lrr  # a
            ΔxₐLxxΔxₐ    .= Δxₐ[istep]' ∘₁ LxxΔxₐ        .+ outXA.Lxr' ∘₁ Δxₐ[istep] .+ outXA.Lrr  # aa  # TODO test symmetry
            La♯         .+= ΔxₐLxxΔx  .- LaxΔx  .- LxΔxₐ                                          # a
            Laa♯        .+= ΔxₐLxxΔxₐ .- LaxΔxₐ .- LaxΔxₐ'                                        # aa   
        
        end # istep

        # update A
        ΔA .= Laa♯\La♯
        decrement!(state[istep],1,ΔA,Adofgr) 
        ΔA²,La²  = sum(ΔA.^2),sum(La.^2)

        # backward sweep
        for istep = nstep:-1:1
            δx        .= Δxₐ[istep] ∘₁ ΔA
            Newmarkβdecrement!{ORDER}(state[istep],δx ,Xdofgr,outXA.c,false,buffer...)
            Δx[istep].+= δx
            Lx♯       .= Lx[istep] - LxxΔx[istep] - Lax[istep]' ∘₁ ΔA 
            ΔΛ        .= Lλx[istep]'\Lx♯
            decrement!{ORDER}(state[istep],1,ΔΛ ,Λdofgr) 
        end

        # Aiter convergence
        if ΔA²≤cΔA² && La²≤cLa² 
            verbose && @printf "    SweepXA converged in %3d A-iterations. |ΔA|=%7.1e |La|=%7.1e\n" iAiter √(δA²) √(La²)
            cAiter = iAiter
            break#out of the iAiter loop
        end
        iAiter==maxAiter && muscadeerror(@sprintf("no convergence of SweepXA after %3d A-iterations |ΔA|=%g / %g, |La|=%g / %g",istep,iAiter,√(ΔA²),maxΔa,√(La²)^2,maxLa))

    end 
    verbose && @printf "\n    nel=%d, nXdof=%d, nstep=%d, nAiter, ΣnXiter=%d\n" getnele(model) getndof(Xdofgr) length(time) cAiter cXiter cXiter/length(time)
    return
end
