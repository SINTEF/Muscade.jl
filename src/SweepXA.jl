### Assembler

mutable struct AssemblySweepXA{ORDER} <: Assembly
    # up
    Lλ         :: 𝕣1  
    Lx         :: 𝕣1  
    Lr         :: 𝕣0   
    La         :: 𝕣1  
    Lλx        :: Sparse𝕣2 
    Lλa        :: 𝕣2 
    Lxx        :: Sparse𝕣2 
    Lxr        :: 𝕣1 
    Lrr        :: 𝕣0 
    Lax        :: 𝕣2 
    Lar        :: 𝕣1 
    Laa        :: 𝕣2 

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
    Lλa                = asmfullmat!(view(asm, 5,:),view(asm,1,:),view(asm,3,:),nXdof,nAdof)  
    Lxx                = asmmat!(view(asm, 6,:),view(asm,2,:),view(asm,2,:),nXdof,nXdof)
    Lxr                = asmvec!(view(asm, 7,:),Xdofgr,dis) 
    Lrr                = 𝕣0()
    Lax                = asmfullmat!(view(asm, 8,:),view(asm,3,:),view(asm,2,:),nAdof,nXdof)  
    Lar                = asmvec!(view(asm, 9,:),Adofgr,dis)  
    Laa                = asmfullmat!(view(asm,10,:),view(asm,3,:),view(asm,3,:),nAdof,nAdof)

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
3)  write specific addiff for ElementCost++
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
    a₁,b₁             = out.c.a₁,out.c.b₁
    x,x′,x″,λ         = ∂0(X),∂1(X),∂2(X),∂0(Λ)
    δΛ,δX,δA          = reδ{2}((;Λ=λ,X=x,A),(;Λ=scale.Λ,X=scale.X,A=scale.A)) 
    iΛ,iX,iA,Nz       = revariate_indices(λ,x,A) 
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
function addin!{Amission}(out::AssemblySweepXA,asm,iele,scale,eleobj::Acost,A::SVector{Na},dbg) where{Na,Amission} # addin Atarget element
    d      = revariate{2}((;A),(;A=scale.A)) # careful: revariate returns a NamedTuple
    ø      = nothing
    C,_    = lagrangian(eleobj,ø,ø,ø,d.A,ø,ø ,dbg)
    ∇ₐC    = ∂{2,Na}(C)
    add_value!(out.La ,asm[ 3],iele,∇ₐC)
    add_∂!{1 }(out.Laa,asm[10],iele,∇ₐC)
end

function showstates(state)
        a = state[1].A
        x = [state[istep].X[1][1] for istep=1:length(state)-2]
        λ = [state[istep].Λ[1][1] for istep=1:length(state)-2]
        @show λ,x,a 
end


struct   Newmarkβincrement!{ORDER} end
function Newmarkβincrement!{2}(state,Δx ,Xdofgr,c,firstiter, a,b,x′,x″,Δx′,Δx″) # x′, x″ are just mutable memory, neither input nor output.
    a₁,a₂,a₃,b₁,b₂,b₃ = c.a₁,c.a₂,c.a₃,c.b₁,c.b₂,c.b₃

    if firstiter
        getdof!(state,1,x′,Xdofgr) 
        getdof!(state,2,x″,Xdofgr) 
        a       .= a₂*x′.+ a₃*x″ 
        b       .= b₂*x′.+ b₃*x″
        Δx′     .= a₁*Δx .- a
        Δx″     .= b₁*Δx .- b
    else
        Δx′     .= a₁*Δx 
        Δx″     .= b₁*Δx 
    end
    increment!(state,1,Δx ,Xdofgr)
    increment!(state,2,Δx′,Xdofgr)
    increment!(state,3,Δx″,Xdofgr)
end
function Newmarkβincrement!{1}(state,Δx ,Xdofgr,c,_,Δx′,args...)
    Δx′      .= c.a₁*Δx            
    increment!(state,1,Δx ,Xdofgr)
    increment!(state,2,Δx′,Xdofgr)
end
function Newmarkβincrement!{0}(state,Δx ,Xdofgr,args...)
    increment!(state,1,Δx ,Xdofgr)
end
function Newmarkβcoefficients(order,Δt,β,γ)
    if     order==0 (a₁=0.      , a₂=0. , a₃=0.         , b₁=0.        , b₂=0.      , b₃=0.  )
    elseif order==1 (a₁=1/Δt    , a₂=0  , a₃=0.         , b₁=0.        , b₂=0.      , b₃=0.  )
    elseif order==2 (a₁=γ/(β*Δt), a₂=γ/β, a₃=(γ/2β-1)*Δt, b₁=1/(β*Δt^2), b₂=1/(β*Δt), b₃=1/2β) # γ, as in Newmark's β and γ
    end
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

    state            = OffsetVector{State{ORDER+1,ORDER+1,1}}(0,nstep+1)
    pstate[]         = state 

    state[0] = s     = State{ORDER+1,ORDER+1,1}(copy(initialstate)) 
    tfinal           = time[nstep] + time[1] - state[0].time
    for istep        = 1:nstep # share A (and U, which won't be touched), and set any time derivatives to zero
        state[istep] = State{ORDER+1,ORDER+1,1}(time[istep],(deepcopy(s.Λ[1]),),(deepcopy(s.X[1]),),s.U,s.A,s.SP,s.model,s.dis)
    end 
    state[nstep+1]   = State{ORDER+1,ORDER+1,1}(tfinal     ,(deepcopy(s.Λ[1]),),(deepcopy(s.X[1]),),s.U,s.A,s.SP,s.model,s.dis)
    buffer           = ntuple(i->𝕣1(undef,nXdof), 6)  
    ΔXₐ              = [𝕣2(undef,nXdof,nAdof) for istep=1:nstep]     
    ΔX               = [𝕣1(undef,nXdof      ) for istep=1:nstep]    
    La♯              =  𝕣1(undef,nAdof      )
    Laa♯             =  𝕣2(undef,nAdof,nAdof)
    Lx♯              =  𝕣1(undef,nXdof      ) 
    ΔΛ               =  𝕣1(undef,nXdof      )
    ΔA               =  𝕣1(undef,nAdof      )
    δXᵃ              =  𝕣1(undef,nXdof      )
    LxxΔx            =  𝕣1(undef,nXdof      )
    LxxΔxₐ           =  𝕣2(undef,nXdof,nAdof)
    LxΔxₐ            =  𝕣1(undef,nAdof      )
    LaxΔxₐ           =  𝕣2(undef,nAdof,nAdof)
    LaxΔx            =  𝕣1(undef,nAdof      )
    ΔxₐLxxΔx         =  𝕣1(undef,nAdof      )
    ΔxₐLxxΔxₐ        =  𝕣2(undef,nAdof,nAdof)

    local Lλx # Lλx scopes the function, although it's going to be initialised in a nested scope
    @printf "As received"
    showstates(state)
    # warming up
    for istep        = 1:nstep
        Δt⁻          = state[istep  ].time-state[istep-1].time
        Δt⁻≤0 && ORDER>0 && muscadeerror(@sprintf("Time step length not strictly positive at istep=%3d",istep))
        c⁻           = Newmarkβcoefficients(ORDER,Δt⁻,β,γ)
        outX.c       = c⁻  

        for ider ∈ 1:ORDER+1
            state[istep].X[ider] .= state[istep-1].X[ider]   
        end

        # std Newmark-β
        for iXiter   = 1:maxXiter
            firstXiter = iXiter==1 
            if ORDER==2 && firstXiter assemble!{:newmark}(outX,asmX,dis,model,state[istep],(dbg...,solver=:SweepXA,phase=:warmup,step=step,iXiter=iXiter))
            else                      assemble!{:iter   }(outX,asmX,dis,model,state[istep],(dbg...,solver=:SweepXA,phase=:warmup,step=step,iXiter=iXiter))
            end
            try if  firstXiter Lλx = lu(outX.Lλx) 
            else               lu!(Lλx, outX.Lλx) 
            end catch;         muscadeerror(@sprintf("Lλx matrix factorization failed at warm-up istep=%i, iXiter=%i",istep,iXiter)) end
            ΔX[1]      .= Lλx\-outX.Lλ
            Newmarkβincrement!{ORDER}(state[istep],ΔX[1] ,Xdofgr,outX.c,firstXiter,buffer...) 
            δx²,Lλ²     = sum(ΔX[1].^2),sum(outX.Lλ.^2)
            cXiter     += 1
            if δx²≤cΔx² && Lλ²≤cLλ² 
                verbose && @printf "    At warm-up, step %3d converged in %3d X-iterations. |ΔX|=%7.1e |Lλ|=%7.1e\n" istep iXiter √(δx²) √(Lλ²)
                break#out of the iXiter loop
            end
            iXiter==maxXiter && muscadeerror(@sprintf("no X-convergence at warm-up, istep=%3d after %3d X-iterations |ΔX|=%g / %g, |Lλ|=%g / %g",istep,iXiter,√(δx²),maxΔx,√(Lλ²)^2,maxLλ))
        end
    end

    @printf "After warmup"
    showstates(state)

    for iAiter = 1:maxAiter
        nz = 2*nstep+1
        M = zeros(nz,nz)
        V = zeros(nz)




        @printf("--- iAiter = %i ---\n",iAiter)
        assembleA!{:ok}(outXA,asmXA,dis,model,state[0],(dbg...,solver=:SweepXA,phase=:Acost,iAiter=iAiter))
        La♯             .= outXA.La   
        Laa♯            .= outXA.Laa  

        @show :assembleA,La♯

        M[nz,nz] = Laa♯[1,1]
        V[nz]    = La♯[1]

        # forward sweep
        for istep        = 1:nstep
        @printf("-- istep = %i --\n",istep)
            Δt⁻          = state[istep  ].time-state[istep-1].time
            Δt⁻≤0 && ORDER>0 && muscadeerror(@sprintf("Time step length not strictly positive at istep=%3d",istep))
            c⁻           = Newmarkβcoefficients(ORDER,Δt⁻,β,γ)
            outXA.c      = c⁻  
            # for ider ∈ 1:ORDER+1
            #     state[istep].X[ider] .= state[istep-1].X[ider]   
            # end

            # sensitivity
         #   @show state[istep].X[1]
            # if ORDER==2 assemble!{:newmark}(outXA,asmXA,dis,model,state[istep],(dbg...,solver=:SweepXA,phase=:sensitivity,iAiter=iAiter,step=step))
            # else        assemble!{:iter   }(outXA,asmXA,dis,model,state[istep],(dbg...,solver=:SweepXA,phase=:sensitivity,iAiter=iAiter,step=step))
            # end
            assemble!{:iter   }(outXA,asmXA,dis,model,state[istep],(dbg...,solver=:SweepXA,phase=:sensitivity,iAiter=iAiter,step=step))
            try if iAiter==1  Lλx = lu(outXA.Lλx) 
            else              lu!(Lλx, outXA.Lλx) 
            end catch;        muscadeerror(@sprintf("Lλx matrix factorization failed at iAiter=%3d, istep=%i, iXiter=%i",iAiter,istep,iXiter)) end

            ΔX[ istep] .= Lλx\-outXA.Lλ  # increment since after the X-iterations
            ΔXₐ[istep] .= Lλx\-outXA.Lλa 
 
            # TODO causing allocations here?
            LxxΔx        .=                          outXA.Lxx  ∘₁ ΔX[ istep] .+ outXA.Lr   # x
            LxxΔxₐ       .=                          outXA.Lxx  ∘₁ ΔXₐ[istep] .+ outXA.Lxr  # xa
            LxΔxₐ        .=                          outXA.Lx   ∘₁ ΔXₐ[istep] .+ outXA.Lr   # a
            LaxΔxₐ       .=                          outXA.Lax  ∘₁ ΔXₐ[istep] .+ outXA.Lar  # aa 
            LaxΔx        .=                          outXA.Lax  ∘₁ ΔX[ istep] .+ outXA.Lar  # a
            ΔxₐLxxΔx     .= ΔXₐ[istep]' ∘₁ LxxΔx  .+ outXA.Lxr' ∘₁ ΔX[ istep] .+ outXA.Lrr  # a
            ΔxₐLxxΔxₐ    .= ΔXₐ[istep]' ∘₁ LxxΔxₐ .+ outXA.Lxr' ∘₁ ΔXₐ[istep] .+ outXA.Lrr  # aa  # is symmetric

            # @show outXA.Lxx
            # @show outXA.Lax
            # @show outXA.Lx
            @show :before,La♯
            # @show ΔxₐLxxΔx
            # @show LaxΔx
            # @show LxΔxₐ

            M[nstep+istep,nstep+istep] = outXA.Lxx[ 1,1]
            M[nstep+istep,nz         ] = outXA.Lax'[1,1]
            M[nz         ,nstep+istep] = outXA.Lax[ 1,1]
            V[nstep+istep            ] = outXA.Lx[  1  ]
            @show V[nz                     ] = La♯[  1  ]

            La♯         .+= ΔxₐLxxΔx  .+ LaxΔx  .+ LxΔxₐ                                    # a
            Laa♯        .+= ΔxₐLxxΔxₐ .+ LaxΔxₐ .+ LaxΔxₐ'                                  # aa   
#            @show ΔX[ istep],ΔXₐ[istep]
#            @show Matrix(outXA.Lxx),Matrix(outXA.Lax),outXA.Lx,outXA.Lλa,outXA.Lλ
#            @show ΔxₐLxxΔx,LaxΔx,LxΔxₐ



        end # istep

        # update A
        ΔA      .= Laa♯\-La♯
        ΔA²,La²  = sum(ΔA.^2),sum(La♯.^2)
        verbose && @printf "    In A-iteration %3d, |ΔA|=%7.1e |La♯|=%7.1e\n" iAiter √(ΔA²) √(La²)

        # backward sweep 
        for istep = nstep:-1:1
            Δt⁻          = state[istep  ].time-state[istep-1].time
            Δt⁺          = state[istep+1].time-state[istep  ].time
            c⁻           = Newmarkβcoefficients(ORDER,Δt⁻,β,γ)    
            c⁺           = Newmarkβcoefficients(ORDER,Δt⁺,β,γ)
            outXA.c      = c⁺  # TODO optimize: outXA is overkill, but I need Lax
            for ider ∈ 1:ORDER+1
                state[istep].Λ[ider] .= state[istep+1].Λ[ider]   
            end

            if ORDER==2  assemble!{:newmark}(outXA,asmXA,dis,model,state[istep],(dbg...,solver=:SweepXA,phase=:backward,iAiter=iAiter,step=step))
            else         assemble!{:iter   }(outXA,asmXA,dis,model,state[istep],(dbg...,solver=:SweepXA,phase=:backward,iAiter=iAiter,step=step))
            end

#            @show ΔXₐ[istep] ∘₁ ΔA    
            ΔX[istep]  .+= ΔXₐ[istep] ∘₁ ΔA # double sign swap here!!!        

            LxxΔx       .=                    outXA.Lxx  ∘₁ ΔX[istep] .+ outXA.Lr   
            Lx♯         .= outXA.Lx + LxxΔx + outXA.Lax' ∘₁ ΔA 
            ΔΛ          .= outXA.Lλx'\-Lx♯  

            #@show ΔΛ,ΔX[ istep],ΔA

            M[nz         ,      istep] = outXA.Lλa'[1,1]
            M[      istep,nstep+istep] = outXA.Lλx[ 1,1]
            M[      istep,nz         ] = outXA.Lλa[ 1,1]
            M[nstep+istep,      istep] = outXA.Lλx'[1,1]
            V[      istep            ] = outXA.Lλ[  1  ]

            @show ΔΛ
            @show ΔX[istep]
            @show ΔA


            Newmarkβincrement!{ORDER}(state[istep],ΔX[istep],Xdofgr,c⁻,false,buffer...) 
            Newmarkβincrement!{ORDER}(state[istep],ΔΛ       ,Λdofgr,c⁺,false,buffer...) 
        end
        increment!(state[1],1,ΔA,Adofgr) # state[i].A === state[j].A

     #   @printf "After Aiter\n"



        #@show M
        #@show V
        @show M\-V
        @show M
        @show V

        showstates(state)

        # Aiter convergence
        if ΔA²≤cΔa² && La²≤cLa² 
            verbose && @printf "    SweepXA converged in %3d A-iterations. |ΔA|=%7.1e |La|=%7.1e\n" iAiter √(ΔA²) √(La²)
            cAiter = iAiter
            break#out of the iAiter loop
        end
        iAiter==maxAiter && muscadeerror(@sprintf("no convergence of SweepXA after %3d A-iterations |ΔA|=%g / %g, |La|=%g / %g",iAiter,√(ΔA²),maxΔa,√(La²)^2,maxLa))

    end 
    verbose && @printf "\n    nel=%d, nXdof=%d, nstep=%d, nAiter=%d, ΣnXiter=%d, mean(nXiter)=%d\n" getnele(model) getndof(Xdofgr) length(time) cAiter cXiter cXiter/length(time)
    return
end
