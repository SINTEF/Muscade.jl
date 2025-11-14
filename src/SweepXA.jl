### Assembler

mutable struct AssemblySweepXA{OX,NDX} <: Assembly
    # up
    Lλ        :: 𝕣1  
    Lx        :: 𝕣1  
    Lr        :: 𝕣0   
    La        :: 𝕣1  
    Lλx       :: Sparse𝕣2 
    Lλa       :: 𝕣2 
    Lxx       :: Sparse𝕣2 
    Lxr       :: 𝕣1 
    Lrr       :: 𝕣0 
    Lax       :: 𝕣2 
    Lar       :: 𝕣1 
    Laa       :: 𝕣2 
    # down
    c         :: Newmarkβcoefficients{OX}
    XorΛ      :: Ref{NTuple{NDX,𝕣1}}
end   

function prepare(::Type{AssemblySweepXA{OX}},model,dis) where{OX}
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

    out                = AssemblySweepXA{OX,OX+1}(Lλ,Lx,Lr,La,Lλx,Lλa,Lxx,Lxr,Lrr,Lax,Lar,Laa, Newmarkβcoefficients{OX}(),Ref{NTuple{OX+1,𝕣1}}()) 
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
end

#=        TODO
solver
write specific adiff for ElementCost
SweepXA for order 0 and 1
Multi load cases        
=#


function addin!{:Xsweep}(out::AssemblySweepXA{2},asm,iele,scale,eleobj,Λ,X::NTuple{Nxder,<:SVector{Nx}},U,A::SVector{Na},t,Δt,SP,dbg) where{Nxder,Nx,Na}
    a₁,a₂,a₃,b₁,b₂,b₃ = out.c.a₁,out.c.a₂,out.c.a₃,out.c.b₁,out.c.b₂,out.c.b₃
    x,x′,x″,λ         = ∂0(X),∂1(X),∂2(X),∂0(Λ)
    x⁻,x′⁻,x″⁻        = X[4],X[5],X[6]
    δΛ,δX,δA,δr       = reδ{2}((;Λ=λ,X=x,A,r=0.),(;Λ=scale.Λ,X=scale.X,A=scale.A,r=1.)) 
    iΛ,iX,iA,ir,Nz    = revariate_indices(λ,x,A,0.) 
    a                 = a₁*(x⁻.-x) + (a₂-1)*x′⁻ +     a₃*x″⁻ + x′      
    b                 = b₁*(x⁻.-x) +     b₂*x′⁻ + (b₃-1)*x″⁻ + x″       
    vx                = x     +    δX
    vx′               = x′    + a₁*δX + a*δr  
    vx″               = x″    + b₁*δX + b*δr 
    L,FB              = getlagrangian(eleobj,λ+δΛ,(vx,vx′,vx″),U,A+δA,t,SP,dbg)
    ∇L                = ∂{2,Nz}(L)
    add_value!(      out.Lλ , asm[ 1], iele, ∇L, iΛ    ;Δt)  # Lλ  = R    
    add_∂!{1,:minus}(out.Lλ , asm[ 1], iele, ∇L, iΛ, ir;Δt)  # Lλ -=   C⋅a + M⋅b   
    add_value!(      out.Lx , asm[ 2], iele, ∇L, iX    ;Δt)     
    add_value!(      out.Lr ,                ∇L, ir    ;Δt)     
    add_value!(      out.La , asm[ 3], iele, ∇L, iA    ;Δt)             
    add_∂!{1       }(out.Lλx, asm[ 4], iele, ∇L, iΛ, iX;Δt)  # Lλx = K + a₁C + b₁M - there is no Lλr
    add_∂!{1       }(out.Lλa, asm[ 5], iele, ∇L, iΛ, iA;Δt)    
    add_∂!{1       }(out.Lxx, asm[ 6], iele, ∇L, iX, iX;Δt)  
    add_∂!{1       }(out.Lxr, asm[ 7], iele, ∇L, iX, ir;Δt) 
    add_∂!{1       }(out.Lrr,                ∇L, ir, ir;Δt)   
    add_∂!{1       }(out.Lax, asm[ 8], iele, ∇L, iA, iX;Δt)  
    add_∂!{1       }(out.Lar, asm[ 9], iele, ∇L, iA, ir;Δt)  
    add_∂!{1       }(out.Laa, asm[10], iele, ∇L, iA, iA;Δt)  
end
function addin!{:Λsweep}(out::AssemblySweepXA{2},asm,iele,scale,eleobj,Λ,X::NTuple{Nxder,<:SVector{Nx}},U,A::SVector{Na},t,Δt,SP,dbg) where{Nxder,Nx,Na}
    a₁,a₂,a₃,b₁,b₂,b₃ = out.c.a₁,out.c.a₂,out.c.a₃,out.c.b₁,out.c.b₂,out.c.b₃
    x,x′,x″,λ,λ′,λ″   = ∂0(X),∂1(X),∂2(X),∂0(Λ),∂1(Λ),∂2(Λ)
    λ⁺,λ′⁺,λ″⁺        = X[4],X[5],X[6]
    δX,δr             = reδ{1}((;X=x,r=0.),(;X=scale.X,r=1.)) 
    iX,ir,Nz          = revariate_indices(x,0.) 
    a                 = a₁*(λ⁺.-λ) + (a₂-1)*λ′⁺ +     a₃*λ″⁺ + λ′      
    b                 = b₁*(λ⁺.-λ) +     b₂*λ′⁺ + (b₃-1)*λ″⁺ + λ″      
    vx                = x  
    vx′               = x′  + a*δr  
    vx″               = x″  + b*δr 
    L,FB              = getlagrangian(eleobj,λ+δΛ,(vx,vx′,vx″),U,A+δA,t,SP,dbg)
    # add_∂!{1}(out.Lx , asm[ 2], iele, L, idvec,iX    ;Δt)    
    # add_∂!{1}(out.Lr ,                L, idvec,ir    ;Δt)     
    ∇L                = ∂{1,Nz}(L)
    add_value!(      out.Lx , asm[ 2], iele, ∇L, iX    ;Δt)    
    add_value!(      out.Lr ,                ∇L, ir    ;Δt)     
end
function addin!{:Acost}(out::AssemblySweepXA,asm,iele,scale,eleobj::Acost,A::SVector{Na},dbg) where{Na} 
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


struct   Newmarkβsweepincrement!{OX} end

function Newmarkβsweepincrement!{2}(state,state⁻,dx ,Xdofgr,c, a,b,x,x′,x″,x⁻,x′⁻,x″⁻,dx′,dx″) # a,b,... are just mutable memory, neither input nor output.
    a₁,a₂,a₃,b₁,b₂,b₃ = c.a₁,c.a₂,c.a₃,c.b₁,c.b₂,c.b₃
    getdof!(state ,0,x  ,Xdofgr) 
    getdof!(state ,1,x′ ,Xdofgr) 
    getdof!(state ,2,x″ ,Xdofgr) 
    getdof!(state⁻,0,x⁻ ,Xdofgr) 
    getdof!(state⁻,1,x′⁻,Xdofgr) 
    getdof!(state⁻,2,x″⁻,Xdofgr) 
    a   .= a₁*(x⁻.-x) + (a₂-1)*x′⁻ +     a₃*x″⁻ + x′      
    b   .= b₁*(x⁻.-x) +     b₂*x′⁻ + (b₃-1)*x″⁻ + x″       
    dx′ .= a₁*dx .- a
    dx″ .= b₁*dx .- b
    increment!(state,1,dx ,Xdofgr)
    increment!(state,2,dx′,Xdofgr)
    increment!(state,3,dx″,Xdofgr)
end
function Newmarkβsweepincrement!{1}(state,dx ,Xdofgr,c,dx′,args...)
    a₁,a₂ = c.a₁,c.a₂
    getdof!(state ,0,x  ,Xdofgr) 
    getdof!(state ,1,x′ ,Xdofgr) 
    getdof!(state⁻,0,x⁻ ,Xdofgr) 
    getdof!(state⁻,1,x′⁻,Xdofgr) 
    a   .= a₁*(x⁻.-x) + (a₂-1)*x′⁻ + x′      
    dx′ .= a₁*dx .- a
    increment!(state,1,dx ,Xdofgr)
    increment!(state,2,dx′,Xdofgr)
end
function Newmarkβsweepincrement!{0}(state,dx ,Xdofgr,args...)
    increment!(state,1,dx ,Xdofgr)
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
struct        SweepXA{OX} <: AbstractSolver end
function solve(SX::Type{SweepXA{OX}},pstate,verbose,dbg;
                    time::AbstractRange{𝕣},
                    initialstate::State,
                    β::𝕣=1/4,γ::𝕣=1/2,
                    maxXiter::ℤ=50,maxΔx::ℝ=1e-5,maxLλ::ℝ=∞,
                    maxAiter::ℤ=50,maxΔa::ℝ=1e-5,maxLa::ℝ=∞) where{OX}

    model,dis        = initialstate.model,initialstate.dis
    outX ,asmX ,       Xdofgr          = prepare(AssemblySweepX{ OX},model,dis)  
    outXA,asmXA,Λdofgr,Xdofgr,Adofgr   = prepare(AssemblySweepXA{OX},model,dis)  
    nXdof            = getndof(Xdofgr)
    nAdof            = getndof(Adofgr)
    nstep            = length(time)
    Δt               = step(time)
    outX.c=outXA.c   = Newmarkβcoefficients{OX}(Δt,β,γ)  
    cΔx²,cLλ²        = maxΔx^2,maxLλ^2
    cΔa²,cLa²        = maxΔa^2,maxLa^2
    cXiter           = 0
    cAiter           = 0

    state            = OffsetVector{State{OX+1,OX+1,1}}(0,nstep+1)
    pstate[]         = view(state.a,2:nstep+1) 

    state[0] = s     = State{OX+1,OX+1,1}(copy(initialstate)) 
    state[0].time    = time[1]-Δt

    for istep        = 1:nstep # share A (and U, which won't be touched), and set any time derivatives to zero
        state[istep] = State{OX+1,OX+1,1}(time[istep]   ,(deepcopy(s.Λ[1]),),(deepcopy(s.X[1]),),s.U,s.A,s.SP,s.model,s.dis)
    end 
    state[nstep+1]   = State{OX+1,OX+1,1}(time[nstep]+Δt,(deepcopy(s.Λ[1]),),(deepcopy(s.X[1]),),s.U,s.A,s.SP,s.model,s.dis)
    buffer           = ntuple(i->𝕣1(undef,nXdof),10)  
    Lλx              = Vector{LU𝕣}(undef,nstep) 
    ΔXₐ              = [𝕣2(undef,nXdof,nAdof) for istep=1:nstep]     
    ΔX               = [𝕣1(undef,nXdof      ) for istep=1:nstep]   
    Lax              = [𝕣2(undef,nAdof,nXdof) for istep=1:nstep]
    Lxx              = Vector{Sparse𝕣2}(undef,nstep)
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
    # warming up: std Newmark-β
    for istep        = 1:nstep
        for iXiter   = 1:maxXiter
            firstXiter = iXiter==1 
            if    firstXiter assemble!{:step}(outX,asmX,dis,model,state[istep],Δt,(dbg...,solver=:SweepXA,phase=:warmup,step=step,iXiter=iXiter))
            else             assemble!{:iter}(outX,asmX,dis,model,state[istep],Δt,(dbg...,solver=:SweepXA,phase=:warmup,step=step,iXiter=iXiter))
            end
            try if  firstXiter Lλx[istep] = lu(outX.Lλx) 
            else               lu!(Lλx[istep], outX.Lλx) 
            end catch;         muscadeerror(@sprintf("Lλx matrix factorization failed at warm-up istep=%i, iXiter=%i",istep,iXiter)) end
            ΔX[1]      .= Lλx[istep]\-outX.Lλ
            Newmarkβincrement!{OX}(state[istep],ΔX[1] ,Xdofgr,outX.c,firstXiter,buffer...) 
            δx²,Lλ²     = sum(ΔX[1].^2),sum(outX.Lλ.^2)
            cXiter     += 1
            if δx²≤cΔx² && Lλ²≤cLλ² 
                break#out of the iXiter loop
            end
            iXiter==maxXiter && muscadeerror(@sprintf("no X-convergence at warm-up, istep=%3d after %3d X-iterations |ΔX|=%g / %g, |Lλ|=%g / %g",istep,iXiter,√(δx²),maxΔx,√(Lλ²)^2,maxLλ))
        end
    end

    A =  state[1].A

    # main part
    for iAiter = 1:maxAiter
        assembleA!{:Acost}(outXA,asmXA,dis,model,state[0],(dbg...,solver=:SweepXA,phase=:Acost,iAiter=iAiter))
        La♯             .= outXA.La   
        Laa♯            .= outXA.Laa  

        # forward sweep
        for istep        = 1:nstep
            Λ             =  state[istep].Λ 
            X             = (state[istep].X...,state[istep-1].X...) 
            U             =  state[istep].U 
            t             =  state[istep].time            
            assemble!{:Xsweep}(outXA,asmXA,dis,model,Λ,X,U,A,t,Δt,(dbg...,solver=:SweepXA,phase=:sensitivity,iAiter=iAiter,step=step))
            try if iAiter==1  Lλx[istep] = lu(outXA.Lλx) 
            else              lu!(Lλx[istep], outXA.Lλx) 
            end catch;        muscadeerror(@sprintf("Lλx matrix factorization failed at iAiter=%3d, istep=%i, iXiter=%i",iAiter,istep,iXiter)) end

            ΔX[ istep] .= Lλx[istep]\-outXA.Lλ  
            ΔXₐ[istep] .= Lλx[istep]\-outXA.Lλa 
 
            Lxx[istep]    =                     copy(outXA.Lxx) 
            Lax[istep]   .=                          outXA.Lax 
            LxxΔx        .=                          outXA.Lxx  ∘₁ ΔX[ istep] .+ outXA.Lr   # x
            LxxΔxₐ       .=                          outXA.Lxx  ∘₁ ΔXₐ[istep] .+ outXA.Lxr  # xa
            LxΔxₐ        .=                          outXA.Lx   ∘₁ ΔXₐ[istep] .+ outXA.Lr   # a
            LaxΔxₐ       .=                          outXA.Lax  ∘₁ ΔXₐ[istep] .+ outXA.Lar  # aa 
            LaxΔx        .=                          outXA.Lax  ∘₁ ΔX[ istep] .+ outXA.Lar  # a
            ΔxₐLxxΔx     .= ΔXₐ[istep]' ∘₁ LxxΔx  .+ outXA.Lxr' ∘₁ ΔX[ istep] .+ outXA.Lrr  # a
            ΔxₐLxxΔxₐ    .= ΔXₐ[istep]' ∘₁ LxxΔxₐ .+ outXA.Lxr' ∘₁ ΔXₐ[istep] .+ outXA.Lrr  # aa  # is symmetric
            La♯         .+= ΔxₐLxxΔx  .+ LaxΔx  .+ LxΔxₐ     + outXA.La                     # a
            Laa♯        .+= ΔxₐLxxΔxₐ .+ LaxΔxₐ .+ LaxΔxₐ'   + outXA.Laa                    # aa   
        end # istep

        # update A
        ΔA              .= Laa♯\-La♯
        ΔA²,La²          = sum(ΔA.^2),sum(La♯.^2)
        verbose && @printf "    In A-iteration %3d, |ΔA|=%7.1e |La♯|=%7.1e\n" iAiter √(ΔA²) √(La²)

        # backward sweep 
        for istep = nstep:-1:1
            Λ             = (state[istep].Λ...,state[istep+1].Λ...) 
            X             = (state[istep].X...,state[istep+1].X...) 
            U             = (state[istep].U...,state[istep+1].U...) 
            t             = state[istep].time            
            assemble!{:Λsweep}(outXA,asmXA,dis,model,Λ,X,U,A,t,Δt,(dbg...,solver=:SweepXA,phase=:backward,iAiter=iAiter,step=step)) # need Lx, Lr backwards
            ΔX[istep]   .+= ΔXₐ[istep] ∘₁ ΔA         
            LxxΔx        .=                    Lxx[istep]  ∘₁ ΔX[istep] .+ outXA.Lr   
            Lx♯          .= outXA.Lx + LxxΔx + Lax[istep]' ∘₁ ΔA 
            ΔΛ           .= Lλx[istep]'\-Lx♯  
        end

        # updates
        for istep = 1:nstep
            Newmarkβsweepincrement!{OX}(state[istep],state[istep-1],ΔX[istep],Xdofgr,outX.c,buffer...) 
            Newmarkβsweepincrement!{OX}(state[istep],state[istep+1],ΔΛ       ,Λdofgr,outX.c,buffer...) 
        end
        increment!(state[1],1,ΔA,Adofgr) # ∀ i,j  state[i].A === state[j].A

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