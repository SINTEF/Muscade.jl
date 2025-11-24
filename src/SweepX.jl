### Assembler

struct Newmarkβcoefficients{OX}
    a₁::𝕣
    a₂::𝕣 
    a₃::𝕣
    b₁::𝕣
    b₂::𝕣
    b₃::𝕣
    Δt::𝕣    
end
Newmarkβcoefficients{0}(Δt,_,_)          = Newmarkβcoefficients{0}(0.      ,0.  ,0.         ,0.        ,0.      ,0.  ,Δt)
Newmarkβcoefficients{1}(Δt,_,γ)          = Newmarkβcoefficients{1}(1/(γ*Δt),1/γ ,0.         ,0.        ,0.      ,0.  ,Δt)
Newmarkβcoefficients{2}(Δt,β,γ)          = Newmarkβcoefficients{2}(γ/(β*Δt),γ/β ,(γ/2β-1)*Δt,1/(β*Δt^2),1/(β*Δt),1/2β,Δt)
Newmarkβcoefficients{O}(      ) where{O} = Newmarkβcoefficients{O}(0.      ,0.  ,0.         ,0.        ,0.      ,0.  ,0.)

mutable struct AssemblySweepX{OX,Tλ,Tλx} <: Assembly
    # up
    Lλ        :: Tλ                
    Lλx       :: Tλx
    # down
    c         :: Newmarkβcoefficients{OX}
end   
function prepare(::Type{AssemblySweepX{OX}},model,dis) where{OX}
    Xdofgr             = allXdofs(model,dis)  # dis: the model's disassembler
    ndof               = getndof(Xdofgr)
    narray,neletyp     = 2,getneletyp(model)
    asm                = Matrix{𝕫2}(undef,narray,neletyp)  # asm[iarray,ieletyp][ieledof,iele]
    Lλ                 = asmvec!(view(asm,1,:),Xdofgr,dis) 
    Lλx                = asmmat!(view(asm,2,:),view(asm,1,:),view(asm,1,:),ndof,ndof) 
    out                = AssemblySweepX{OX,typeof(Lλ),typeof(Lλx)}(Lλ,Lλx,Newmarkβcoefficients{OX}()) 
    return out,asm,Xdofgr
end
function zero!(out::AssemblySweepX) 
    zero!(out.Lλ)
    zero!(out.Lλx)
end
# jump over elements without Xdofs in a SweepX analysis, for all orders, all missions
addin!{:step}(out::AssemblySweepX{0},asm,iele,scale,eleobj,Λ,X::NTuple{Nxder,<:SVector{0}},U,A,t,Δt,SP,dbg) where{Nxder} = return
addin!{:iter}(out::AssemblySweepX{0},asm,iele,scale,eleobj,Λ,X::NTuple{Nxder,<:SVector{0}},U,A,t,Δt,SP,dbg) where{Nxder} = return
addin!{:step}(out::AssemblySweepX{1},asm,iele,scale,eleobj,Λ,X::NTuple{Nxder,<:SVector{0}},U,A,t,Δt,SP,dbg) where{Nxder} = return
addin!{:iter}(out::AssemblySweepX{1},asm,iele,scale,eleobj,Λ,X::NTuple{Nxder,<:SVector{0}},U,A,t,Δt,SP,dbg) where{Nxder} = return
addin!{:step}(out::AssemblySweepX{2},asm,iele,scale,eleobj,Λ,X::NTuple{Nxder,<:SVector{0}},U,A,t,Δt,SP,dbg) where{Nxder} = return
addin!{:iter}(out::AssemblySweepX{2},asm,iele,scale,eleobj,Λ,X::NTuple{Nxder,<:SVector{0}},U,A,t,Δt,SP,dbg) where{Nxder} = return
function addin!{:step}(out::AssemblySweepX{2},asm,iele,scale,eleobj,Λ,X::NTuple{Nxder,<:SVector{Nx}},U,A,t,Δt,SP,dbg) where{Nxder,Nx}
    a₁,a₂,a₃,b₁,b₂,b₃ = out.c.a₁,out.c.a₂,out.c.a₃,out.c.b₁,out.c.b₂,out.c.b₃
    x,x′,x″    = ∂0(X),∂1(X),∂2(X)
    δX,δr      = reδ{1}((;X=x,r=0.),(;X=scale.X,r=1.))
    a          = a₂*x′ + a₃*x″
    b          = b₂*x′ + b₃*x″
    vx         = x  +    δX
    vx′        = x′ + a₁*δX + a*δr 
    vx″        = x″ + b₁*δX + b*δr
    Lλ,FB      = getresidual(eleobj,(vx,vx′,vx″),U,A,t,SP,dbg)
    Lλ         = Lλ .* scale.X
    add_value!(       out.Lλ ,asm[1],iele,Lλ             )  # rhs  = R    
    add_∂!{1,:minus}( out.Lλ ,asm[1],iele,Lλ,1:Nx,(Nx+1,))  # rhs +=  -C⋅a -M⋅b 
    add_∂!{1       }( out.Lλx,asm[2],iele,Lλ,1:Nx,1:Nx   )  # Mat  =  K + a₁C + b₁M
end
function addin!{:iter}(out::AssemblySweepX{2},asm,iele,scale,eleobj,Λ,X::NTuple{Nxder,<:SVector{Nx}},U,A,t,Δt,SP,dbg) where{Nxder,Nx} 
    a₁,b₁      = out.c.a₁,out.c.b₁
    δX         = δ{1,Nx,𝕣}(scale.X)
    Lλ,FB      = getresidual(eleobj,(∂0(X)+δX, ∂1(X)+a₁*δX, ∂2(X)+b₁*δX),U,A,t,SP,dbg)
    Lλ         = Lλ .* scale.X
    add_value!(out.Lλ ,asm[1],iele,Lλ          )
    add_∂!{1}( out.Lλx,asm[2],iele,Lλ,1:Nx,1:Nx)
end
function addin!{:step}(out::AssemblySweepX{1},asm,iele,scale,eleobj,Λ,X::NTuple{Nxder,<:SVector{Nx}},U,A,t,Δt,SP,dbg) where{Nxder,Nx}
    a₁,a₂      = out.c.a₁,out.c.a₂
    x,x′       = ∂0(X),∂1(X)
    δX,δr      = reδ{1}((;X=x,r=0.),(;X=scale.X,r=1.))
    a          = a₂*x′
    vx         = x  +    δX   
    vx′        = x′ + a₁*δX + a*δr  
    Lλ,FB      = getresidual(eleobj,(vx,vx′),U,A,t,SP,dbg)
    Lλ         = Lλ .* scale.X
    add_value!(out.Lλ ,asm[1],iele,Lλ                    )  # rhs  = R    
    add_∂!{1,:minus}( out.Lλ ,asm[1],iele,Lλ,1:Nx,(Nx+1,))  # rhs +=  -C⋅a 
    add_∂!{1}( out.Lλx,asm[2],iele,Lλ,1:Nx,1:Nx          )  # Mat  = K + C/Δt 
end
function addin!{:iter}(out::AssemblySweepX{1},asm,iele,scale,eleobj,Λ,X::NTuple{Nxder,<:SVector{Nx}},U,A,t,Δt,SP,dbg) where{Nxder,Nx}
    a₁         = out.c.a₁
    δX         = δ{1,Nx,𝕣}(scale.X)
    Lλ,FB      = getresidual(eleobj,(∂0(X)+δX, ∂1(X)+a₁*δX),U,A,t,SP,dbg)
    Lλ         = Lλ .* scale.X
    add_value!(out.Lλ ,asm[1],iele,Lλ           )  # rhs  = R    
    add_∂!{1}( out.Lλx,asm[2],iele,Lλ,1:Nx,1:Nx )  # Mat  = K + C/Δt 
end

function addin!{Both}(out::AssemblySweepX{0},asm,iele,scale,eleobj,Λ,X::NTuple{Nxder,<:SVector{Nx}},U,A,t,Δt,SP,dbg) where{Both,Nxder,Nx} 
    δX         = δ{1,Nx,𝕣}(scale.X)
    Lλ,FB      = getresidual(eleobj,(∂0(X)+δX,),U,A,t,SP,dbg)
    Lλ         = Lλ .* scale.X
    add_value!(out.Lλ ,asm[1],iele,Lλ)
    add_∂!{1}( out.Lλx,asm[2],iele,Lλ)
end

struct   Newmarkβdecrement!{OX} end
function Newmarkβdecrement!{2}(state,Δx ,Xdofgr,c,firstiter, a,b,x′,x″,Δx′,Δx″,args...) # x′, x″ are just mutable memory, neither input nor output.
    a₁,a₂,a₃,b₁,b₂,b₃ = c.a₁,c.a₂,c.a₃,c.b₁,c.b₂,c.b₃

    if firstiter
        getdof!(state,1,x′,Xdofgr) 
        getdof!(state,2,x″,Xdofgr) 
        a       .= a₂*x′.+ a₃*x″ 
        b       .= b₂*x′.+ b₃*x″
        Δx′     .= a₁*Δx .+ a
        Δx″     .= b₁*Δx .+ b
    else
        Δx′     .= a₁*Δx 
        Δx″     .= b₁*Δx 
    end
    decrement!(state,1,Δx ,Xdofgr)
    decrement!(state,2,Δx′,Xdofgr)
    decrement!(state,3,Δx″,Xdofgr)
end
function Newmarkβdecrement!{1}(state,Δx ,Xdofgr,c,firstiter, a,x′,Δx′,args...)
    a₁,a₂ = c.a₁,c.a₂

    if firstiter
        getdof!(state,1,x′,Xdofgr) 
        a       .= a₂*x′
        Δx′     .= a₁*Δx .+ a
    else
        Δx′     .= a₁*Δx 
    end
    decrement!(state,1,Δx ,Xdofgr)
    decrement!(state,2,Δx′,Xdofgr)
end
function Newmarkβdecrement!{0}(state,Δx ,Xdofgr,args...)
    decrement!(state,1,Δx ,Xdofgr)
end

"""
	SweepX{OX}

A non-linear, time domain solver, that solves the problem time-step by time-step.
Only the `X`-dofs of the model are solved for, while `U`-dofs and `A`-dofs are unchanged.

- `SweepX{0}` is Newton-Raphson. 
- `SweepX{1}` is a first order variant of Newmark-β with Newton-Raphson iterations. 
- `SweepX{2}` is Newmark-β, with Newton-Raphson iterations.

IMPORTANT NOTE: Muscade does not allow elements to have state variables, for example, plastic strain,
or shear-free position for dry friction.  Where the element implements such physics, this 
is implemented by introducing the state as a degree of freedom of the element, and solving
for its evolution, *even in a quasi-static problem*, requires the use of `OX≥1`.

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
- `β=1/4`,`γ=1/2`     parameters to the Newmark-β algorithm. 
                      `β` is dummy if `OX<2`.
                      `γ` is dummy if `OX<1`.
- `maxiter=50`        maximum number of equilibrium iterations at each step.
- `maxΔx=1e-5`        convergence criteria: norm of `X`. 
- `maxLλ=∞`           convergence criteria: norm of the residual. 
- `saveiter=false`    set to true so that output `states` contains the state
                      at the iteration of the last step analysed.  Useful to study
                      a step that fails to converge. 

# Output

A vector of length equal to that of the named input argument `time` containing the states at the time steps.                       

See also: [`solve`](@ref), [`initialize!`](@ref), [`findlastassigned`](@ref), [`study_singular`](@ref), [`DirectXUA`](@ref), [`FreqXU`](@ref)  
"""
struct        SweepX{OX} <: AbstractSolver end
function solve(SX::Type{SweepX{OX}},pstate,verbose,dbg;
                    time::AbstractVector{𝕣},
                    initialstate::State,
                    β::𝕣=1/4,γ::𝕣=1/2,
                    maxiter::ℤ=50,maxΔx::ℝ=1e-5,maxLλ::ℝ=∞,
                    saveiter::𝔹=false) where{OX}
                    
    model,dis        = initialstate.model,initialstate.dis
    out,asm,Xdofgr   = prepare(AssemblySweepX{OX},model,dis)  
    nXdof            = getndof(Xdofgr)
    buffer           = ntuple(i->𝕣1(undef,nXdof), 6)  
    citer            = 0
    cΔx²,cLλ²        = maxΔx^2,maxLλ^2
    state            = State{1,OX+1,1}(copy(initialstate)) 
    states           = allocate(pstate,Vector{typeof(state)}(undef,saveiter ? maxiter : length(time))) # states is not a return argument of this function.  Hence it is not lost in case of exception
    local Lλx # declare Lλx to scope the function, without having to actualy initialize the variable
    for (step,t)     ∈ enumerate(time)
        oldt         = state.time
        state.time   = t
        Δt           = t-oldt
        Δt ≤ 0 && OX>0 && muscadeerror(@sprintf("Time step length not strictly positive at step=%3d",step))
        out.c        = Newmarkβcoefficients{OX}(Δt,β,γ)
        for iiter    = 1:maxiter
            citer   += 1
            firstiter = iiter==1
            if   firstiter assemble!{:step}(out,asm,dis,model,state,Δt,(dbg...,solver=:SweepX,step=step,iiter=iiter))
            else           assemble!{:iter}(out,asm,dis,model,state,Δt,(dbg...,solver=:SweepX,step=step,iiter=iiter))
            end
            try if step==1  && firstiter  Lλx = lu(out.Lλx) # here we do not write "local Lλx", so we refer to the variable defined outside the loops (we do not shadow Lλx)
            else                          lu!(Lλx, out.Lλx) 
            end catch;    muscadeerror(@sprintf("matrix factorization failed at step=%i, iiter=%i",step,iiter)) end
            Δx       = Lλx\out.Lλ
            Δx²,Lλ²  = sum(Δx.^2),sum(out.Lλ.^2)
            Newmarkβdecrement!{OX}(state,Δx ,Xdofgr,out.c,firstiter,buffer...)
 
            verbose && saveiter && @printf("        iteration %3d, γ= %7.1e\n",iiter,γ)
            saveiter && (states[iiter]=State(state.time,state.Λ,deepcopy(state.X),state.U,state.A,state.SP,model,dis))
            if Δx²≤cΔx² && Lλ²≤cLλ² 
                verbose && @printf "    step %3d converged in %3d iterations. |Δx|=%7.1e |Lλ|=%7.1e\n" step iiter √(Δx²) √(Lλ²)
                ~saveiter && (states[step]=State(state.time,state.Λ,deepcopy(state.X),state.U,state.A,state.SP,model,dis))
                break#out of the iiter loop
            end
            iiter==maxiter && muscadeerror(@sprintf("no convergence of step %3d after %3d iterations |Δx|=%g / %g, |Lλ|=%g / %g",step,iiter,√(Δx²),maxΔx,√(Lλ²)^2,maxLλ))
        end
    end
    verbose && @printf "\n    nel=%d, ndof=%d, nstep=%d, niter=%d, niter/nstep=%5.2f\n" getnele(model) getndof(Xdofgr) length(time) citer citer/length(time)
    return
end
