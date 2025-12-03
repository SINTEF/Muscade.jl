### Assembler

mutable struct AssemblySweepXA{OX,NDX} <: Assembly
    # up
    La        :: 𝕣1  
    Lλa       :: 𝕣2 
    Laa       :: 𝕣2 
    # down
    c         :: Newmarkβcoefficients{OX}
end   

function prepare(::Type{AssemblySweepXA{OX}},model,dis) where{OX}
    Xdofgr             = allXdofs(model,dis) 
    Adofgr             = allAdofs(model,dis)
    nXdof  = nΛdof     = getndof(Xdofgr)
    nAdof              = getndof(Adofgr)
    narray,neletyp     = 3,getneletyp(model)
    asm                = Matrix{𝕫2}(undef,narray,neletyp)  # asm[iarray,ieletyp][ieledof,iele]
    asmx               = Vector{𝕫2}(undef       ,neletyp)  # asmx[ieletyp][ieledof,iele]
    Lx                 = asmvec!(asmx          ,Xdofgr,dis)
    La                 = asmvec!(view(asm, 1,:),Adofgr,dis)
    Lλa                = asmfullmat!(view(asm, 2,:),asmx,view(asm,1,:),nΛdof,nAdof)  
    Laa                = asmfullmat!(view(asm, 3,:),view(asm,1,:),view(asm,1,:),nAdof,nAdof)

    out                = AssemblySweepXA{OX,OX+1}(La,Lλa,Laa, Newmarkβcoefficients{OX}()) 
    return out,asm,Xdofgr,Adofgr
end
function zero!(out::AssemblySweepXA) 
    zero!(out.La )
    zero!(out.Lλa)
    zero!(out.Laa)
end

#=        TODO
solver
write specific adiff for ElementCost
SweepXA for order 0 and 1
Multi load cases        
=#
function Newmarkβdecrement!{OX}(Xₐ,ΔXₐ,Xgr,Agr,c) where{OX}
    # xₐ    -=    Δxₐ
    # xₐ′   -= a₁*Δxₐ
    # xₐ″   -= b₁*Δxₐ
    f  = (1.,c.a₁,c.b₁)  
    nA = size(Xₐ[1],2)
    for ider = 1:OX+1
        for igA ∈ eachindex(Agr.iA)
            iA,jA,sA = Agr.iA[igA], Agr.jA[igA], 1 / Agr.scaleA[igA] # inverse scaleA, because Xₐ is dX / dA 
            for igX ∈ eachindex(Xgr.iX)
                iX,jX,sX = Xgr.iX[igX], Xgr.jX[igX], Xgr.scaleX[igX]
                Xₐ[ider][iX,iA] -= f[ider] .* ΔXₐ[jX,jA] * sX * sA
            end
        end 
    end
end


# # Ignore no_second_order - we will only use it with length(δr)==1, promise.
# function getresidual_with_higher_order(eleobj::Eleobj,X,U,A,t,SP,dbg,req) where{Eleobj} 
#     R,FB,eleres = getresidual(eleobj,hasresidual(eleobj),haslagrangian(eleobj),Val(false),X,U,A,t,SP,dbg,req) 
#     hasnan(R,FB) && muscadeerror((dbg...,t=t,SP=SP ),@sprintf("residual(%s,...) returned NaN in R, FB or derivatives",Eleobj))  
#     return R,FB,eleres
# end
# function getresidual_with_higher_order(eleobj::Eleobj,X,U,A,t,SP,dbg,req...) where{Eleobj} 
#     R,FB     = getresidual(eleobj,hasresidual(eleobj),haslagrangian(eleobj),Val(false),X,U,A,t,SP,dbg    ) 
#     hasnan(R,FB) && muscadeerror((dbg...,t=t,SP=SP ),@sprintf("residual(%s,...) returned NaN in R, FB or derivatives",Eleobj))  
#     return R,FB
# end

#### special case for the computation of Rₐ in SweepXA, add a SMatrix into a Matrix.  No 'build' or 'split'.
# function add_value!(out::𝕣0,a,ia::𝕫;Δt=idmult)  # Lr, scalar in Newmark-β context
#     out[] += VALUE(a[ia])*Δt
# end
function add_value!(out::𝕣2,asm,iele,a::SMatrix ; Δt=idmult) 
    for i ∈ eachindex(a)
        iout = asm[i,iele]
        if iout≠0 
            out[iout]+=VALUE(a[i])*Δt 
        end
    end
end   
add_∂!{         P,S,T}(out::𝕣2,asm, iele, a::SMatrix{Na,Ma,R        } ; Δt=idmult) where{P,S,T,Na,Ma,R} = nothing # if a is 𝕣2 or P does not match a
function add_∂!{P,S,T}(out::𝕣2,asm, iele, a::SMatrix{Na,Ma,∂ℝ{P,1,R}} ; Δt=idmult) where{P,S,T,R,Na,Ma} 
    for i ∈ eachindex(a)
        iout = asm[i,iele]
        if iout≠0
            if     S==:plus   out[iout]+=a[i].dx[1]*Δt  
            elseif S==:minus  out[iout]-=a[i].dx[1]*Δt  
            else   muscadeerror((;S=S),"Illegal value of parameter S")    
            end
        end
    end
end   

abstract type WithXₐ               end
abstract type MissionQₐ♯ <: WithXₐ end
abstract type MissionXₐ  <: WithXₐ end
function addin!{MissionXₐ}(out::AssemblySweepXA{OX},asm,iele,scale,eleobj,Λ,X::NTuple{NXder,<:SVector{Nx}},U,A::SVector{Na},Xₐ::NTuple{NXder,<:SMatrix{Nx,Na}},t,Δt,SP,dbg) where{OX,NXder,Nx,Na}
    @assert NXder == OX+1                                                    
    δA                = δ{1,Na,𝕣}(scale.A)
    vX                = ntuple(ider->X[ider] + Xₐ[ider] ∘₁ δA, NXder)
    vA                =              A       +             δA
    R,FB              = getresidual(eleobj,vX,U,vA,t,SP,dbg) 
    add_∂!{1}( out.Lλa ,asm[2],iele,R)  
end
function addin!{MissionQₐ♯}(out::AssemblySweepXA{OX},asm,iele,scale,eleobj,Λ,X::NTuple{NXder,<:SVector{Nx}},U,A::SVector{Na},Xₐ::NTuple{NXder,<:SMatrix{Nx,Na}},t,Δt,SP,dbg) where{OX,NXder,Nx,Na}
    @assert NXder == OX+1
    λ                 = ∂0(Λ)
    (δA,)             = reδ{2}((;A),(;A=scale.A)) 
    vX                = ntuple(ider->X[ider] + Xₐ[ider] ∘₁ δA, NXder)
    vA                =              A       +             δA
    L,FB              = getlagrangian(eleobj,λ,vX,U,vA,t,SP,dbg) # TODO jump over elements with residual.  "getcost"
    # REPRISE L is a pack of zeros. inputs vX, vA are good
    ∇L                = ∂{2,Na}(L)
    add_value!(out.La , asm[1], iele, ∇L ; Δt)             
    add_∂!{1 }(out.Laa, asm[3], iele, ∇L ; Δt)  
    if dbg.istep==10 && dbg.ieletyp==1
        @show dbg
        @show typeof(eleobj)
        @show Xₐ
        @show vX
        @show vA
        @show L
    end
end

function addin!{:Acost}(out::AssemblySweepXA,asm,iele,scale,eleobj::Acost,A::SVector{Na},dbg) where{Na} 
    d      = revariate{2}((;A),(;A=scale.A)) # careful: revariate returns a NamedTuple
    ø      = nothing
    C,_    = lagrangian(eleobj,ø,ø,ø,d.A,ø,ø ,dbg)
    ∇ₐC    = ∂{2,Na}(C)
    add_value!(out.La ,asm[1],iele,∇ₐC)
    add_∂!{1 }(out.Laa,asm[3],iele,∇ₐC)
end
# special assembly, passing and disassembling Xₐ
function assemble!{mission}(out::Assembly,asm,dis,model,state,Δt,Xₐ,dbg) where{mission <: WithXₐ }
    zero!(out)
    for ieletyp = 1:lastindex(model.eleobj)
        eleobj  = model.eleobj[ieletyp]
        assemble_!{mission}(out,view(asm,:,ieletyp),dis.dis[ieletyp],eleobj,state,state.time,Δt,state.SP,Xₐ,(dbg...,ieletyp=ieletyp))
    end
end
assemble_!{         mission}(out::Assembly,asm,dis,eleobj::Vector{<:Acost},state::State{nΛder,nXder,nUder},t,Δt,SP,Xₐ,dbg) where{mission <: WithXₐ,nΛder,nXder,nUder} = nothing
function assemble_!{mission}(out::Assembly,asm,dis,eleobj                 ,state::State{nΛder,nXder,nUder},t,Δt,SP,Xₐ,dbg) where{mission <: WithXₐ,nΛder,nXder,nUder}
    for iele  = 1:lastindex(eleobj)
        index = dis.index[iele]
        Λe    = NTuple{nΛder}(λ[index.X] for λ∈state.Λ)
        Xe    = NTuple{nXder}(x[index.X] for x∈state.X)
        Ue    = NTuple{nUder}(u[index.U] for u∈state.U)
        Ae    = state.A[index.A]
        Xₐe   = NTuple{nXder}(xₐ[index.X,index.A] for xₐ∈Xₐ)
        addin!{mission}(out,asm,iele,dis.scale,eleobj[iele],Λe,Xe,Ue,Ae,Xₐe, t,Δt,SP,(dbg...,iele=iele)) # defined by solver.  Called for each element. But the asm that is passed
    end                                                                              # is of the form asm[iarray][i,iele], because addin! will add to all arrays in one pass
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
- `maxXiter=50`        maximum number of equilibrium iterations at each step.
- `maxΔx=1e-5`        convergence criteria: norm of `X`. 
- `maxLλ=∞`           convergence criteria: norm of the residual. 
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
                    time::AbstractVector{𝕣},
                    initialstate::State,
                    β::𝕣=1/4,γ::𝕣=1/2,
                    maxXiter::ℤ=50,maxΔx::ℝ=1e-5,maxLλ::ℝ=∞,
                    maxAiter::ℤ=50,maxΔa::ℝ=1e-5,maxLa::ℝ=∞) where{OX}

    model,dis        = initialstate.model,initialstate.dis
    outX ,asmX ,Xdofgr          = prepare(AssemblySweepX{ OX},model,dis)  
    outXA,asmXA,Xdofgr,Adofgr   = prepare(AssemblySweepXA{OX},model,dis)  
    nXdof            = getndof(Xdofgr)
    nAdof            = getndof(Adofgr)
    buffer           = ntuple(i->𝕣1(undef,nXdof), 6)  
    cΔX²,cLλ²        = maxΔx^2,maxLλ^2
    cΔA²,cLa²        = maxΔa^2,maxLa^2
    cXiter           = 0

    state            = State{1,OX+1,1}(copy(initialstate)) 
    state.Λ[1]      .= 0.
    states           = allocate(pstate,Vector{typeof(state)}(undef,length(time))) 

    La♯              = 𝕣1(undef,nAdof      )
    Laa♯             = 𝕣2(undef,nAdof,nAdof)
    Xₐ               = ntuple(i->𝕣2(undef,nXdof,nAdof),OX+1)
    local Lλx # declare Lλx to scope the function, without having to actualy initialize the variable

    Ra = 𝕣1(undef,length(time))       ### dbg
    deltaXa = 𝕣1(undef,length(time))  ### dbg
    Xa = 𝕣1(undef,length(time))       ### dbg
    Va = 𝕣1(undef,length(time))       ### dbg
    Aa = 𝕣1(undef,length(time))       ### dbg
    ΔXs = 𝕣1(undef,length(time))      ### dbg

    # main part
    for iAiter = 1:maxAiter

        # time independant Acost
        assembleA!{:Acost}(outXA,asmXA,dis,model,state,(dbg...,solver=:SweepXA,phase=:Acost,iAiter=iAiter))
        La♯              .= outXA.La   
        Laa♯             .= outXA.Laa  
        zero!(Xₐ)

        for (istep,t)    ∈ enumerate(time)
            oldt         = state.time
            state.time   = t
            Δt           = t-oldt
            Δt ≤ 0 && OX>0 && muscadeerror(@sprintf("Time step length not strictly positive at step=%3d",istep))
            outX.c = outXA.c = Newmarkβcoefficients{OX}(Δt,β,γ)

            # step and iterations

            for iXiter   = 1:maxXiter
                cXiter  += 1
                firstiter = iXiter==1
                if   firstiter assemble!{:step}(outX,asmX,dis,model,state,Δt,(dbg...,solver=:SweepXA,mission=:step,iAiter=iAiter,istep=istep,iXiter=iXiter))
                else           assemble!{:iter}(outX,asmX,dis,model,state,Δt,(dbg...,solver=:SweepXA,mission=:iter,iAiter=iAiter,istep=istep,iXiter=iXiter))
                end
                try if istep==1  && firstiter  Lλx = lu(outX.Lλx) # here we do not write "local Lλx", so we refer to the variable defined outside the loops (we do not shadow Lλx)
                else                           lu!(Lλx, outX.Lλx) 
                end catch;    muscadeerror(@sprintf("matrix factorization failed at iAiter=%i step=%i, iXiter=%i",iAiter,istep,iXiter)) end
                ΔX       = Lλx\outX.Lλ

                if iXiter == 1
                    ΔXs[istep] = ΔX[1] ### dbg
                end

                ΔX²,Lλ²  = sum(ΔX.^2),sum(outX.Lλ.^2)
                Newmarkβdecrement!{OX}(state,ΔX ,Xdofgr,outX.c,firstiter,buffer...)

                if ΔX²≤cΔX² && Lλ²≤cLλ² 
                    #verbose && @printf "        step %3d converged in %3d iterations. |Δx|=%7.1e |Lλ|=%7.1e\n" istep iXiter √(ΔX²) √(Lλ²)
                    states[istep] = State(state.time,state.Λ,deepcopy(state.X),state.U,state.A,state.SP,model,dis)
                    break#out of the iXiter loop
                end
                iXiter==maxXiter && muscadeerror(@sprintf("no convergence of step %3d for iAiter %3d after %3d iterations |Δx|=%g / %g, |Lλ|=%g / %g",istep,iAiter,iXiter,√(ΔX²),maxΔx,√(Lλ²)^2,maxLλ))
            end

            # sensitivity

            a₁,a₂,a₃,b₁,b₂,b₃ = outXA.c.a₁,outXA.c.a₂,outXA.c.a₃,outXA.c.b₁,outXA.c.b₂,outXA.c.b₃
            if OX≥2 b₂♯,b₃♯   = b₂/(1-a₂), a₃/(1-a₂)+b₃       end
            if OX≥1 Xₐ[2]   .-= a₂  .* Xₐ[2] .+ a₃ .* Xₐ[3]   end #         xₐ′-= aₐ
            if OX≥2 Xₐ[3]   .-= b₂♯ .* Xₐ[2] .+ b₃♯ .*Xₐ[3]   end # same as xₐ″-= bₐ but in place
            assemble!{MissionXₐ}(outXA,asmXA,dis,model,state,Δt,Xₐ,(dbg...,solver=:SweepXA,mission=:Xₐ,iAiter=iAiter,istep=istep))
            ΔXₐ       = Lλx\outXA.Lλa 
            Newmarkβdecrement!{OX}(Xₐ,ΔXₐ ,Xdofgr,Adofgr,outXA.c)

            Ra[istep] = outXA.Lλa[1,1]  ### dbg   

            assemble!{MissionQₐ♯}(outXA,asmXA,dis,model,state,Δt,Xₐ,(dbg...,solver=:SweepXA,mission=:Qₐ♯,iAiter=iAiter,istep=istep))
            La♯         .+= outXA.La    
            Laa♯        .+= outXA.Laa 
            

             
            deltaXa[istep] = -ΔXₐ[1,1]  ### dbg   
            Xa[istep] = Xₐ[1][1]  ### dbg   
            Va[istep] = Xₐ[2][1]  ### dbg   
            Aa[istep] = Xₐ[3][1]  ### dbg   

        end # istep

        Xs = [s.X[1][1] for s∈states] ### dbg   
        Vs = [s.X[2][1] for s∈states] ### dbg   

        if iAiter==1 ### dbg   
            fig      = Figure(size = (1000,800))
            axeX      = Axis(fig[1,1])
            lines!(axeX,time,Xs,color=:black)
           #lines!(axeX,time,Vs,color=:grey)
           #lines!(axeX,time,Ra,color=:blue)
           #lines!(axeX,time,120*deltaXa,color=:green)
           lines!(axeX,time,Xs+ .1*Xa,color=:red)
          # lines!(axeX,time,-10*Va,color=:magenta)
           #lines!(axeX,time,-10*Aa,color=:cyan)
           #lines!(axeX,time,ΔXs*-10,color=:orange)
           display(fig)
        end

 
        # update A
        ΔA               = Laa♯\La♯  # TODO try catch
        @show Laa♯[1,1],La♯[1],ΔA[1]  ### dbg
        ΔA²,La²          = sum(ΔA.^2),sum(La♯.^2)
        decrement!(state,1,ΔA,Adofgr) 
        verbose && @printf "    In A-iteration %3d, |ΔA|=%7.1e |La♯|=%7.1e\n" iAiter √(ΔA²) √(La²)

        # Aiter convergence
        if ΔA²≤cΔA² && La²≤cLa² 
            verbose && @printf "    SweepXA converged in %3d A-iterations. |ΔA|=%7.1e / %g |La|=%7.1e / %g\n" iAiter √(ΔA²) maxΔa √(La²) maxLa
            break#out of the iAiter loop
        end
        iAiter==maxAiter && muscadeerror(@sprintf("no convergence of SweepXA after %3d A-iterations |ΔA|=%g / %g, |La|=%g / %g",iAiter,√(ΔA²),maxΔa,√(La²)^2,maxLa))

        # reset state to initial conditions
        state.time = initialstate.time
        for i=1:min(OX+1,length(initialstate.X))
            state.X[i]     .= initialstate.X[i]
        end
        for i= length(initialstate.X)+1:OX+1
            state.X[i]     .= 0.
        end
    end 
    verbose && @printf "\n    nel=%d, nXdof=%d, nstep=%d, ΣnXiter=%d, mean(nXiter)=%d\n" getnele(model) getndof(Xdofgr) length(time) cXiter cXiter/length(time)

    return
end