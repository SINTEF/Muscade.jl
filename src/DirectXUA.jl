const λxua = 1:4
const sym  = (λ=1,x=2,u=3,a=4)

# We make a distinction between nΛder==nAder==1, nXder=length(X), nUder=length(U) on the one hand, and mΞder ≤ nΞder.  This allows
# 1) to freeze A for XU algo (or any class)
# 2) not to compute cost on U′ or U′′ if these costs are known to be zero (same with X)                                      

mutable struct AssemblyDirect{T1,T2}  <:Assembly
    L1    :: T1
    L2    :: T2
end  
struct AssemblerDirect{Mder}
    vec :: Matrix{𝕫2}
    mat :: Matrix{𝕫2}
end
function prepare(::Type{AssemblyDirect},model,dis,mder) 
    dofgr              = (allΛdofs(model,dis),allXdofs(model,dis),allUdofs(model,dis),allAdofs(model,dis))
    ndof               = getndof.(dofgr)
    neletyp            = getneletyp(model)
    vec                = Matrix{𝕫2}(undef,4,neletyp)
    mat                = Matrix{𝕫2}(undef,16,neletyp)
    asm                = AssemblerDirect{:full,mder}(vec,mat)
    L1                 = [asmvec!(view(asm.vec,α  ,:),dofgr[α],dis)                                 for ider=1:mder[α]               ] 
    L2                 = [asmmat!(view(asm.mat,α,β,:),view(asm,α,:),view(asm,β,:),ndof[α],ndof[β])  for ider=1:mder[α],jder=1:mder[β]]
    out                = AssemblyDirect(L1,L2)
    return out,asm#,Ydofgr,Adofgr
end
function zero!(out::AssemblyDirect)
    for α∈λxua 
        zero!.(out.L1[α])
        for β∈λxua
            zero!.(out.L2[α,β])
        end
    end
end
function addin!(out::AssemblyDirect,asm::AssemblerDirect{Mder},iele,scale,eleobj,Λ::SVector{Nx},X::NTuple{nXder,SVector{Nx,T}},
                                             U::NTuple{nUder,SVector{Nu,T}},A::SVector{Na},t,SP,dbg) where{nXder,nUder,Nx,Nu,Na,Mder,T} 

    ndof  = (Nx  ,Nx   ,Nu   ,Na  )

    nder  = (1   ,nXder,nUder,1   )
    V     = ((Λ,),X    ,U    ,(A,)) # does this trigger copying?
    p     = 0
    V∂    = ntuple(4) do α
                ntuple(nder[α]) do ider 
                    X∂ᵢ = ider>Mder[α] ? V[α][ider] : SVector{Nx}(  ∂²ℝ{1,Nz}(V[α][ider][idof],p+ix)   for idof=1:ndof[α]) # type stable?
                    p  += Nx
                    X∂ᵢ
                end
            end
    
    L,FB    = getlagrangian(eleobj, V∂[1][1],V∂[2],V∂[3],V∂[4][1],t,SP,dbg)
 
    ∇L      = ∂{2,Nz}(L)
    pα      = 0
    for α∈λxua, i=1:Mder[α]
        iα  = pα+(1:ndof[α])
        pα += ndof[α]
        add_value!(out.L1[α][i] ,asm.vec[α],iele,∇L,iα)
        pβ      = 0
        for β∈λxua, j=1:Mder[β]
            iβ  = pβ+(1:ndof[β])
            pβ += ndof[β]
            add_∂!{1}( out.L2[α,β][i,j],asm.mat[α,β],iele,∇L,iα,iβ)
        end
    end
end

######################

mutable struct AssemblyDirectLine  <:Assembly
    ming  :: 𝕣
    minλ  :: 𝕣
    Σλg   :: 𝕣
    npos  :: 𝕫
end  
struct AssemblerDirectLine end
prepare(::Type{AssemblyDirectLine} = AssemblyDirectLine(La,∞,∞,0.,0),AssemblerDirectLine()
function zero!(out::AssemblyDirectLine)
    out.ming = ∞    
    out.minλ = ∞
    out.Σλg  = 0.
    out.npos = 0    
end
function addin!(out::AssemblyDirectLine,asm::AssemblerDirectLine,iele,scale,eleobj,Λ,X,U,A,t,SP,dbg) 
    L,FB    = getlagrangian(eleobj, Λ,X,U,A,t,SP,dbg)
    if hasfield(typeof(FB),:mode) && FB.mode==:positive
        out.ming   = min(out.ming,FB.g)
        out.minλ   = min(out.minλ,FB.λ)
        out.Σλg   += FB.g*FB.λ
        out.npos  += 1
    end
end



"""
	DirectXUA

A non-linear direct solver for optimisation FEM.

An analysis is carried out by a call with the following syntax:

```
initialstate    = initialize!(model)
setdof!(initialstate,1.;class=:U,field=:λcsr)
stateX          = solve(SweepX{0}  ;initialstate=initialstate,time=[0.,1.])
stateXUA        = solve(DirectXUA;initialstate=stateX)
```

The interior point algorithm requires a starting point that is
strictly primal feasible (at all steps, all inequality constraints must have 
positive gaps) and strictly dual feasible (at all steps, all associated Lagrange 
multipliers must be strictly positive). Note the use of `setdof!` in the example
above to ensure dual feasibility.

# Named arguments
- `dbg=(;)`           a named tuple to trace the call tree (for debugging)
- `verbose=true`      set to false to suppress printed output (for testing)
- `silenterror=false` set to true to suppress print out of error (for testing) 
- `initialstate`      a vector of `state`s, one for each load case in the optimization problem, 
                      obtained from one or several previous `SweepX` analyses
- `maxiter=50`        maximum number of Newton-Raphson iterations 
- `maxΔa=1e-5`        "outer" convergence criteria: a norm on the scaled `A` increment 
- `maxΔy=1e-5`        "inner" convergence criteria: a norm on the scaled `Y=[ΛXU]` increment 
- `saveiter=false`    set to true so that the output `state` is a vector (over the Aiter) of 
                      vectors (over the steps) of `State`s of the model (for debugging 
                      non-convergence). 
- `maxLineIter=50`    maximum number of iterations in the linear search that ensure interior points   
- `β=0.5`             `β∈]0,1[`. In the line search, if conditions are not met, then a new line-iteration is done
                      with `s *= β` where  `β→0` is a hasty backtracking, while `β→1` stands its ground.            
- `γfac=0.5`          `γfac∈[0,1[`. At each iteration, the barrier parameter γ is taken as `γ = (∑ⁿᵢ₌₁ λᵢ gᵢ)/n*γfac` where
                      `(∑ⁿᵢ₌₁ λᵢ gᵢ)/n` is the complementary slackness, and `n` the number of inequality constraints.
- `γbot=1e-8`         `γ` will not be reduced to under the original complementary slackness divided by `γbot`,
                      to avoid conditioning problems.                                               

# Output

A vector of length equal to that of `initialstate` containing the state of the optimized model at each of these steps.                       

See also: [`solve`](@ref), [`SweepX`](@ref), [`setdof!`](@ref) 
"""
struct DirectXUA{NA,ND} <: AbstractSolver end 
function solve(::Type{DirectXUA{NA,ND}},pstate,verbose::𝕓,dbg;initialstate::Vector{<:State},
    maxiter::ℤ=50,maxΔy::ℝ=1e-5,maxΔa::ℝ=1e-5,
    saveiter::𝔹=false,
    maxLineIter::ℤ=50,β::𝕣=.5,γfac::𝕣=.5,γbot::𝕣=1e-8) where{NA,ND}

    model,dis             = initialstate[begin].model,initialstate[begin].dis
    out,asm,Ydofgr,Adofgr = prepare(AssemblyDirect    ,model,dis)
    out2,asm2             = prepare(AssemblyDirectLine,model,dis)

    cΔy²,cΔa²             = maxΔy^2,maxΔa^2
    nX,nU,nA              = getndof(model,(:X,:U,:A))
    nstep                 = length(initialstate)
    nV                    = nstep*(2*nX+nU) + nA
    nblock                = nstep + 1
    ΣLa                   = Vector{𝕣}(undef,nA   )

    # block                 = Matrix{SparseMatrixCSC{𝕣,𝕫}}(undef,nblock,nblock)
    # for step ∈ eachindex(initialstate)
    #     block[step  ,step  ]  = out.Lyy
    #     block[step  ,nblock]  = out.Lya
    #     block[nblock,step  ]  = out.Lay
    #     block[nblock,nblock]  = out.Laa
    # end
    i                     = 𝕫1(undef,4*length(initialstate))
    j                     = 𝕫1(undef,4*length(initialstate))
    v                     = Vector{typeof(out.Lyy)}(undef,4*length(initialstate))
    for step ∈ eachindex(initialstate)
        i[4step-3],j[4step-3],v[4step-3] = step  ,step  ,out.Lyy
        i[4step-2],j[4step-2],v[4step-2] = step  ,nblock,out.Lya
        i[4step-1],j[4step-1],v[4step-1] = nblock,step  ,out.Lay
        i[4step-0],j[4step-0],v[4step-0] = nblock,nblock,out.Laa
    end
    block = SparseBlocks(v,i,j)
    Lvv,blkasm            = prepare(block)
    Lv                    = 𝕣1(undef,nV)


    states                = [State{1,1,1}(i,(γ=0.,)) for i ∈ initialstate]
    if saveiter
        statess           = Vector{Vector{State{1,1,1,typeof((γ=0.,))}}}(undef,maxiter) 
        pstate[]          = statess
    else
        pstate[]          = states    
    end    

    Δy²                   = Vector{𝕣 }(undef,nstep)

    Σλg,npos              = 0.,0
    for (step,state)   ∈ enumerate(states) 
        assemble!(out2,asm2,dis,model,state,(dbg...,solver=:DirectXUA,phase=:preliminary,step=step))
        out2.ming ≤ 0 && muscadeerror(@sprintf("Initial point is not strictly primal-feasible at step %3d",step))
        out2.minλ ≤ 0 && muscadeerror(@sprintf("Initial point is not strictly dual-feasible at step %3d"  ,step))
        Σλg  += out2.Σλg
        npos += out2.npos
    end    
    γ = γ₀ = Σλg/max(1,npos)*γfac

    local LU
    for iter              = 1:maxiter
        verbose && @printf("    iteration %3d, γ=%g\n",iter,γ)

        zero!(Lvv)
        zero!(Lv )
        for (step,state)   ∈ enumerate(states)
            state.SP = (γ=γ ,)
            assemble!(out,asm,dis,model,state,(dbg...,solver=:DirectXUA,step=step,iter=iter))
            addin!(Lvv,out.Lyy,blkasm,step  ,step  )
            addin!(Lvv,out.Lya,blkasm,step  ,nblock)
            addin!(Lvv,out.Lay,blkasm,nblock,step  )
            addin!(Lvv,out.Laa,blkasm,nblock,nblock) # while A is step indep, Laa and La can be step dep
            addin!(Lv ,out.Ly ,blkasm,step         )
            addin!(Lv ,out.La ,blkasm,nblock       )
        end   

        try if iter==1 LU = lu(Lvv) 
        else           lu!(LU ,Lvv)
        end catch; muscadeerror(@sprintf("Lvv matrix factorization failed at iter=%i",iter));end
        Δv               = LU\Lv 

        Δa               = getblock(Δv,blkasm,nblock)
        Δa²              = sum(Δa.^2)
        for (step,state)   ∈ enumerate(states)
            Δy           = getblock(Δv,blkasm,step  )
            Δy²[step]    = sum(Δy.^2)
            decrement!(state,0,Δy,Ydofgr)
            decrement!(state,0,Δa,Adofgr)
        end    
        
        s  = 1.  
        local  Σλg,npos 
        for iline = 1:maxLineIter
            ΣLa              .= 0   
            minλ,ming         = ∞,∞
            Σλg,npos          = 0.,0
            for (step,state)  ∈ enumerate(states)
                assemble!(out2,asm2,dis,model,state,(dbg...,solver=:DirectXUAstepwise,phase=:linesearch,iter=iter,iline=iline,step=step))
                ΣLa         .+= out2.La 
                minλ          = min(minλ,out2.minλ)
                ming          = min(ming,out2.ming)
                Σλg          += out2.Σλg
                npos         += out2.npos
            end
            if minλ>0 && ming>0 
                verbose && @printf("    %3d line-iterations\n",iline)
                break#out of line search
            end
            iline==maxLineIter && muscadeerror(@sprintf("Line search failed at iter=%3d, iline=%3d, s=%7.1e",iter,iline,s))
            Δs                = s*(β-1)
            s                += Δs
            for (step,state)  ∈ enumerate(states)
                decrement!(state,0,Δs*getblock(Δv,blkasm,step),Ydofgr)
                decrement!(state,0,Δs*Δa                      ,Adofgr)
            end
        end
        γ                     = max(Σλg/max(1,npos)*γfac, γ₀*γbot)

        if saveiter
            statess[iter]     = copy.(states) 
        end

        if all(Δy².≤cΔy²)  && Δa²≤cΔa²  
            verbose && @printf("\n    DirectXUA converged in %3d iterations.\n",iter)
            verbose && @printf(  "    maxₜ(|ΔY|)=%7.1e  |ΔA|=%7.1e  \n",√(maximum(Δy²)),√(Δa²) )
            verbose && @printf(  "    nel=%d, nvariables=%d, nstep=%d, niter=%d\n",getnele(model),nV,nstep,iter)
            break#out of iter
        end
        iter<maxiter || muscadeerror(@sprintf("no convergence after %3d iterations. |ΔY|=%7.1e  |ΔA|=%7.1e \n",iter,√(maximum(Δy²)),√(Δa²)))
    end
    return
end


