const λxua = 1:4



mutable struct AssemblyDirectΛXU_A{Tλ,Tx,Tu,Ta,Tλa,Txa,Tua,Taa,Tλx,Txx,Tux,Tax,Tλu,Txu,Tuu,Tau,Tλa,Txa,Tua,Taa}  <:Assembly
    Lλ    :: Tλ
    Lx    :: Tx
    Lu    :: Tu
    La    :: Ta
    Lλa   :: Tλa
    Lxa   :: Txa
    Lua   :: Tua
    Laa   :: Taa
    Lλx   :: Tλx
    Lxx   :: Txx
    Lux   :: Tux
    Lax   :: Tax
    Lλu   :: Tλu
    Lxu   :: Txu
    Luu   :: Tuu
    Lau   :: Tau
    Lλa   :: Tλa
    Lxa   :: Txa
    Lua   :: Tua
    Laa   :: Taa
end   
function prepare(::Type{AssemblyDirectΛXU_A},model,dis,????) 
    Λdofgr             = allΛdofs(model,dis)
    Xdofgr             = allXdofs(model,dis)
    Udofgr             = allUdofs(model,dis)
    Adofgr             = allAdofs(model,dis)
    nX,nU,nA           = getndof(Xdofgr),getndof(Udofgr),getndof(Adofgr)
    narray,neletyp     = 20,getneletyp(model)
    asm                = Matrix{𝕫2}(undef,narray,neletyp)  
    Lλ                 = [asmvec!(view(asm, 1,:),Λdofgr,dis)                         for ider=1:nΛder             ] 
    Lx                 = [asmvec!(view(asm, 2,:),Xdofgr,dis)                         for ider=1:nXder             ] 
    Lu                 = [asmvec!(view(asm, 3,:),Udofgr,dis)                         for ider=1:nUder             ] 
    La                 = [asmvec!(view(asm, 4,:),Adofgr,dis)                         for ider=1:nAder             ] 
    Lλλ                = [asmmat!(view(asm, 5,:),view(asm,1,:),view(asm,1,:),nX,nX)  for ider=1:nΛder,jder=1:nΛder]
    Lxλ                = [asmmat!(view(asm, 6,:),view(asm,2,:),view(asm,1,:),nX,nX)  for ider=1:nXder,jder=1:nΛder]
    Luλ                = [asmmat!(view(asm, 7,:),view(asm,3,:),view(asm,1,:),nU,nX)  for ider=1:nUder,jder=1:nΛder]
    Laλ                = [asmmat!(view(asm, 8,:),view(asm,4,:),view(asm,1,:),nA,nX)  for ider=1:nAder,jder=1:nΛder]
    Lλx                = [asmmat!(view(asm, 9,:),view(asm,1,:),view(asm,2,:),nX,nX)  for ider=1:nΛder,jder=1:nXder]
    Lxx                = [asmmat!(view(asm,10,:),view(asm,2,:),view(asm,2,:),nX,nX)  for ider=1:nXder,jder=1:nXder]
    Lux                = [asmmat!(view(asm,11,:),view(asm,3,:),view(asm,2,:),nU,nX)  for ider=1:nUder,jder=1:nXder]
    Lax                = [asmmat!(view(asm,12,:),view(asm,4,:),view(asm,2,:),nA,nX)  for ider=1:nAder,jder=1:nXder]
    Lλu                = [asmmat!(view(asm,13,:),view(asm,1,:),view(asm,3,:),nX,nU)  for ider=1:nΛder,jder=1:nUder]
    Lxu                = [asmmat!(view(asm,14,:),view(asm,2,:),view(asm,3,:),nX,nU)  for ider=1:nXder,jder=1:nUder]
    Luu                = [asmmat!(view(asm,15,:),view(asm,3,:),view(asm,3,:),nU,nU)  for ider=1:nUder,jder=1:nUder]
    Lau                = [asmmat!(view(asm,16,:),view(asm,4,:),view(asm,3,:),nA,nU)  for ider=1:nAder,jder=1:nUder]
    Lλa                = [asmmat!(view(asm,17,:),view(asm,1,:),view(asm,4,:),nX,nA)  for ider=1:nΛder,jder=1:nAder]
    Lxa                = [asmmat!(view(asm,18,:),view(asm,2,:),view(asm,4,:),nX,nA)  for ider=1:nXder,jder=1:nAder]
    Lua                = [asmmat!(view(asm,19,:),view(asm,3,:),view(asm,4,:),nU,nA)  for ider=1:nUder,jder=1:nAder]
    Laa                = [asmmat!(view(asm,20,:),view(asm,4,:),view(asm,4,:),nA,nA)  for ider=1:nAder,jder=1:nAder]

    out                = AssemblyDirectΛXU_A(Lλ,Lx,Lu,La,Lλa,Lxa,Lua,Laa,Lλx,Lxx,Lux,Lax,Lλu,Lxu,Luu,Lau,Lλa,Lxa,Lua,Laa)
    return out,asm,Ydofgr,Adofgr
end
function zero!(out::AssemblyDirectΛXU_A)
   zero!.(out.Lλ )                
   zero!.(out.Lx )                
   zero!.(out.Lu )                
   zero!.(out.La )                
   zero!.(out.Lλλ) 
   zero!.(out.Lxλ)               
   zero!.(out.Luλ)               
   zero!.(out.Laλ)               
   zero!.(out.Lλx) 
   zero!.(out.Lxx)               
   zero!.(out.Lux)               
   zero!.(out.Lax)               
   zero!.(out.Lλu) 
   zero!.(out.Lxu)               
   zero!.(out.Luu)               
   zero!.(out.Lau)               
   zero!.(out.Lλa) 
   zero!.(out.Lxa)               
   zero!.(out.Lua)               
   zero!.(out.Laa)
end
function addin!(out::AssemblyDirectΛXU_A,asm,iele,scale,eleobj::E,Λ::SVector{Nx},X::NTuple{nXder,<:SVector{Nx}},
                                         U::NTuple{nUder,<:SVector{Nu}},A::SVector{Na},t,SP,dbg) where{E,nXder,nUder,Nx,Nu,Na} 
    # We make a distinction between nΛder==nAder==1, nXder=length(X), nUder=length(U) on the one hand, and mΞder ≤ nΞder.  This allows
    # 1) to freeze A for XU algo (or any class)
    # 2) not to compute cost on U′ or U′′ if these costs are known to be zero (same with X)                                      
    mΛder,mXder,mUder,mAder = 1,Nder,Nder,1 

    Λ∂ = nΛder==0 ? Λ : SVector{Nx}(  ∂²ℝ{1,Nz}(Λ[   iλ],  iλ)   for iλ=1:Nx)
    n       = nΛder*Nx
    X∂      = ntuple(Nder) do i 
        X∂ᵢ =           SVector{Nx}(  ∂²ℝ{1,Nz}(X[i][ix],n+ix)   for ix=1:Nx) 
        n  += Nx
    end
    U∂      = ntuple(Nder) do i 
        U∂ᵢ =           SVector{Nu}(  ∂²ℝ{1,Nz}(U[i][iu],n+iu)   for iu=1:Nu) 
        n  += Nu
    end
    A∂ = nAder==0 ? A : SVector{Na}(  ∂²ℝ{1,Nz}(A[   ia],n+ia)   for ia=1:Na)

    L,FB    = getlagrangian(eleobj, Λ∂,X∂,U∂,A∂,t,SP,dbg)
  
    ∇L      = ∂{2,Nz}(L)
    for α∈λxua, i=1:Nder[α]
        ip = 
        add_value!(out.L1[α][i] ,asm.vec[α],iele,∇L,ip   )
        for β∈λxua, j=1:Nder[β]
            jp = 
            add_∂!{1}( out.L2[α,β][ider,jder],asm.mat[α,β],iele,∇L,ip,jp)
        end
    end

end

###--------------------- ASMDirectXUAstepwiseline: for line search

mutable struct AssemblyDirectΛXU_Aline{Ty,Ta} <:Assembly
    Ly    :: Ty
    La    :: Ta
    ming  :: 𝕣
    minλ  :: 𝕣
    Σλg   :: 𝕣
    npos  :: 𝕫
end   
function prepare(::Type{AssemblyDirectΛXU_Aline},model,dis,wantA,Nder) 
    Ydofgr             = allΛXUdofs(model,dis)
    Adofgr             = wantA ? allAdofs(model,dis) : nodofs(model,dis)
    narray,neletyp     = 2,getneletyp(model)
    asm                = Matrix{𝕫2}(undef,narray,neletyp)  
    Ly                 = [asmvec!(view(asm,1,:),Ydofgr,dis)  for ider=0:Nder] 
    La                 =  asmvec!(view(asm,2,:),Adofgr,dis) 
    out                = AssemblyDirectΛXU_Aline(Ly,La,∞,∞,0.,0) # sequential
    return out,asm,Ydofgr,Adofgr
end
function zero!(out::AssemblyDirectΛXU_Aline)
    zero!.(out.Ly)
    zero!( out.La)
    out.ming = ∞    
    out.minλ = ∞
    out.Σλg  = 0.
    out.npos = 0    
end
function addin!(out::AssemblyDirectΛXU_Aline,asm,iele,scale,eleobj::E,Λ,X::NTuple{nXder,<:SVector{Nx}},
                                              U::NTuple{nUder,<:SVector{Nu}},A::SVector{Na},t,SP,dbg) where{E,nXder,Nx,nUder,Nu,Na}
    Ny              = 2Nx+Nu                           # Y=[Λ;X;U]   
    Nz              = 2Nx+Nu+Na                        # Z = [Y;A]=[Λ;X;U;A]       
    scaleZ          = SVector(scale.Λ...,scale.X...,scale.U...,scale.A...)
    ΔZ              = δ{1,Nz,𝕣}(scaleZ)                 
    iλ,ix,iu,ia     = gradientpartition(Nx,Nx,Nu,Na) # index into element vectors ΔZ and Lz
    ΔΛ,ΔX,ΔU,ΔA     = view(ΔZ,iλ),view(ΔZ,ix),view(ΔZ,iu),view(ΔZ,ia) 
    L,FB            = getlagrangian(eleobj, ∂0(Λ)+ΔΛ, (∂0(X)+ΔX,),(∂0(U)+ΔU,),A+ΔA,t,SP,dbg)
    ∇L              = ∂{1,Nz}(L)
    add_value!(out.Ly ,asm[1],iele,∇L,1:Ny)
    add_value!(out.La ,asm[2],iele,∇L,ia  )
    if hasfield(typeof(FB),:mode) && FB.mode==:positive
        out.ming   = min(out.ming,VALUE(FB.g))
        out.minλ   = min(out.minλ,VALUE(FB.λ))
        out.Σλg   += VALUE(FB.g)*VALUE(FB.λ)
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
    out,asm,Ydofgr,Adofgr = prepare(AssemblyDirectΛXU_A    ,model,dis)
    out2,asm2,_     ,_    = prepare(AssemblyDirectΛXU_Aline,model,dis)

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


