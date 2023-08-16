
mutable struct AssemblyStaticΛXU_A{Ty,Ta,Tyy,Tya,Tay,Taa}  <:Assembly
    Ly    :: Ty
    La    :: Ta
    Lyy   :: Tyy 
    Lya   :: Tya 
    Lay   :: Tay 
    Laa   :: Taa
end   
function prepare(::Type{AssemblyStaticΛXU_A},model,dis) 
    Ydofgr             = allΛXUdofs(model,dis)
    Adofgr             = allAdofs(  model,dis)
    nY,nA              = getndof(Ydofgr),getndof(Adofgr)
    narray,neletyp     = 6,getneletyp(model)
    asm                = Matrix{𝕫2}(undef,narray,neletyp)  
    Ly                 = asmvec!(view(asm,1,:),Ydofgr,dis) 
    La                 = asmvec!(view(asm,2,:),Adofgr,dis) 
    Lyy                = asmmat!(view(asm,3,:),view(asm,1,:),view(asm,1,:),nY,nY) 
    Lya                = asmmat!(view(asm,4,:),view(asm,1,:),view(asm,2,:),nY,nA) 
    Lay                = asmmat!(view(asm,5,:),view(asm,2,:),view(asm,1,:),nA,nY) 
    Laa                = asmmat!(view(asm,6,:),view(asm,2,:),view(asm,2,:),nA,nA)  
    out                = AssemblyStaticΛXU_A(Ly,La,Lyy,Lya,Lay,Laa)
    return out,asm,Ydofgr,Adofgr
end
function zero!(out::AssemblyStaticΛXU_A)
    zero!(out.Ly )
    zero!(out.La )
    zero!(out.Lyy)
    zero!(out.Lya)
    zero!(out.Lay)
    zero!(out.Laa)
end
function add!(out1::AssemblyStaticΛXU_A,out2::AssemblyStaticΛXU_A) 
    add!(out1.Ly,out2.Ly)
    add!(out1.La,out2.La)
    add!(out1.Lyy,out2.Lyy)
    add!(out1.Lya,out2.Lya)
    add!(out1.Lay,out2.Lay)
    add!(out1.Laa,out2.Laa)
end
function addin!(out::AssemblyStaticΛXU_A,asm,iele,scale,eleobj::E,Λ,X::NTuple{Nxder,<:SVector{Nx}},
                                         U::NTuple{Nuder,<:SVector{Nu}},A::SVector{Na},t,SP,dbg) where{E,Nxder,Nx,Nuder,Nu,Na} # TODO make Nx,Nu,Na types
    Ny              = 2Nx+Nu                           # Y=[Λ;X;U]   
    Nz              = 2Nx+Nu+Na                        # Z = [Y;A]=[Λ;X;U;A]       
    scaleZ          = SVector(scale.Λ...,scale.X...,scale.U...,scale.A...)
    ΔZ              = variate{2,Nz}(δ{1,Nz,𝕣}(scaleZ),scaleZ)                 
    iλ,ix,iu,ia     = gradientpartition(Nx,Nx,Nu,Na) # index into element vectors ΔZ and Lz
    iy              = 1:Ny  
    ΔΛ,ΔX,ΔU,ΔA     = view(ΔZ,iλ),view(ΔZ,ix),view(ΔZ,iu),view(ΔZ,ia) # TODO Static?
    L,χn,FB         = getlagrangian(implemented(eleobj)...,eleobj, ∂0(Λ)+ΔΛ, (∂0(X)+ΔX,),(∂0(U)+ΔU,),A+ΔA,t,nothing,nothing,SP,dbg)
    ∇L              = ∂{2,Nz}(L)
    add_value!(out.Ly ,asm[1],iele,∇L,iy   )
    add_value!(out.La ,asm[2],iele,∇L,ia   )
    add_∂!{1}( out.Lyy,asm[3],iele,∇L,iy,iy)
    add_∂!{1}( out.Lya,asm[4],iele,∇L,iy,ia)
    add_∂!{1}( out.Lay,asm[5],iele,∇L,ia,iy)
    add_∂!{1}( out.Laa,asm[6],iele,∇L,ia,ia)
end

###--------------------- ASMStaticXUAstepwiseline: for line search

mutable struct AssemblyStaticΛXU_Aline{Ty,Ta} <:Assembly
    Ly    :: Ty
    La    :: Ta
    ming  :: 𝕣
    minλ  :: 𝕣
    Σλg   :: 𝕣
    npos  :: 𝕫
end   
function prepare(::Type{AssemblyStaticΛXU_Aline},model,dis) 
    Ydofgr             = allΛXUdofs(model,dis)
    Adofgr             = allAdofs(  model,dis)
    narray,neletyp     = 2,getneletyp(model)
    asm                = Matrix{𝕫2}(undef,narray,neletyp)  
    Ly                 = asmvec!(view(asm,1,:),Ydofgr,dis) 
    La                 = asmvec!(view(asm,2,:),Adofgr,dis) 
    out                = AssemblyStaticΛXU_Aline(Ly,La,∞,∞,0.,0) # sequential
    return out,asm,Ydofgr,Adofgr
end
function zero!(out::AssemblyStaticΛXU_Aline)
    zero!(out.Ly)
    zero!(out.La)
    out.ming = ∞    
    out.minλ = ∞
    out.Σλg  = 0.
    out.npos = 0    
end
function add!(out1::AssemblyStaticΛXU_Aline,out2::AssemblyStaticΛXU_Aline) 
    add!(out1.Ly,out2.Ly)
    add!(out1.La,out2.La)
    out1.ming = min(out1.ming,out2.ming)
    out1.minλ = min(out1.minλ,out2.minλ)
    out1.Σλg += out2.Σλg
    out1.npos+= out2.npos
end
function addin!(out::AssemblyStaticΛXU_Aline,asm,iele,scale,eleobj::E,Λ,X::NTuple{Nxder,<:SVector{Nx}},
                                              U::NTuple{Nuder,<:SVector{Nu}},A::SVector{Na},t,SP,dbg) where{E,Nxder,Nx,Nuder,Nu,Na}
    Ny              = 2Nx+Nu                           # Y=[Λ;X;U]   
    Nz              = 2Nx+Nu+Na                        # Z = [Y;A]=[Λ;X;U;A]       
    scaleZ          = SVector(scale.Λ...,scale.X...,scale.U...,scale.A...)
    ΔZ              = δ{1,Nz,𝕣}(scaleZ)                 
    iλ,ix,iu,ia     = gradientpartition(Nx,Nx,Nu,Na) # index into element vectors ΔZ and Lz
    ΔΛ,ΔX,ΔU,ΔA     = view(ΔZ,iλ),view(ΔZ,ix),view(ΔZ,iu),view(ΔZ,ia) # TODO Static?
    L,χn,FB         = getlagrangian(implemented(eleobj)...,eleobj, ∂0(Λ)+ΔΛ, (∂0(X)+ΔX,),(∂0(U)+ΔU,),A+ΔA,t,nothing,nothing,SP,dbg)
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
	StaticXUA

A non-linear static solver for optimisation FEM.
The current algorithm does not handle element memory. 

An analysis is carried out by a call with the following syntax:

```
initialstate    = initialize!(model)
stateX          = solve(StaticX  ;initialstate=initialstate,time=[0.,1.])
stateXUA        = solve(StaticXUA;initialstate=stateX)
```

# Named arguments
- `dbg=(;)`           a named tuple to trace the call tree (for debugging)
- `verbose=true`      set to false to suppress printed output (for testing)
- `silenterror=false` set to true to suppress print out of error (for testing) 
- `initialstate`      a vector of `state`s, one for each load case in the optimization problem, 
                      obtained from one or several previous `StaticX` analyses
- `maxAiter=50`       maximum number of "outer" Newton-Raphson iterations over `A` 
- `maxΔa=1e-5`        "outer" convergence criteria: a norm on the scaled `A` increment 
- `maxLa=∞`           "outer" convergence criteria: a norm on the scaled `La` residual
- `maxYiter=0`        maximum number of "inner" Newton-Raphson iterations over `X` 
                      and `U` for every value of `A`.  Experience so far is that these inner
                      iterations do not increase performance, so the default is "no inner 
                      iterations".   
- `maxΔy=1e-5`        "inner" convergence criteria: a norm on the scaled `Y=[XU]` increment 
- `maxLy=∞`           "inner" convergence criteria: a norm on the scaled `Ly=[Lx,Lu]` residual
- `saveiter=false`    set to true so that the output `state` is a vector (over the Aiter) of 
                      vectors (over the steps) of `State`s of the model (for debugging 
                      non-convergence). 
- `γ0=1.`             an initial value of the barrier coefficient for the handling of contact
                      using an interior point method
- `γfac1=0.5`         at each iteration, the barrier parameter γ is multiplied 
- `γfac2=100.`        by γfac1*exp(-min(αᵢ)/γfac2)^2), where αᵢ is computed by the i-th
                      interior point savvy element as αᵢ=abs(λ-g)/γ                                               

# Output

A vector of length equal to that of `initialstate` containing the state of the optimized model at each of these steps.                       

See also: [`solve`](@ref), [`StaticX`](@ref) 
"""
struct StaticXUA <: AbstractSolver end 
function solve(::Type{StaticXUA},pstate,verbose::𝕓,dbg;initialstate::Vector{<:State},
    maxiter::ℤ=50,maxΔy::ℝ=1e-5,maxΔa::ℝ=1e-5,
    saveiter::𝔹=false,γ0::𝕣=1.,γfac1::𝕣=.5,γfac2::𝕣=100.)

    model,dis             = initialstate[begin].model,initialstate[begin].dis
    out,asm,Ydofgr,Adofgr = prepare(AssemblyStaticΛXU_A,model,dis)

    cΔy²,cΔa²             = maxΔy^2,maxΔa^2
    nX,nU,nA              = getndof(model,(:X,:U,:A))
    nstep                 = length(initialstate)
    nV                    = nstep*(2*nX+nU) + nA
    nblock                = nstep + 1

    block                 = Matrix{SparseMatrixCSC{𝕣,𝕫}}(undef,nblock,nblock)
    for step ∈ eachindex(initialstate)
        block[step  ,step  ]  = out.Lyy
        block[step  ,nblock]  = out.Lya
        block[nblock,step  ]  = out.Lay
        block[nblock,nblock]  = out.Laa
    end
    Lvv,blkasm            = blocksparse(block)
    Lv                    = 𝕣1(undef,nV)

    if saveiter; states   = allocate(pstate,Vector{Vector{State{1,1,1,typeof((γ=0.,))}}}(undef,maxiter)) 
    else         state    = allocate(pstate,[State{1,1,1}(i,(γ=γ0,)) for i ∈ initialstate]) # deepcopy dofs from initstate (including A) 
    end    
    Δy²                   = Vector{𝕣 }(undef,nstep)

    local LU
    for iter              = 1:maxiter
        verbose && @printf "    iteration %3d\n" iter
        if saveiter
            states[iter]  = [State{1,1,1}(i,(γ=0.,)) for i ∈ (iter==1 ? initialstate : states[iter-1])]
            state         = states[iter]
        end
        zero!(Lvv)
        zero!(Lv )
        for (step,s)   ∈ enumerate(state)
            assemble!(out,asm,dis,model,s,(dbg...,solver=:StaticXUA,step=step,iter=iter))
            addin!(Lvv,out.Lyy,blkasm,step  ,step  )
            addin!(Lvv,out.Lya,blkasm,step  ,nblock)
            addin!(Lvv,out.Lay,blkasm,nblock,step  )
            addin!(Lvv,out.Laa,blkasm,nblock,nblock) # while A is step indep, Laa and La can be step dep
            addin!(Lv ,out.Ly ,blkasm,step         )
            addin!(Lv ,out.La ,blkasm,nblock       )
        end   

        try if iter==1 LU = lu(Lvv) 
        else           lu!(LU ,Lvv)
        end catch; muscadeerror(@sprintf("Lvv matrix factorization failed at iAiter=%i",iAiter));end

        Δv               = LU\Lv 

        Δa               = getblock(Δv,blkasm,nblock)
        Δa²              = sum(Δa.^2)
        for (step,s)   ∈ enumerate(state)
            Δy           = getblock(Δv,blkasm,step  )
            Δy²[step]    = sum(Δy.^2)
            decrement!(s,0,Δy,Ydofgr)
            decrement!(s,0,Δa,Adofgr)
            s.SP = (γ= s.SP.γ* γfac1*exp(-(out.α/γfac2)^2),)
        end    
        
        if all(Δy².≤cΔy²)  && Δa²≤cΔa²  
            verbose && @printf "\n    StaticXUA converged in %3d iterations.\n" iter
            verbose && @printf "    maxₜ(|ΔY|)=%7.1e  |ΔA|=%7.1e  \n" √(maximum(Δy²)) √(Δa²) 
            break#out of iter
        end
        iter<maxiter || muscadeerror(@sprintf("no convergence after %3d iterations. |ΔY|=%7.1e  |ΔA|=%7.1e \n",iAiter,√(maximum(Δy²)),√(Δa²)))
    end
    verbose && @printf "\n    nel=%d, ndof=%d, nstep=%d, niter=%d\n" getnele(model) getndof(Adofgr) nstep cAiter
    verbose && @printf "\n    nYiter=%d, nYiter/(nstep*nAiter)=%5.2f\n" cYiter cYiter/nstep/cAiter
    return
end


