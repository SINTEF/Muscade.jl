
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
    L,χn,FB         = getlagrangian(eleobj, ∂0(Λ)+ΔΛ, (∂0(X)+ΔX,),(∂0(U)+ΔU,),A+ΔA,t,nothing,SP,dbg)
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
    L,χn,FB         = getlagrangian(eleobj, ∂0(Λ)+ΔΛ, (∂0(X)+ΔX,),(∂0(U)+ΔU,),A+ΔA,t,nothing,SP,dbg)
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
setdof!(initialstate,1.;class=:U,field=:λcsr)
stateX          = solve(StaticX  ;initialstate=initialstate,time=[0.,1.])
stateXUA        = solve(StaticXUA;initialstate=stateX)
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
                      obtained from one or several previous `StaticX` analyses
- `maxiter=50`        maximum number of Newton-Raphson iterations 
- `maxΔa=1e-5`        "outer" convergence criteria: a norm on the scaled `A` increment 
- `maxΔy=1e-5`        "inner" convergence criteria: a norm on the scaled `Y=[ΛXU]` increment 
- `saveiter=false`    set to true so that the output `state` is a vector (over the Aiter) of 
                      vectors (over the steps) of `State`s of the model (for debugging 
                      non-convergence). 
- `maxLineIter=50`    maximum number of iterations in the linear search that ensure interior points   
- `α=0.1`             α∈[0,1]. Besides primal and dual feasibility, the line search ensures that
                      `|Lvline|≤|Lv|*(1-α*s)` where `|Lv|` and `|Lvline|` are the gradients of the 
                      Lagrangian, repsectively at the Newton step and the line search. `s∈]0,1]`
                      is the step reduction factor in the line search. `α→0` is lenient, `α→1` stringent.
- `β=0.5`             `β∈]0,1[`. In the line search, if conditions are not met, then a new line-iteration is done
                      with `s *= β` where  `β→0` is a hasty backtracking, while `β→1` stands its ground.            
- `γfac=0.5`          `γfac∈[0,1[`. At each iteration, the barrier parameter γ is taken as `γ = (∑ⁿᵢ₌₁ λᵢ gᵢ)/n*γfac` where
                      `(∑ⁿᵢ₌₁ λᵢ gᵢ)/n` is the complementary slackness, and `n` the number of inequality constraints.
- `γbot=1e-8`         `γ` will not be reduced to under the original complementary slackness divided by `γbot`,
                      to avoid conditioning problems.                                               

# Output

A vector of length equal to that of `initialstate` containing the state of the optimized model at each of these steps.                       

See also: [`solve`](@ref), [`StaticX`](@ref), [`setdof!`](@ref) 
"""
struct StaticXUA <: AbstractSolver end 
function solve(::Type{StaticXUA},pstate,verbose::𝕓,dbg;initialstate::Vector{<:State},
    maxiter::ℤ=50,maxΔy::ℝ=1e-5,maxΔa::ℝ=1e-5,
    saveiter::𝔹=false,
    maxLineIter::ℤ=50,α::𝕣=.1,β::𝕣=.5,γfac::𝕣=.5,γbot::𝕣=1e-8)

    model,dis             = initialstate[begin].model,initialstate[begin].dis
    out,asm,Ydofgr,Adofgr = prepare(AssemblyStaticΛXU_A    ,model,dis)
    out2,asm2,_     ,_    = prepare(AssemblyStaticΛXU_Aline,model,dis)

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
        assemble!(out2,asm2,dis,model,state,(dbg...,solver=:StaticXUA,phase=:preliminary,step=step))
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
            assemble!(out,asm,dis,model,state,(dbg...,solver=:StaticXUA,step=step,iter=iter))
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
            Lv²line,minλ,ming = 0.,∞,∞
            Σλg,npos          = 0.,0
            for (step,state)  ∈ enumerate(states)
                assemble!(out2,asm2,dis,model,state,(dbg...,solver=:StaticXUAstepwise,phase=:linesearch,iter=iter,iline=iline,step=step))
                ΣLa         .+= out2.La    
                Lv²line      += sum(out2.Ly.^2) 
                minλ          = min(minλ,out2.minλ)
                ming          = min(ming,out2.ming)
                Σλg          += out2.Σλg
                npos         += out2.npos
            end
            Lv²line          += sum(ΣLa.^2)
            if minλ>0 && ming>0 && Lv²line≤sum(Lv.^2)*(1-α*s)^2   
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
            statess[iter]     = deepcopy(states) # TODO this deepcopies the model!
        end

        if all(Δy².≤cΔy²)  && Δa²≤cΔa²  
            verbose && @printf("\n    StaticXUA converged in %3d iterations.\n",iter)
            verbose && @printf(  "    maxₜ(|ΔY|)=%7.1e  |ΔA|=%7.1e  \n",√(maximum(Δy²)),√(Δa²) )
            verbose && @printf(  "    nel=%d, nvariables=%d, nstep=%d, niter=%d\n",getnele(model),nV,nstep,iter)
            break#out of iter
        end
        iter<maxiter || muscadeerror(@sprintf("no convergence after %3d iterations. |ΔY|=%7.1e  |ΔA|=%7.1e \n",iter,√(maximum(Δy²)),√(Δa²)))
    end
    return
end


