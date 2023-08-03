
mutable struct AssemblyStaticΛXU_A{Ty,Ta,Tyy,Tya,Taa}  <:Assembly
    Ly    :: Ty
    La    :: Ta
    Lyy   :: Tyy 
    Lya   :: Tya 
    Laa   :: Taa
end   
function prepare(::Type{AssemblyStaticΛXU_A},model,dis) 
    Ydofgr             = allΛXUdofs(model,dis)
    Adofgr             = allAdofs(  model,dis)
    nY,nA              = getndof(Ydofgr),getndof(Adofgr)
    narray,neletyp     = 5,getneletyp(model)
    asm                = Matrix{𝕫2}(undef,narray,neletyp)  
    Ly                 = asmvec!(view(asm,1,:),Ydofgr,dis) 
    La                 = asmvec!(view(asm,2,:),Adofgr,dis) 
    Lyy                = asmmat!(view(asm,3,:),view(asm,1,:),view(asm,1,:),nY,nY) 
    Lya                = asmfullmat!(view(asm,4,:),view(asm,1,:),view(asm,2,:),nY,nA) 
    Laa                = asmfullmat!(view(asm,5,:),view(asm,2,:),view(asm,2,:),nA,nA)  
    out                = AssemblyStaticΛXU_A(Ly,La,Lyy,Lya,Laa)
    return out,asm,Ydofgr,Adofgr
end
function zero!(out::AssemblyStaticΛXU_A)
    zero!(out.Ly )
    zero!(out.La )
    zero!(out.Lyy)
    zero!(out.Lya)
    zero!(out.Laa)
end
function add!(out1::AssemblyStaticΛXU_A,out2::AssemblyStaticΛXU_A) 
    add!(out1.Ly ,out2.Ly )
    add!(out1.La ,out2.La )
    add!(out1.Lyy,out2.Lyy)
    add!(out1.Lya,out2.Lya)
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
    add_∂!{1}( out.Laa,asm[5],iele,∇L,ia,ia)
end

###--------------------- ASMstaticXUAline: for line search

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
    maxAiter::ℤ=50,maxΔy::ℝ=1e-5,maxLy::ℝ=∞,maxΔa::ℝ=1e-5,maxLa::ℝ=∞,
    saveiter::𝔹=false,
    maxLineIter::ℤ=50,α::𝕣=.1,β::𝕣=.5,γfac::𝕣=.5)

    model,dis          = initialstate[begin].model,initialstate[begin].dis
    out1,asm1,Ydofgr,Adofgr = prepare(AssemblyStaticΛXU_A    ,model,dis)
    out2,asm2,_     ,_      = prepare(AssemblyStaticΛXU_Aline,model,dis)
    if saveiter
        statess        = allocate(pstate,Vector{Vector{State{1,1,1,typeof((γ=0.,))}}}(undef,maxAiter)) 
        states = [State{1,1,1}(i,(γ=0.,)) for i ∈ initialstate]
    else
        states         = allocate(pstate,[State{1,1,1}(i,(γ=0.,)) for i ∈ initialstate]) 
    end    
    cΔy²,cLy²,cΔa²,cLa²= maxΔy^2,maxLy^2,maxΔa^2,maxLa^2
    nA,nStep           = getndof(model,:A),length(initialstate)
    La                 = Vector{𝕣 }(undef,nA   )
    Laa                = Matrix{𝕣 }(undef,nA,nA)
    Δy                 = Vector{𝕣1}(undef,nStep)
    y∂a                = Vector{𝕣2}(undef,nStep)
    Δy²,Ly²,La²        = copies(3,Vector{𝕣}(undef,nStep))
    cAiter,cYiter      = 0,0
    local facLyy, Δa

    Σλg,npos           = 0.,0
    for (step,state)   ∈ enumerate(states) 
        assemble!(out2,asm2,dis,model,state,(dbg...,solver=:StaticXUA,phase=:preliminary,step=step))
        out2.ming ≤ 0 && muscadeerror(@sprintf("Initial point is not strictly primal-feasible at step %3d",step))
        out2.minλ ≤ 0 && muscadeerror(@sprintf("Initial point is not strictly dual-feasible at step %3d"  ,step))
        Σλg  += out2.Σλg
        npos += out2.npos
    end    
    γ = Σλg/npos * γfac
    for state ∈ states
        state.SP = (γ=γ,)
    end

    for iAiter          = 1:maxAiter
        if saveiter
            statess[iAiter] = [State{1,1,1}(i,(γ=0.,)) for i ∈ (iAiter==1 ? initialstate : statess[iAiter-1])]
            states          = statess[iAiter]
        end
        verbose && @printf "    A-iteration %3d\n" iAiter

        La            .= 0
        Laa           .= 0
        for (step,state)   ∈ enumerate(states)
            assemble!(out1,asm1,dis,model,state,(dbg...,solver=:StaticXUA,phase=:direction,iAiter=iAiter,step=step))
            try if iAiter==1 && step==1
                facLyy = lu(out1.Lyy) 
            else
                lu!(facLyy,out1.Lyy)
            end catch; muscadeerror(@sprintf("Lyy matrix factorization failed at step=%i, iAiter=%i",step,iAiter));end
            Δy[ step]  = -(facLyy\out1.Ly)  
            y∂a[step]  = -(facLyy\out1.Lya) 
            La       .+= out1.La  + out1.Lya' * Δy[ step]  
            Laa      .+= out1.Laa + out1.Lya' * y∂a[step]
            Ly²[step]  = sum(out1.Ly.^2) 
            La²[step]  = sum(out1.La.^2) 
        end   
        Lz² = sum(La²)+sum(Ly²)

        try 
            Δa         = -(Laa\La) 
        catch; muscadeerror(@sprintf("Laa\\La solution failed at iAiter=%i",iAiter));end
        Δa²            = sum(Δa.^2)

        for (step,state)   ∈ enumerate(states)
            Δy[step] .+= y∂a[step] * Δa
            Δy²[step]  = sum(Δy[step].^2)
            increment!(state,0,Δy[step],Ydofgr)
            increment!(state,0,Δa,Adofgr)
        end    
        
        s  = 1.    
        for iline = 1:maxLineIter
            Lz²line,minλ,ming = 0.,∞,∞
            for (step,state)   ∈ enumerate(states)
                assemble!(out2,asm2,dis,model,state,(dbg...,solver=:StaticXUA,phase=:linesearch,iAiter=iAiter,iline=iline,step=step))
                Lz²line += sum(out2.Ly.^2) + sum(out2.La.^2)
                minλ = min(minλ,out2.minλ)
                ming = min(ming,out2.ming)
            end
            @show iAiter,iline,s,Lz²-Lz²line
#            @show minλ > 0 , ming > 0 , Lz²line ≤ Lz²*(1-α*s)^2
#            @show Lz²line , sum(La²),sum(Ly²)
#            minλ > 0 && ming > 0 && Lz²line ≤ Lz²*(1-α*s)^2 && break#out of line search
            minλ > 0 && ming > 0  && break#out of line search
            iline==maxLineIter && muscadeerror(@sprintf("Line search failed at iAiter=%3d, iline=%3d, s=%7.1e",iAiter,iline,s))
            Δs = s*(β-1)
            s += Δs
            for (step,state)   ∈ enumerate(states)
                increment!(state,0,Δs*Δy[step],Ydofgr)
                increment!(state,0,Δs*Δa      ,Adofgr)
            end
        end

        if all(Δy².*s^2 .≤cΔy²) && all(Ly².≤cLy²) && Δa²*s^2 .≤cΔa² && all(La².≤cLa²) 
            cAiter    = iAiter
            verbose && @printf "\n    StaticXUA converged in %3d A-iterations.\n" iAiter
            verbose && @printf "    maxₜ(|ΔY|)=%7.1e  maxₜ(|∇L/∂Y|)=%7.1e  |ΔA|=%7.1e  |∇L/∂A|=%7.1e\n" √(maximum(Δy²)) √(maximum(Ly²)) √(Δa²) √(La²)
            break#out of iAiter
        end
        iAiter==maxAiter && muscadeerror(@sprintf("no convergence at iAiter=%3d, |ΔY|=%7.1e |Ly|=%7.1e |ΔA|=%7.1e |La|=%7.1e\n",iAiter,√(maximum(Δy²)),√(maximum(Ly²)),√(Δa²),√(La²)))

        for state ∈ states
            state.SP     = (γ=state.SP.γ * γfac,)
        end
    end
    verbose && @printf "\n    nel=%d, ndof=%d, nstep=%d, nAiter=%d\n" getnele(model) getndof(Adofgr) nStep cAiter
    verbose && @printf "\n    nYiter=%d, nYiter/(nstep*nAiter)=%5.2f\n" cYiter cYiter/nStep/cAiter
    return
end


