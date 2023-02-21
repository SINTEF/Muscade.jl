
mutable struct OUTstaticΛXU_A{Ty,Ta,Tyy,Tya,Taa}  
    Ly    :: Ty
    La    :: Ta
    Lyy   :: Tyy 
    Lya   :: Tya 
    Laa   :: Taa
    α     :: 𝕣
end   
function prepare(::Type{OUTstaticΛXU_A},model,dis) 
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
    out                = OUTstaticΛXU_A(Ly,La,Lyy,Lya,Laa,0.)
    return out,asm,Adofgr,Ydofgr
end
function zero!(out::OUTstaticΛXU_A)
    zero!(out.Ly )
    zero!(out.La )
    zero!(out.Lyy)
    zero!(out.Lya)
    zero!(out.Laa)
    out.α = ∞    
end
function addin!(out::OUTstaticΛXU_A,asm,iele,scale,eleobj::E,Λ,X::NTuple{Nxdir,<:SVector{Nx}},
                                                               U::NTuple{Nudir,<:SVector{Nu}},A::SVector{Na}, t,γ,dbg) where{E,Nxdir,Nx,Nudir,Nu,Na} # TODO make Nx,Nu,Na types
    Ny              = 2Nx+Nu                           # Y=[Λ;X;U]   
    Nz              = 2Nx+Nu+Na                        # Z = [Y;A]=[Λ;X;U;A]       
    ΔZ              = variate{2,Nz}(δ{1,Nz,𝕣}())                 
    iλ,ix,iu,ia     = gradientpartition(Nx,Nx,Nu,Na) # index into element vectors ΔZ and Lz
    iy              = 1:Ny  
    ΔΛ,ΔX,ΔU,ΔA     = view(ΔZ,iλ),view(ΔZ,ix),view(ΔZ,iu),view(ΔZ,ia) # TODO Static?
    L,α             = scaledlagrangian(scale,eleobj, Λ+ΔΛ, (∂0(X)+ΔX,),(∂0(U)+ΔU,),A+ΔA, t,γ,dbg)
    ∇L              = ∂{2,Nz}(L)
    add_value!(out.Ly ,asm[1],iele,∇L,iy   )
    add_value!(out.La ,asm[2],iele,∇L,ia   )
    add_∂!{1}( out.Lyy,asm[3],iele,∇L,iy,iy)
    add_∂!{1}( out.Lya,asm[4],iele,∇L,iy,ia)
    add_∂!{1}( out.Laa,asm[5],iele,∇L,ia,ia)
    out.α           = min(out.α,α)
end

#------------------------------------

mutable struct OUTstaticΛXU{Ty,Tyy}  
    Ly    :: Ty
    Lyy   :: Tyy 
    α     :: 𝕣
end   
function prepare(::Type{OUTstaticΛXU},model,dis) 
    Ydofgr             = allΛXUdofs(model,dis)
    nY                 = getndof(Ydofgr)
    narray,neletyp     = 2,getneletyp(model)
    asm                = Matrix{𝕫2}(undef,narray,neletyp)  
    Ly                 = asmvec!(view(asm,1,:),Ydofgr,dis) 
    Lyy                = asmmat!(view(asm,2,:),view(asm,1,:),view(asm,1,:),nY,nY) 
    out                = OUTstaticΛXU(Ly,Lyy,0.)
    return out,asm,Ydofgr
end
function zero!(out::OUTstaticΛXU)
    zero!(out.Ly )
    zero!(out.Lyy)
    out.α = ∞    
end
function addin!(out::OUTstaticΛXU,asm,iele,scale,eleobj::E,Λ,X::NTuple{Nxdir,<:SVector{Nx}},
                                                             U::NTuple{Nudir,<:SVector{Nu}},A, t,γ,dbg) where{E,Nxdir,Nx,Nudir,Nu}
    Ny              = 2Nx+Nu                           # Y=[Λ;X;U]  TODO compile time? 
    if Ny==0; return end # don't waste time on Acost elements...    
    ΔY              = variate{2,Ny}(δ{1,Ny,𝕣}())                 
    iλ,ix,iu,_      = gradientpartition(Nx,Nx,Nu,0) # index into element vectors ΔY and Ly
    ΔΛ,ΔX,ΔU        = view(ΔY,iλ),view(ΔY,ix),view(ΔY,iu)
    L,α             = scaledlagrangian(scale,eleobj, Λ+ΔΛ, (∂0(X)+ΔX,),(∂0(U)+ΔU,),A, t,γ,dbg)
    ∇L              = ∂{2,Ny}(L)
    add_value!(out.Ly ,asm[1],iele,∇L)
    add_∂!{1}( out.Lyy,asm[2],iele,∇L)
    out.α          = min(out.α,α)
end

#------------------------------------

function staticXUA(pstate,dbg;initialstate::Vector{State},
    maxAiter::ℤ=50,maxYiter::ℤ=0,maxΔy::ℝ=1e-5,maxLy::ℝ=∞,maxΔa::ℝ=1e-5,maxLa::ℝ=∞,γ0::𝕣=1.,γfac1::𝕣=.5,γfac2::𝕣=100.,verbose::𝕓=true)

    model,dis          = initialstate[begin].model,initialstate[begin].dis
    out1,asm1,Ydofgr   = prepare(OUTstaticΛXU  ,model,dis)
    out2,asm2,Adofgr,_ = prepare(OUTstaticΛXU_A,model,dis)
    state              = allocate(pstate,deepcopy.(initialstate)) 
    cΔy²,cLy²,cΔa²,cLa²= maxΔy^2,maxLy^2,maxΔa^2,maxLa^2
    nA                 = getndof(model,:A)
    La                 = Vector{𝕣 }(undef,nA   )
    Laa                = Matrix{𝕣 }(undef,nA,nA)
    Δy                 = Vector{𝕣1}(undef,length(state))
    y∂a                = Vector{𝕣2}(undef,length(state))
    Δy²,Ly²            = Vector{𝕣 }(undef,length(state)),Vector{𝕣}(undef,length(state))
    γ                  = γ0
    local facLyy
    local facLyys
    for iAiter          = 1:maxAiter
        verbose && @printf "    A-iteration %3d\n" iAiter
        La            .= 0
        Laa           .= 0
        for step     ∈ eachindex(state)
            for iYiter = 1:maxYiter
                assemble!(out1,asm1,dis,model,state[step], γ,(dbg...,solver=:StaticXUA,step=step,iYiter=iYiter))
                try if iAiter==1 && step==1 && iYiter==1
                    facLyys = lu(out1.Lyy) 
                else
                    lu!(facLyys,out1.Lyy) 
                end catch; muscadeerror(@sprintf("Incremental Y-solution failed at step=%i, iAiter=%i, iYiter",step,iAiter,iYiter)) end
                Δy[ step]  = facLyys\out1.Ly
                decrement!(state[step],0,Δy[ step],Ydofgr)
                Δy²s,Ly²s = sum(Δy[step].^2),sum(out2.Ly.^2)
                if Δy²s≤cΔy² && Ly²s≤cLy² 
                    verbose && @printf "        step % i Y-converged in %3d Y-iterations:   |ΔY|=%7.1e  |∇L/∂Y|=%7.1e\n" step iYiter √(Δy²s) √(Ly²s)
                    break#out of iYiter
                end
                iYiter==maxYiter && muscadeerror(@sprintf("no Y-convergence after %3d Y-iterations. |ΔY|=%7.1e |Ly|=%7.1e\n",iYiter,√(Δy²s),√(Ly²s)))
            end
            assemble!(out2,asm2,dis,model,state[step], γ,(dbg...,solver=:StaticXUA,step=step,iAiter=iAiter))
            try if iAiter==1 && step==1
                facLyy = lu(out2.Lyy) 
            else
                lu!(facLyy,out2.Lyy)
            end catch; muscadeerror(@sprintf("Incremental solution failed at step=%i, iAiter=%i",step,iAiter));end
            Δy[ step]  = facLyy\out2.Ly  
            y∂a[step]  = facLyy\out2.Lya 
            La       .+= out2.La  - out2.Lya' * Δy[ step]  
            Laa      .+= out2.Laa - out2.Lya' * y∂a[step]
            Δy²[step],Ly²[step] = sum(Δy[step].^2),sum(out2.Ly.^2)
        end    
        Δa             = Laa\La 
        Δa²,La²        = sum(Δa.^2),sum(La.^2)
        for step       ∈ eachindex(state)
            ΔY         = Δy[step] - y∂a[step] * Δa
            decrement!(state[step],0,ΔY,Ydofgr)
            decrement!(state[step],0,Δa,Adofgr)
        end    
        γ             *= γfac1*exp(-(out2.α/γfac2)^2)
        if all(Δy².≤cΔy²) && all(Ly².≤cLy²) && Δa².≤cΔa² && La².≤cLa² 
            verbose && @printf "\n    StaticXUA converged in %3d A-iterations.\n" iAiter
            verbose && @printf "    maxₜ(|ΔY|)=%7.1e  maxₜ(|∇L/∂Y|)=%7.1e  |ΔA|=%7.1e  |∇L/∂A|=%7.1e\n" √(maximum(Δy²)) √(maximum(Ly²)) √(Δa²) √(La²)
            break#out of iAiter
        end
        iAiter==maxAiter && muscadeerror(@sprintf("no convergence after %3d A-iterations. |ΔY|=%7.1e |Ly|=%7.1e |ΔA|=%7.1e |La|=%7.1e\n",iAiter,√(maximum(Δy²)),√(maximum(Ly²)),√(Δa²),√(La²)))
    end
    return
end


