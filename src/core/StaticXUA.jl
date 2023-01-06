
struct OUTstaticΛXU_A  
    Ly    :: 𝕣1
    La    :: 𝕣1
    Lyy   :: SparseMatrixCSC{𝕣,𝕫} 
    Lya   :: SparseMatrixCSC{𝕣,𝕫} 
    Laa   :: SparseMatrixCSC{𝕣,𝕫} # TODO make this a full matrix?
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
    Lya                = asmmat!(view(asm,4,:),view(asm,1,:),view(asm,2,:),nY,nA) 
    Laa                = asmmat!(view(asm,5,:),view(asm,2,:),view(asm,2,:),nA,nA)  
    out                = OUTstaticΛXU_A(Ly,La,Lyy,Lya,Laa)
    return out,asm,Adofgr,Ydofgr
end
function zero!(out::OUTstaticΛXU_A)
    out.Ly        .= 0
    out.La        .= 0
    out.Lyy.nzval .= 0
    out.Lya.nzval .= 0
    out.Laa.nzval .= 0
end
function addin!(out,asm,iele,scale,eleobj,Λ,X,U,A, t,ε,dbg) 
    Nx,Nu,Na        = length(X[1]),length(U[1]),length(A) # in the element
    Nz              = 2Nx+Nu+Na                           # Z = [Y;A]=[Λ;X;U;A]       
    ΔZ              = variate{2,Nz}(δ{1,Nz,𝕣}())                 
    iλ,ix,iu,ia     = gradientpartition(Nx,Nx,Nu,Na) # index into element vectors ΔZ and Lz
    ΔΛ,ΔX,ΔU,ΔA     = view(ΔZ,iλ),view(ΔZ,ix),view(ΔZ,iu),view(ΔZ,ia) # TODO Static?
    L               = scaledlagrangian(scale,eleobj, Λ+ΔΛ, (∂0(X)+ΔX,),(∂0(U)+ΔU,),A+ΔA, t,ε,dbg)
    Lz,Lzz          = value_∂{1,Nz}(∂{2,Nz}(L)) 
    iy              = 1:(2Nx+Nu)  
    @show 
    addin!(out.Ly       ,asm[1],iele,view(Lz,iy))
    addin!(out.La       ,asm[2],iele,view(Lz,ia))
    addin!(out.Lyy.nzval,asm[3],iele,view(Lzz,iy,iy))
    addin!(out.Lya.nzval,asm[4],iele,view(Lzz,iy,ia))
    addin!(out.Laa.nzval,asm[5],iele,view(Lzz,ia,ia))  
end

#------------------------------------

function staticXUA(pstate,dbg;model::Model,time::AbstractVector{𝕣},
    initial::State=State(model,Disassembler(model)),
    maxiter::ℤ=50,maxΔy::ℝ=1e-5,maxLy::ℝ=∞,maxΔa::ℝ=1e-5,maxLa::ℝ=∞,verbose::𝕓=true)

    verbose && @printf "    staticXUA solver\n\n"
    dis                = initial.dis
    out,asm,Adofgr,Ydofgr = prepare(OUTstaticΛXU_A,model,dis)
    cΔy²,cLy²,cΔa²,cLa²= maxΔy^2,maxLy^2,maxΔa^2,maxLa^2
    state              = allocate(pstate,[settime(deepcopy(initial),t) for t∈time]) 
    nA                 = getndof(model,:A)
    La                 = Vector{𝕣 }(undef,nA   )
    Laa                = Matrix{𝕣 }(undef,nA,nA)
    Δy                 = Vector{𝕣1}(undef,length(time))
    y∂a                = Vector{𝕣2}(undef,length(time))
    Δy²,Ly²            = Vector{𝕣 }(undef,length(time)),Vector{𝕣}(undef,length(time))
    for iiter          = 1:maxiter
        verbose && @printf "    A-iteration %3d\n" iiter
        La            .= 0
        Laa           .= 0
        for step     ∈ eachindex(time)
            assemble!(out,asm,dis,model,state[step], 0.,(dbg...,solver=:StaticXUA,step=step))
            Δy[ step]  = try out.Lyy\out.Ly  catch; muscadeerror(@sprintf("Incremental solution failed at step=%i, iiter=%i",step,iiter)) end
            y∂a[step]  = try out.Lyy\Matrix(out.Lya) catch; muscadeerror(@sprintf("Incremental solution failed at step=%i, iiter=%i",step,iiter)) end
            La       .+= out.La  - out.Lya' * Δy[ step]  # TODO is it correct to add out.La and out.Laa nstep times?
            Laa      .+= out.Laa - out.Lya' * y∂a[step]
            Δy²[step],Ly²[step] = sum(Δy[step].^2),sum(out.Ly.^2)
        end    
        Δa             = Laa\La 
        Δa²,La²        = sum(Δa.^2),sum(La.^2)
        for step       ∈ eachindex(time)
            ΔY         = Δy[step] - y∂a[step] * Δa
            decrement!(state[step],0,ΔY,Ydofgr)
            decrement!(state[step],0,Δa,Adofgr)
        end    
        if all(Δy².≤cΔy²) && all(Ly².≤cLy²) && Δa².≤cΔa² && La².≤cLa² 
            verbose && @printf "\n    StaticXUA converged in %3d A-iterations.\n" iiter
            verbose && @printf "    maxₜ(|ΔY|)=%7.1e  maxₜ(|∂L/∂Y|)=%7.1e  |ΔA|=%7.1e  |∂L/∂A|=%7.1e\n" √(maximum(Δy²)) √(maximum(Ly²)) √(Δa²) √(La²)
            break#out of the iiter loop
        end
        iiter==maxiter && muscadeerror(@sprintf("no convergence after %3d iterations. |Δy|=%7.1e |Ly|=%7.1e |Δa|=%7.1e |La|=%7.1e\n",iiter,√(maximum(Δy²)),√(maximum(Ly²)),√(Δa²),√(La²)))
    end
    return
end


