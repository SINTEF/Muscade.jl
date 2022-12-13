struct AllAdofs <: DofGroup 
    scale :: 𝕣1
end
function AllAdofs(model::Model,dis)
    scale  = Vector{𝕣}(undef,getndof(model,:A))
    for di ∈ dis
        for d ∈ di
            scale[d.index.A] = d.scale.A
        end
    end
    return AllAdofs(scale)
end
function decrement!(s::State,a::𝕣1,gr::AllAdofs) 
    s.A .-= a.*gr.scale
end
Base.getindex(s::State,gr::AllAdofs) = s.A./gr.scale # not used by solver
getndof(gr::AllAdofs) = length(gr.scale)

struct AllΛXUdofs <: DofGroup 
    Λscale :: 𝕣1
    Xscale :: 𝕣1
    Uscale :: 𝕣1
    nX     :: 𝕣
    nU     :: 𝕣
end
function AllΛXUdofs(model::Model,dis)
    nX     = getndof(model,:X)
    nU     = getndof(model,:U)
    Λscale = Vector{𝕣}(undef,nX)
    Xscale = Vector{𝕣}(undef,nX)
    Uscale = Vector{𝕣}(undef,nU)
    for di ∈ dis
        for d ∈ di
            Λscale[d.index.X] = d.scale.Λ
            Xscale[d.index.X] = d.scale.X
            Uscale[d.index.U] = d.scale.U
        end
    end
    return AllΛXUdofs(Λscale,Xscale,Uscale,nX,nU)
end
function decrement!(s::State,y::𝕣1,gr::AllΛXUdofs) 
    nX,nU = length(s.X[1]),length(s.U[1])
    s.Λ    .-= y[    1: nX   ].*gr.Λscale
    s.X[1] .-= y[ nX+1:2nX   ].*gr.Xscale
    s.U[1] .-= y[2nX+1:2nX+nU].*gr.Uscale
end
function Base.getindex(s::State,gr::AllΛXUdofs) # not used by solver
    nX,nU = length(s.X[1]),length(s.U[1])
    y = 𝕣1(undef,2nX+nU)
    y[    1: nX   ] = s.Λ    ./gr.Λscale
    y[ nX+1:2nX   ] = s.X[1] ./gr.Xscale
    y[2nX+1:2nX+nU] = s.U[1] ./gr.Uscale
    return y
end
getndof(gr::AllΛXUdofs) = length(2gr.nX+gr.nU)

struct ASMstaticΛXU_A <: Assembler 
    dis   :: Vector{Any}          # naïve version! 
    Ly    :: 𝕣1
    La    :: 𝕣1
    Lyy   :: SparseMatrixCSC{𝕣,𝕫} 
    Lya   :: 𝕣2   
    Laa   :: SparseMatrixCSC{𝕣,𝕫} 
    nX    :: 𝕫
end #  
spa(a,n) = sparse(Int64[],Int64[],Float64[],a,n)
function ASMstaticΛXU_A(model::Model,dis) 
    nX,nU,nA = getndof(model,(:X,:U,:A))
    return ASMstaticΛXU_A(dis,zeros(2nX+nU),zeros(nA),spa(2nX+nU,2nX+nU),𝕣2(undef,2nX+nU,nA),spa(nA,nA),nX)
end
function zero!(asm::ASMstaticΛXU_A)
    asm.Ly  .= 0
    asm.La  .= 0
    asm.Lyy .= 0
    asm.Lya .= 0
    asm.Laa .= 0
end
function addin!(asm::ASMstaticΛXU_A,scale,ieletyp,iele,eleobj,Λ,X,U,A, t,ε,dbg) 
    Nx,Nu,Na        = length(X[1]),length(U[1]),length(A) # in the element
    Nz              = 2Nx+Nu+Na                           # Z = [Y;A]=[Λ;X;U;A]       
    iλ,ix,iu,ia     = 1:Nx, Nx+1:2Nx, 2Nx+1:2Nx+Nu, 2Nx+Nu+1:2Nx+Nu+Na # index into element vectors ΔZ and Lz
    iy              = 1:2Nx+Nu           
    ΔZ              = variate{2,Nz}(δ{1,Nz,𝕣}())                 
    ΔΛ,ΔX,ΔU,ΔA     = view(ΔZ,iλ),view(ΔZ,ix),view(ΔZ,iu),view(ΔZ,ia) # TODO Static?

    L               = scaledlagrangian(scale,eleobj, Λ+ΔΛ, (∂0(X)+ΔX,),(∂0(U)+ΔU,),A+ΔA, t,ε,dbg)
    Lz,Lzz          = value_∂{1,Nz}(∂{2,Nz}(L)) 
    i               = asm.dis[ieletyp][iele].index
    iY              = Vector([i.X;i.X.+asm.nX;i.U.+2asm.nX]) # index of element dofs into model Ly
    iA              = Vector(i.A)                          # index of element dofs into model La
    asm.La[iA]     += Lz[ia]  
    asm.Ly[iY]     += Lz[iy]  
    asm.Laa[iA,iA] += Lzz[ia,ia]
    asm.Lya[iY,iA] += Lzz[iy,ia]
    asm.Lyy[iY,iY] += Lzz[iy,iy]
end
function staticXUA(pstate,dbg;model::Model,time::AbstractVector{𝕣},
    initial::State=State(model,Disassembler(model)),
    maxiter::ℤ=50,maxΔy::ℝ=1e-5,maxLy::ℝ=∞,maxΔa::ℝ=1e-5,maxLa::ℝ=∞,verbose::𝕓=true)

    verbose && @printf "    staticXUA solver\n\n"
    dis                = initial.dis
    asm                = ASMstaticΛXU_A(model,dis)
    Adofgr             = AllAdofs(  model,dis)
    Ydofgr             = AllΛXUdofs(model,dis)
    cΔy²,cLy²,cΔa²,cLa²= maxΔy^2,maxLy^2,maxΔa^2,maxLa^2
    state              = allocate(pstate,[settime(deepcopy(initial),t) for t∈time]) 
    nA                 = getndof(model,:A)
    La                 = Vector{𝕣}(undef,nA   )
    Laa                = Matrix{𝕣}(undef,nA,nA)
    Δy                 = Vector{𝕣1}(undef,length(time))
    y∂a                = Vector{𝕣2}(undef,length(time))
    Δy²,Ly²            = Vector{𝕣}(undef,length(time)),Vector{𝕣}(undef,length(time))
    for iiter          = 1:maxiter
        verbose && @printf "    A-iteration %3d\n" iiter
        La            .= 0
        Laa           .= 0
        for step     ∈ eachindex(time)
            assemble!(asm,dis,model,state[step], 0.,(dbg...,solver=:StaticXUA,step=step))
            Δy[ step]  = try asm.Lyy\asm.Ly  catch; muscadeerror(@sprintf("Incremental solution failed at step=%i, iiter=%i",step,iiter)) end
            y∂a[step]  = try asm.Lyy\asm.Lya catch; muscadeerror(@sprintf("Incremental solution failed at step=%i, iiter=%i",step,iiter)) end
            La       .+= asm.La  + asm.Lya' * Δy[ step]
            Laa      .+= asm.Laa + asm.Lya' * y∂a[step]
            Δy²[step],Ly²[step] = sum(Δy[step].^2),sum(asm.Ly.^2)
        end    
#        @show Laa
#        @show La
        Δa             = Laa\La 
#        @show Δa
        Δa²,La²        = sum(Δa.^2),sum(La.^2)
#        @show √Δa²
        for step       ∈ eachindex(time)
            ΔY         = Δy[step] - y∂a[step] * Δa
            decrement!(state[step],ΔY,Ydofgr)
            decrement!(state[step],Δa,Adofgr)
        end    
#        @show state[1].A
        if all(Δy².≤cΔy²) && all(Ly².≤cLy²) && Δa².≤cΔa² && La².≤cLa² 
            verbose && @printf "\n    StaticXUA converged in %3d A-iterations.\n" iiter
            verbose && @printf "    |Δy|=%7.1e |Ly|=%7.1e |Δa|=%7.1e |La|=%7.1e\n" √(maximum(Δy²)) √(maximum(Ly²)) √(Δa²) √(La²)
            break#out of the iiter loop
        end
        iiter==maxiter && muscadeerror(@sprintf("no convergence after %3d iterations. |Δy|=%7.1e |Ly|=%7.1e |Δa|=%7.1e |La|=%7.1e\n",iiter,√(maximum(Δy²)),√(maximum(Ly²)),√(Δa²),√(La²)))
    end
end


