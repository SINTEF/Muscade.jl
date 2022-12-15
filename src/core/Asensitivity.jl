struct ASMsensitivityA <: Assembler 
    dis   :: Vector{Any}          # naïve version! 
    La    :: 𝕣1
    Lλx   :: SparseMatrixCSC{𝕣,𝕫}   
    Lλa   :: SparseMatrixCSC{𝕣,𝕫} 
    nX    :: 𝕫
    nA    :: 𝕫
end #  
function ASMsensitivityA(model::Model,dis) 
    nX,nA = getndof(model,(:X,:A))
    return ASMsensitivityA(dis,zeros(nA),spa(nX,nX),spa(nX,nA),nX,nA)
end
function zero!(asm::ASMsensitivityA)
    asm.La  .= 0
    asm.Lλa .= 0
    asm.Lλx .= 0
end
function addin!(asm::ASMsensitivityA,scale,ieletyp,iele,eleobj,Λ,X,U,A, t,ε,dbg) 
    Nx,Na           = getndof(typeof(eleobj),:X),getndof(typeof(eleobj),:A) # in the element
    Nz              = 2Nx+Na                                  
    iλ,ix,ia        = 1:Nx, Nx+1:2Nx ,2Nx+1:2Nx+Na
    ΔZ              = variate{2,Nz}(δ{1,Nz,𝕣}())                 
    ΔΛ,ΔX,ΔA        = view(ΔZ,iλ),view(ΔZ,ix),view(ΔZ,ia) # TODO Static?
    L               = scaledlagrangian(scale,eleobj, Λ+ΔΛ, (∂0(X)+ΔX,),U,A+ΔA, t,ε,dbg)
    Lz,Lzz          = value_∂{1,Nz}(∂{2,Nz}(L)) 
    i               = asm.dis[ieletyp][iele].index
    iΛ = iX         = Vector(i.X)
    iA              = Vector(i.A)                          # index of element dofs into model La
    asm.La[iA]     += Lz[ia]  
    asm.Lλx[iΛ,iX] += Lzz[iλ,ix]
    asm.Lλa[iΛ,iA] += Lzz[iλ,ia]
end

get(v::ℝ ,i) = v
get(v::ℝ1,i) = v[i]
function Asensitivity(pstate,dbg;model::Model,time::𝕣=0.,initial::State=State(model,Disassembler(model);time), Δa=1.,pJa=Ref{Any}(),verbose::𝕓)
    verbose && @printf "    Asensitivity\n"
    nA                 = getndof(model,:A)
    dis                = initial.dis
    asm                = ASMsensitivityA(model,dis)
    Xdofgr             = AllXdofs(model,dis)
    state              = allocate(pstate,[deepcopy(initial) for iA=1:nA]) 
    assemble!(asm,dis,model,initial, 0.,(dbg...,solver=:Asensitivity))
    for iA = 1:nA
        ΔR = Vector(asm.Lλa[:,iA]*get(Δa,iA))
        ΔX = try  asm.Lλx\ΔR catch; muscadeerror(@sprintf("Sensitivity failed for iA=%i",iA)) end
        decrement!(state[iA], ΔX, Xdofgr)
    end
    allocate(pJa,asm.La-initial.Λ∘₁asm.Lλa) 
    return
end