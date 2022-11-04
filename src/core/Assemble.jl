using ForwardDiff, DiffResults

# TODO consider Yota.jl
# TODO XOR
# TODO use Holy traits to dispatch on wether Elements have "residual" or "lagrangian" https://www.juliabloggers.com/the-emergent-features-of-julialang-part-ii-traits/ 
# TODO memory management in hessian and gradient, make Y static
# TODO gradient and hessian receives scaled values. The closure unscales them. In that way, Ly and Lyy are correctly scaled. 
# TODO Solvers store scaled states? What is the convention? dofID? Are result scaled?  Maybe the storage is unscaled and compact, 
# with a solver dependent accessor provided for the user (NodalResults)



copies(n,a::T) where{T}    = NTuple{n,T}(deepcopy(a) for i∈1:n)
struct XUA{T,nX,nU,nA} 
    X::SVector{nX,T}
    U::SVector{nU,T}
    A::SVector{nA,T}
end
struct ΛXUA{T,nX,nU,nA} 
    Λ::SVector{nX,T}
    X::SVector{nX,T}
    U::SVector{nU,T}
    A::SVector{nA,T}
end
struct IS{nX,nU,nA} 
    index:: XUA{𝕫,nX,nU,nA}
    scale::ΛXUA{𝕣,nX,nU,nA}
end
function Disassembler(model::Model)
    neletyp          = length(model.eleobj)  
    dis              = Vector{Any}(undef,neletyp)
    for ieletyp      = 1:neletyp
        nele         = length(model.eleobj[ieletyp])  
        E            = eltype(model.eleobj[ieletyp])
        nX,nU,nA     = getndofs(E)
        dis[ieletyp] = Vector{IS{nX,nU,nA}}(undef,nele)
        iX           = 𝕫1(undef,nX)  # working arrays
        sX           = 𝕣1(undef,nX)  
        sΛ           = 𝕣1(undef,nX)  
        iU           = 𝕫1(undef,nU)
        sU           = 𝕣1(undef,nU)
        iA           = 𝕫1(undef,nA)
        sA           = 𝕣1(undef,nA)
        for iele     = 1:nele
            ixdof,iudof,iadof = 0,0,0
            for dofID         ∈ model.ele[ieletyp][iele].dofID
                doftyp        = getdoftyp(model,dofID)
                class,scale   = doftyp.class,doftyp.scale
                if     class == :X
                    ixdof    += 1
                    iX[ixdof] = dofID.idof  
                    sX[ixdof] = scale
                    sΛ[ixdof] = scale * model.Λscale
                elseif class == :U
                    iudof    += 1
                    iU[iudof] = dofID.idof
                    sU[iudof] = scale
                elseif class == :A
                    iadof    += 1
                    iA[iadof] = dofID.idof
                    sA[iadof] = scale
                else
                    muscadeerror("dof class must be :X,:U or :A")
                end
            end
            dis[ieletyp][iele] = IS(XUA{𝕫,nX,nU,nA}(iX,iU,iA),ΛXUA{𝕣,nX,nU,nA}(sΛ,sX,sU,sA))# IS{nX,nU,nA}(XUA(iX,iU,iA),ΛXUA(sΛ,sX,sU,sA))
        end
    end
    return dis
end
function setscaling!(model;scale=nothing,Λscale=nothing)  # scale = (X=(tx=10,rx=1),A=(drag=3.))
    if ~isnothing(scale)
        for doftyp ∈ model.doftyp
            if doftyp.class ∈ keys(scale) && doftyp.field ∈ keys(scale[doftyp.class])
                doftyp.scale = scale[doftyp.class][doftyp.field] # otherwise leave untouched
            end
        end
    end
    if ~isnothing(Λscale)
        model.Λscale = Λscale
    end
    model.disassembler = Disassembler(model) # 
end

abstract type Assembler end
function assemble!(asm::Assembler,model,Λ,X,U,A, t,ε,dbg)
    for ieletyp ∈ eachindex(model.eletyp)
        eleobj  = model.eleobj[ieletyp]
        dis     = model.disassembler[ieletyp]
        assemblesequential!(asm,ieletyp,dis, eleobj,Λ,X,U,A, t,ε,dbg)
    end
end
function assembleequential!(asm::Assembler,ieletyp,dis, eleobj,Λ,X,U,A, t,ε,dbg) 
    for iele  ∈ eachindex(eleobj)
        scale = dis[iele].scale
        index = dis[iele].index
        Λe    =       Λ[index.X].*scale.Λ               
        Xe    = Tuple(x[index.X].*scale.X for x∈X)
        Ue    = Tuple(u[index.U].*scale.U for u∈U)
        Ae    =       A[index.A].*scale.A
        addin!(asm,ieletyp,iele, eleobj[iele],Λe,Xe,Ue,Ae, t,ε,dbg)
    end
end


######
struct ASMstaticX <: Assembler 
    dis   :: Vector{Any}          # naïve version! - just a shallow copy of model.disassembler
    R     :: 𝕣1
    K     :: SparseMatrixCSC{𝕣,𝕫} 
end # for good old static FEM 
function ASMstaticX(model::Model) 
    nX = getndof(model,:X)
    return ASMstaticX(model.disassembler,zeros(nX),sparse(Int64[],Int64[],Float64[],nX,nX))
end
function addin!(asm::ASMstaticX,ieletyp,iele,eleobj,Λ,X,U,A, t,ε,dbg) 
    for iele ∈ eachindex(eleobj)
        closure(X)  = residual(eleobj, Re,X,U,A, t,ε,dbg)
        result      = DiffResults.HessianResult(X)
        result      = ForwardDiff.hessian!(result, closure, X)
        r           = DiffResults.value(result)
        r∂x         = DiffResults.gradient(result)
        i           = asm.dis[ieletyp][iele].index.X    # TODO not type stable!
        asm.R[i  ] += r
        asm.K[i,i] += r∂x                               # TODO very slow!
    end
end



abstract type ASMjointΛXAstatic  <: Assembler end # for XA
function hessian(::Type{ASMjointΛXAstatic}, ele::E,Λ,X,U,A, t,ε,dbg) where{E<:AbstractElement}
    nX,_,nA      = getndofs(E)
    iΛ,iX,iA     = (1:nX) , (1:nX) .+ nX ,  (1:nA) .+ (2nX) 
    closure(Y)   = lagrangian(ele,Y[iΛ],[Y[iX]],U,Y[iA], t,ε,dbg)
    Y            = vcat(Λ,∂0(X),A)
    result       = DiffResults.HessianResult(Y)
    result       = ForwardDiff.hessian!(result, closure, Y)
    return (L=DiffResults.value(result), Ly=DiffResults.gradient(result), Lyy=DiffResults.hessian(result))
end

# For the purpose of testing elements
abstract type ASMseverΛXUAstatic <: Assembler end 
function gradient(::Type{ASMseverΛXUAstatic}, ele::E,Λ,X,U,A, t,ε,dbg) where{E<:AbstractElement}
    nX,nU,nA     = getndofs(E)
    iΛ,iX,iU,iA  = (1:nX) , (1:nX) .+ nX , (1:nU) .+ 2nX , (1:nA) .+ (2nX+nU)         
    closure(Y)   = lagrangian(ele,Y[iΛ],[Y[iX]],[Y[iU]],Y[iA], t,ε,dbg)
    Y            = vcat(Λ,∂0(X),∂0(U),A)
    result       = DiffResults.GradientResult(Y)
    result       = ForwardDiff.gradient!(result,closure,Y)
    Ly           = DiffResults.gradient(result)
    return (L=DiffResults.value(result), Lλ=Ly[iΛ], Lx=Ly[iX], Lu=Ly[iU], La=Ly[iA])
end
