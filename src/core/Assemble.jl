using ForwardDiff, DiffResults

# TODO consider Yota.jl
# TODO XOR
# TODO use Holy traits to dispatch on wether Elements have "residual" or "lagrangian" https://www.juliabloggers.com/the-emergent-features-of-julialang-part-ii-traits/ 
# TODO memory management in hessian and gradient, make Y static
# TODO gradient and hessian receives scaled values. The closure unscales them. In that way, Ly and Lyy are correctly scaled. 
# TODO Solvers store scaled states? What is the convention? dofID? Are result scaled?  Maybe the storage is unscaled and compact, 
# with a solver dependent accessor provided for the user (NodalResults)


# REFACTORING
#
# dofID, for a vector must be a unique, sequential identifier _within_a_class_. In other words, dofs must be stored in model by class.  model.dof.A[dofID] accesses a dof.
# This change "show elements" and other accessors to model.dof.
# This will allow to store a state as an object with separate classes, and provide the user with the key to access data in the state. 
#  
# Redesign disassembler to account for the new form of dofID and to have the structure
# asm[eletypID][iele].dofID.  X|U|A[ieledof]   
# asm[eletypID][iele].scale.Λ|X|U|A[ieledof]


copies(n,a::T) where{T}    = NTuple{n,T}(deepcopy(a) for i∈1:n)
struct Disassembler
    iX  :: Vector{𝕫2}  
    iU  :: Vector{𝕫2}  
    iA  :: Vector{𝕫2}  
    sΛ  :: Vector{𝕫2}
    sX  :: Vector{𝕫2}
    sU  :: Vector{𝕫2}
    sA  :: Vector{𝕫2}
end

function Disassembler(model::Model)
    neletyp          = length(model.eleobj)  
    iX,iU,iA         = copies(3,Vector{𝕫2}(undef,neletyp))
    sΛ,sX,sU,sA      = copies(4,Vector{𝕣2}(undef,neletyp))
    for eletypID     = 1:neletyp
        nele         = length(model.eleobj[eletypID])  
        nX, nU, nA   = getndofs(eltype(model.eleobj[eletypID])) 
        iX[eletypID] = 𝕫2(undef,nX,nele)
        iU[eletypID] = 𝕫2(undef,nU,nele)
        iA[eletypID] = 𝕫2(undef,nA,nele)
        sΛ[eletypID] = 𝕣2(undef,nX,nele)
        sX[eletypID] = 𝕣2(undef,nX,nele)
        sU[eletypID] = 𝕣2(undef,nU,nele)
        sA[eletypID] = 𝕣2(undef,nA,nele)
    end
    for ele               ∈ model.ele
        dofID             = ele.dofID
        eletypID          = ele.eletypID
        iele              = ele.iele
        ixdof,iudof,iadof = 0,0,0
        for dofID         ∈ ele.dofID
            doftyp        = model.doftyp[model.dof[dofID].doftypID]
            class,scale   = doftyp.class,doftyp.scale
            if     class == :X
                ixdof    += 1
                iX[eletypID][ixdof,iele] = dofID
                sX[eletypID][ixdof,iele] = scale
                sΛ[eletypID][ixdof,iele] = scale*model.Λscale
            elseif class == :U
                iudof    += 1
                iU[eletypID][iudof,iele] = dofID
                sU[eletypID][iudof,iele] = scale
            elseif class == :A
                iadof    += 1
                iA[eletypID][iadof,iele] = dofID
                sA[eletypID][iadof,iele] = scale
            else
                muscadeerror("dof class must be :X,:U or :A")
            end
        end
    end
    return Disassembler(iX,iU,iA,sΛ,sX,sU,sA)
end
function finalize!(model;scale=nothing,Λscale=nothing)  # scale = (X=(tx=10,rx=1),A=(drag=3.))
    if ~isnothing(scale)
        for doftyp ∈ model.doftyp
            if doftyp.class ∈ keys(scale) && doftyp.field ∈ keys(scale[doftyp.class])
                doftyp.scale = scale[doftyp.class][doftyp.field]
            end
        end
    end
    if ~isnothing(Λscale)
        model.Λscale = Λscale
    end
    model.disassembler = Disassembler(model)
end

abstract type Assembler end
function assemble!(asm::Assembler,model,Λ,X,U,A, t,ε,dbg)
    for ieletyp ∈ eachindex(model.eletyp)
        eleobj  = model.eleobj[ieletyp]
        dis     = model.disassembler[ieletyp]
        assemblesequential!(asm,ieletyp,dis,eleobj,Λ,X,U,A, t,ε,dbg)
    end
end
function assembleequential!(asm::Assembler,ieletyp,dis,eleobj,Λ,X,U,A, t,ε,dbg) 
    for iele  ∈ eachindex(eleobj)
        scale = dis[iele].scale
        dofID = dis[iele].dofID
        Λe    =       Λ[dofID.Λ].*scale.Λ               
        Xe    = Tuple(x[dofID.X].*scale.X for x∈X)
        Ue    = Tuple(u[dofID.U].*scale.U for u∈U)
        Ae    =       A[dofID.A].*scale.A
        addin!(asm,ieletyp,iele, eleobj[iele],Λe,Xe,Ue,Ae, t,ε,dbg)
    end
end


######
struct ASMstaticX <: Assembler 
    dis :: Disassembler          # naïve version! - just a shallow copy of model.disassembler
    R   :: 𝕣1
    K   :: XXXXX
    dofgr:: DofGroup
end # for good old static FEM 
ASMstaticX(model::Model) = ASMstaticX(model.disassembler,zeros(getndof(model)),XXXX EMPTY SPARSE XXXX)
function addin!(asm::ASMstaticX,ieletyp,iele,eleobj,Λ,X,U,A, t,ε,dbg) 
    for iele ∈ eachindex(eleobj)
        r,r∂x       = residual()

        i           =     
        asm.R[i  ] += r
        asm.K[i,i] += r∂x
    end
end

abstract type ASMjointΛXAstatic  <: Assembler end # for XA
function hessian(::Type{JointΛXAstatic}, ele::E,Λ,X,U,A, t,ε,dbg) where{E<:AbstractElement}
    nX,_,nA      = getndofs(E)
    iΛ,iX,iA     = (1:nX) , (1:nX) .+ nX ,  (1:nA) .+ (2nX) 
    closure(Y)   = lagrangian(ele,Y[iΛ],[Y[iX]],U,Y[iA], t,ε,dbg)
    Y            = vcat(Λ,∂0(X),A)
    result       = DiffResults.HessianResult(Y)
    result       = ForwardDiff.hessian!(result, closure, Y)
    return (L=DiffResults.value(result), Ly=DiffResults.gradient(result), Lyy=DiffResults.hessian(result))
end


abstract type ASMseverΛXUAstatic <: Assembler end # For the purpose of testing elements
function gradient(::Type{SeverΛXUAstatic}, ele::E,Λ,X,U,A, t,ε,dbg) where{E<:AbstractElement}
    nX,nU,nA     = getndofs(E)
    iΛ,iX,iU,iA  = (1:nX) , (1:nX) .+ nX , (1:nU) .+ 2nX , (1:nA) .+ (2nX+nU)         
    closure(Y)   = lagrangian(ele,Y[iΛ],[Y[iX]],[Y[iU]],Y[iA], t,ε,dbg)
    Y            = vcat(Λ,∂0(X),∂0(U),A)
    result       = DiffResults.GradientResult(Y)
    result       = ForwardDiff.gradient!(result,closure,Y)
    Ly           = DiffResults.gradient(result)
    return (L=DiffResults.value(result), Lλ=Ly[iΛ], Lx=Ly[iX], Lu=Ly[iU], La=Ly[iA])
end
