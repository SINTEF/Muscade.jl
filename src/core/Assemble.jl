using ForwardDiff, DiffResults  # phasing this out, though!

# TODO consider Yota.jl
# TODO XOR
# TODO use Holy traits to dispatch on wether Elements have "residual" or "lagrangian" https://www.juliabloggers.com/the-emergent-features-of-julialang-part-ii-traits/ 
# TODO memory management in hessian and gradient, make Y static
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
        iX,iU,iA     =              𝕫1(undef,nX),𝕫1(undef,nU),𝕫1(undef,nA)  # tmp arrays
        sΛ,sX,sU,sA  = 𝕣1(undef,nX),𝕣1(undef,nX),𝕣1(undef,nU),𝕣1(undef,nA)
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
            dis[ieletyp][iele] = IS(XUA{𝕫,nX,nU,nA}(iX,iU,iA),ΛXUA{𝕣,nX,nU,nA}(sΛ,sX,sU,sA))
        end
    end
    return dis
end
function setscale!(model;scale=nothing,Λscale=nothing)  # scale = (X=(tx=10,rx=1),A=(drag=3.))
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
    for ieletyp ∈ eachindex(model.ele)
        eleobj  = model.eleobj[ieletyp]
        dis     = model.disassembler[ieletyp]
        assemblesequential!(asm,ieletyp,dis, eleobj,Λ,X,U,A, t,ε,(dbg...,ieletyp=ieletyp))
    end
end
function assemblesequential!(asm::Assembler,ieletyp,dis, eleobj,Λ,X,U,A, t,ε,dbg) 
    for iele  ∈ eachindex(eleobj)
        scale = dis[iele].scale
        index = dis[iele].index
        Λe    =       Λ[index.X].*scale.Λ               
        Xe    = Tuple(x[index.X].*scale.X for x∈X)
        Ue    = Tuple(u[index.U].*scale.U for u∈U)
        Ae    =       A[index.A].*scale.A
        addin!(asm,ieletyp,iele, eleobj[iele],Λe,Xe,Ue,Ae, t,ε,(dbg...,iele=iele))
    end
end


###### Automatic differentiation and adding in for single elements

# Lagrangian from residual and residual from Lagrangian
function lagrangian(ele::E,δX,X,U,A, t,ε,dbg) where{E<:AbstractElement} 
    TRe   = promote_type(eltype(δX),eltype(X[1]),eltype(U[1]),eltype(A))
    Re    = zeros(TRe,getndof(E,:X)) # TODO this allocates.  Can we allocate at compilation and zero at each call?
    residual(ele,Re,X,U,A, t,ε,dbg)
    return δX ∘₁ Re
end
function residual(ele::E, Re,X,U,A, t,ε,dbg) where{E<:AbstractElement} 
    P            = constants(∂0(X),∂0(U),A,t)
    N            = getndof(E,:X)
    δX           = δ{P,N,𝕣}()                        
    L            = lagrangian(ele,δX,X,U,A, t,ε,dbg)
    Re          .= ∂{P,N}(L)
end

# For the purpose of testing elements: get all the gradients
function gradient(ele::E,Λ,X,U,A, t,ε,dbg) where{E<:AbstractElement}
    P            = constants(Λ,∂0(X),∂0(U),A,t)
    nX,nU,nA     = getndofs(E) # TODO type stability?
    N            = 2nX+nU+nA
    iΛ,iX,iU,iA  = (1:nX) , (1:nX) .+ nX , (1:nU) .+ 2nX , (1:nA) .+ (2nX+nU)  
    ΔY           = δ{P,N,𝕣}()                        
    L            = lagrangian(ele,Λ+ΔY[iΛ],(∂0(X)+ΔY[iX],),(∂0(U)+ΔY[iU],),A+ΔY[iA], t,ε,dbg)
    Ly           = ∂{P,N}(L)
    return (L=value{P}(L), Lλ=Ly[iΛ], Lx=Ly[iX], Lu=Ly[iU], La=Ly[iA])
end

# Static X
struct ASMstaticX <: Assembler 
    dis   :: Vector{Any}          # naïve version! - just a shallow copy of model.disassembler
    R     :: 𝕣1
    K     :: SparseMatrixCSC{𝕣,𝕫} 
end # for good old static FEM 
function ASMstaticX(model::Model) 
    nX       = getndof(model,:X)
    return ASMstaticX(model.disassembler,zeros(nX),sparse(Int64[],Int64[],Float64[],nX,nX))
end
@generated function addin!(asm::ASMstaticX,ieletyp,iele,eleobj::E,Λ,X,U,A, t,ε,dbg)  where{E<:AbstractElement}
    Nx      = length(Λ) 
    ΔX      = δ{1,Nx,𝕣}()                 # NB: precedence==1, because input is not Adiff
    Re      = Vector{∂ℝ{1,Nx,𝕣}}(undef,Nx)  # BUG one memory - common for all CPU threads?
    return quote
        $Re        .= 0.
        residual(eleobj, $Re,(∂0(X)+$ΔX,),U,A, t,ε,dbg)
        i           = asm.dis[ieletyp][iele].index.X    # TODO not type stable!
        asm.R[i  ] += value{1}($Re)
        j           = Vector(i)                               #
        asm.K[j,j] += ∂{1,$Nx}($Re)                 # TODO very slow!
    end
end



