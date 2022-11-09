using ForwardDiff, DiffResults  # phasing this out, though!

# TODO consider Yota.jl
# TODO XOR
# TODO use Holy traits to dispatch on wether Elements have "residual" or "lagrangian" https://www.juliabloggers.com/the-emergent-features-of-julialang-part-ii-traits/ 
# TODO memory management in hessian and gradient, make Y static
# TODO Solvers store scaled states? What is the convention? dofID? Are result scaled?  Maybe the storage is unscaled and compact, 
# with a solver dependent accessor provided for the user (NodalResults)

####### Lagrangian from residual and residual from Lagrangian
# an assembler that calls "lagrangian" will call the element's own method if implemented, or this one, which then calls the element's residual method
function lagrangian(ele::E,δX,X,U,A, t,ε,dbg) where{E<:AbstractElement} 
    TRe   = promote_type(eltype(δX),eltype(X[1]),eltype(U[1]),eltype(A))
    Re    = zeros(TRe,getndof(E,:X)) # TODO this allocates.  Can we allocate at compilation and zero at each call?
    residual(ele,Re,X,U,A, t,ε,dbg)
    return δX ∘₁ Re
end
# an assembler that calls "residual" will call the element's own method if implemented, or this one, which then calls the element's lagrangian method
function residual(ele::E, Re,X,U,A, t,ε,dbg) where{E<:AbstractElement} 
    P            = constants(∂0(X),∂0(U),A,t)
    Nx           = length(∂0(X))
    δX           = δ{P,Nx,𝕣}()                        
    L            = lagrangian(ele,δX,X,U,A, t,ε,dbg)
    Re          .= ∂{P,Nx}(L)
end
# if an element implements neither lagrangian nor residual, the above code will flat-spin recursively

####### For testing: get all the gradients. 
function gradient(ele::E,Λ,X,U,A, t,ε,dbg) where{E<:AbstractElement}
    P            = constants(Λ,∂0(X),∂0(U),A,t)
    nX,nU,nA     = length(Λ),length(∂0(U)),length(A)
    N            = 2nX+nU+nA
    iΛ,iX,iU,iA  = (1:nX) , (1:nX) .+ nX , (1:nU) .+ 2nX , (1:nA) .+ (2nX+nU)  
    ΔY           = δ{P,N,𝕣}()                        
    L            = lagrangian(ele,Λ+ΔY[iΛ],(∂0(X)+ΔY[iX],),(∂0(U)+ΔY[iU],),A+ΔY[iA], t,ε,dbg)
    Ly           = ∂{P,N}(L)
    return (L=value{P}(L), Lλ=Ly[iΛ], Lx=Ly[iX], Lu=Ly[iU], La=Ly[iA])
end

###### scaling functions
function scaledlagrangian(scale,ele::E,Λs,Xs,Us,As, t,ε,dbg) where{E<:AbstractElement}
    Λ     =       Λs.*scale.Λ                 
    X     = Tuple(xs.*scale.X for xs∈Xs)
    U     = Tuple(us.*scale.U for us∈Us)
    A     =       As.*scale.A
    return lagrangian(ele,Λe,Xe,Ue,Ae, t,ε,dbg)
end    
function scaledresidual(scale,ele::E, Re,Xs,Us,As, t,ε,dbg) where{E<:AbstractElement} 
    X     = Tuple(xs.*scale.X for xs∈Xs)
    U     = Tuple(us.*scale.U for us∈Us)
    A     =       As.*scale.A
    residual(ele, Re,X,U,A, t,ε,dbg)
    Re  .*= scale.Λ
end

######## The disassembler
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


######## Generic assembler

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
        scale = dis[iele].scale  # TODO unnecessary replication of "scale": is identical over iele...
        index = dis[iele].index
        Λe    =       Λ[index.X]                 
        Xe    = Tuple(x[index.X] for x∈X)
        Ue    = Tuple(u[index.U] for u∈U)
        Ae    =       A[index.A]
        addin!(asm,scale,ieletyp,iele,eleobj[iele],Λe,Xe,Ue,Ae, t,ε,(dbg...,iele=iele))
    end
end


##### specialised addin

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
length(::StaticArrays.StaticIndexing{StaticArraysCore.SVector{L, Int64}}) where{L} = L
@generated function addin!(asm::ASMstaticX,scale,ieletyp,iele,eleobj::E,Λ,X,U,A, t,ε,dbg)  where{E<:AbstractElement}
    Nx      = length(Λ) 
    ΔX      = δ{1,Nx,𝕣}()                 # NB: precedence==1, input must not Adiff
    Re      = Vector{∂ℝ{1,Nx,𝕣}}(undef,Nx)  # BUG one memory - common for all CPU threads?
    i       = Vector{𝕫         }(undef,Nx)  # BUG one memory - common for all CPU threads?
    return quote
        $Re          .= 0.
        scaledresidual(scale,eleobj, $Re,(∂0(X)+$ΔX,),U,A, t,ε,dbg)
        $i           .= asm.dis[ieletyp][iele].index.X    # TODO not type stable (X is SVector)!
        asm.R[$i   ] += value{1}($Re)
        asm.K[$i,$i] += ∂{1,$Nx}($Re)                     # TODO very slow!
    end
end



