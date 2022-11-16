###### scaling functions
function scaledlagrangian(scale,eleobj::E,Λs,Xs,Us,As, t,ε,dbg) where{E<:AbstractElement}
    Λ     =       Λs.*scale.Λ                 
    X     = Tuple(xs.*scale.X for xs∈Xs)
    U     = Tuple(us.*scale.U for us∈Us)
    A     =       As.*scale.A
    return lagrangian(eleobj,Λe,Xe,Ue,Ae, t,ε,dbg)
end    
function scaledresidual(scale,eleobj::E, Xs,Us,As, t,ε,dbg) where{E<:AbstractElement} 
    X     = Tuple(xs.*scale.X for xs∈Xs)
    U     = Tuple(us.*scale.U for us∈Us)
    A     =       As.*scale.A
    Re    = scale.Λ .* residual(eleobj, X,U,A, t,ε,dbg) 
end

######## The disassembler
copies(n,a::T) where{T}    = NTuple{n,T}(deepcopy(a) for i∈1:n) # TODO move to Dialect.jl
# dis[ieletyp][iele].index.[  X|U|A]
# dis[ieletyp][iele].scale.[Λ|X|U|A]
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
function assemble!(asm::Assembler,dis,model,state,ε,dbg)
    zero!(asm)
    for ieletyp ∈ eachindex(model.eleobj)
        eleobj  = model.eleobj[ieletyp]
        assemblesequential!(asm,ieletyp,dis[ieletyp], eleobj,state,ε,(dbg...,ieletyp=ieletyp))
    end
end
function assemblesequential!(asm::Assembler,ieletyp,dis, eleobj,state,ε,dbg) 
    for iele  ∈ eachindex(eleobj)
        scale = dis[iele].scale  # TODO unnecessary replication of "scale": is identical over iele...
        index = dis[iele].index
        Λe    = state.Λ[index.X]                 
        Xe    = Tuple(x[index.X] for x∈state.X)
        Ue    = Tuple(u[index.U] for u∈state.U)
        Ae    = state.A[index.A]
        addin!(asm,scale,ieletyp,iele,eleobj[iele],Λe,Xe,Ue,Ae, state.t,ε,(dbg...,iele=iele))
    end
end

