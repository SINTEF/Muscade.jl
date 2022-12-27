###### scaling functions
function scaledlagrangian(scale,eleobj::E,Λs,Xs,Us,As, t,ε,dbg) where{E<:AbstractElement}
    Λ     =       Λs.*scale.Λ                 
    X     = Tuple(xs.*scale.X for xs∈Xs)  # TODO Tuple is slow, not typestable
    U     = Tuple(us.*scale.U for us∈Us)
    A     =       As.*scale.A
    L     = lagrangian(eleobj,Λ,X,U,A, t,ε,dbg)
    hasnan(L) && muscadeerror(dbg,"NaN in a Lagrangian or its partial derivatives")
    return L
end    
function scaledresidual(scale,eleobj::E, Xs,Us,As, t,ε,dbg) where{E<:AbstractElement} 
    X     = Tuple(xs.*scale.X for xs∈Xs)
    U     = Tuple(us.*scale.U for us∈Us)
    A     =       As.*scale.A
    R     = scale.Λ .* residual(eleobj, X,U,A, t,ε,dbg) 
    hasnan(R) && muscadeerror(dbg,"NaN in a residual or its partial derivatives")
    return R
end

######## The disassembler
# dis[ieletyp].index.[iele][X|U|A]
# dis[ieletyp].scale.[Λ|X|U|A]
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
struct Disassembler{nX,nU,nA}
    index :: Vector{XUA{𝕫,nX,nU,nA}}
    scale :: ΛXUA{𝕣,nX,nU,nA}
end
function Disassembler(model::Model)
    neletyp          = length(model.eleobj)  
    dis              = Vector{Disassembler}(undef,neletyp)
    for ieletyp      = 1:neletyp
        nele         = length(model.eleobj[ieletyp])  
        E            = eltype(model.eleobj[ieletyp])
        nX,nU,nA     = getndofs(E)
        iX,iU,iA     = 𝕫1(undef,nX),𝕫1(undef,nU),𝕫1(undef,nA)  # tmp arrays
        index        = Vector{XUA{𝕫,nX,nU,nA}}(undef,nele)
        for iele     = 1:nele
            ixdof,iudof,iadof = 0,0,0
            for dofID         ∈ model.ele[ieletyp][iele].dofID
                doftyp        = getdoftyp(model,dofID)
                class         = doftyp.class
                if     class == :X
                    ixdof    += 1
                    iX[ixdof] = dofID.idof  
                elseif class == :U
                    iudof    += 1
                    iU[iudof] = dofID.idof
                elseif class == :A
                    iadof    += 1
                    iA[iadof] = dofID.idof
                else
                    muscadeerror("dof class must be :X,:U or :A")
                end
            end
            index[iele] = XUA{𝕫,nX,nU,nA}(iX,iU,iA)
        end
        sΛ,sX,sU,sA       = 𝕣1(undef,nX),𝕣1(undef,nX),𝕣1(undef,nU),𝕣1(undef,nA)
        ixdof,iudof,iadof = 0,0,0
        for dofID         ∈ model.ele[ieletyp][begin].dofID
            doftyp        = getdoftyp(model,dofID)
            class,scale   = doftyp.class,doftyp.scale
            if     class == :X
                ixdof    += 1
                sX[ixdof] = scale
                sΛ[ixdof] = scale * model.Λscale
            elseif class == :U
                iudof    += 1
                sU[iudof] = scale
            elseif class == :A
                iadof    += 1
                sA[iadof] = scale
            end
        end
        scale             = ΛXUA{𝕣,nX,nU,nA}(sΛ,sX,sU,sA)
        dis[ieletyp]      = Disassembler{nX,nU,nA}(index,scale)
    end
    return dis
end

######## Generic assembler

abstract type Assembler end
function assemble!(asm::Assembler,dis,model,state,ε,dbg)
    zero!(asm)
    for ieletyp ∈ eachindex(model.eleobj)
        eleobj  = model.eleobj[ieletyp]
        assemblesequential!(asm,dis[ieletyp], eleobj,state,ε,(dbg...,ieletyp=ieletyp))
    end
end
function assemblesequential!(asm::Assembler,dis, eleobj,state,ε,dbg) 
    scale     = dis.scale
    for iele  ∈ eachindex(eleobj)
        index = dis.index[iele]
        Λe    = state.Λ[index.X]                 
        Xe    = Tuple(x[index.X] for x∈state.X)
        Ue    = Tuple(u[index.U] for u∈state.U)
        Ae    = state.A[index.A]
        addin!(asm,index,scale,eleobj[iele],Λe,Xe,Ue,Ae, state.time,ε,(dbg...,iele=iele))
    end
end

