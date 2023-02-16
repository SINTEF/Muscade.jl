abstract type AbstractElement  end

# MUST be used by elements to unpack X and U.  Today, the various derivatives are packed into tuples.  Would we use Adiff tomorrow, allowing
# correct computation of e.g. Coriolis terms in beam elements?
∂n(Y,n) = n+1≤lastindex(Y) ? Y[n+1] : zeros(eltype(Y[1]),size(Y[1])...)  # this implementation will be slow if zero is to be returned!
∂0(y)  = ∂n(y,0)
∂1(y)  = ∂n(y,1)
∂2(y)  = ∂n(y,2)

draw(axe,key,out, ::E,args...;kwargs...) where{E<:AbstractElement}    = nothing # by default, an element draws nothing

espyable(    ::Type{E}) where{E<:AbstractElement}  = ()
request2draw(::Type{E}) where{E<:AbstractElement}  = ()
doflist(     ::Type{E}) where{E<:AbstractElement}  = (inod=𝕫[],class=Symbol[],field=Symbol[])
### Not part of element API, not exported by Muscade

getnnod(E::DataType)              = maximum(doflist(E).inod) 
getdoflist(E::DataType)           = doflist(E).inod, doflist(E).class, doflist(E).field
getidof(E::DataType,class)        = findall(doflist(E).class.==class)  
getndof(E::DataType)              = length(doflist(E).inod)
getndof(E::DataType,class)        = length(getidof(E,class))  
getndof(E::DataType,class::Tuple) = ntuple(i->getndof(E,class[i]),length(class))

####### Lagrangian from residual and residual from Lagrangian
# function implemented(::Type{T}) where{T} 
#     return (Val{hasmethod(residual  ,(T,   NTuple,NTuple,𝕣1,𝕣,𝕣,NamedTuple))},
#             Val{hasmethod(lagrangian,(T,𝕣1,NTuple,NTuple,𝕣1,𝕣,𝕣,NamedTuple))})
# end
@generated function implemented(eleobj) 
    r = hasmethod(residual  ,(eleobj,   NTuple,NTuple,𝕣1,𝕣,𝕣,NamedTuple))
    l = hasmethod(lagrangian,(eleobj,𝕣1,NTuple,NTuple,𝕣1,𝕣,𝕣,NamedTuple))
    return :(Val{$r},Val{$l})
end

# if residual or lagrange outputs just one vector or number, this element does not implementinequality constraints, so append minγfac=0.
defminγfac(x::Union{Number,AbstractVector})               = x,0.
defminγfac(x::Tuple)                                      = x

getresidual(          ::Type{<:Val}     ,::Type{<:Val}     ,args...) = muscadeerror(args[end],"No method 'lagrangian' or 'residual' for this element")
getlagrangian(        ::Type{<:Val}     ,::Type{<:Val}     ,args...) = muscadeerror(args[end],"No method 'lagrangian' or 'residual' for this element")

# Go straight
getresidual(          ::Type{Val{true}} ,::Type{<:Val}     ,args...) = defminγfac(residual(  args...))
getlagrangian(        ::Type{<:Val}     ,::Type{Val{true}} ,args...) = defminγfac(lagrangian(args...))

# Swap
# TODO merge the function pairs into one with Julia 1.9
function getresidual(  ::Type{Val{false}},::Type{Val{true}} ,eleobj, X,U,A, t,γ,dbg)  
    P            = constants(∂0(X),∂0(U),A,t)
    Nx           = length(∂0(X))
    δX           = δ{P,Nx,𝕣}()   
    L,minγfac    = defminγfac(lagrangian(eleobj,δX,X,U,A, t,γ,dbg))
    return ∂{P,Nx}(L),minγfac
end
function getresidual(::Type{Val{false}},::Type{Val{true}} ,out,key,eleobj,X,U,A, t,γ,dbg)  
    P            = constants(∂0(X),∂0(U),A,t)
    Nx           = length(∂0(X))
    δX           = δ{P,Nx,𝕣}()   
    L,minγfac    = defminγfac(lagrangian(out,key,eleobj,δX,X,U,A, t,γ,dbg))
    return ∂{P,Nx}(L),minγfac
end
function getlagrangian(::Type{Val{true}} ,::Type{Val{false}},eleobj,δX,X,U,A, t,γ,dbg) 
    (R,minγfac) = defminγfac(residual(eleobj,X,U,A, t,γ,dbg))
    return (δX ∘₁ R,minγfac)
end
function getlagrangian(::Type{Val{true}} ,::Type{Val{false}},out,key,eleobj,δX,X,U,A, t,γ,dbg) 
    (R,minγfac) = defminγfac(residual(out,key,eleobj,X,U,A, t,γ,dbg))
    return (δX ∘₁ R,minγfac)
end

###### scaled functions

function scaledlagrangian(scale,eleobj::E,Λs,Xs::NTuple{Nxder},Us::NTuple{Nuder},As, t,γ,dbg) where{E<:AbstractElement,Nxder,Nuder}
    Λ     =       Λs.*scale.Λ                 
    X     = NTuple{Nxder}(xs.*scale.X for xs∈Xs)  
    U     = NTuple{Nuder}(us.*scale.U for us∈Us)
    A     =       As.*scale.A
    L,minγfac = getlagrangian(implemented(eleobj)...,eleobj,Λ,X,U,A, t,γ,dbg)
    hasnan(L) && muscadeerror((dbg...,eletype=E),"NaN in a Lagrangian or its partial derivatives")
    return L
end    
function scaledresidual(scale,eleobj::E, Xs::NTuple{Nxder},Us::NTuple{Nuder},As, t,γ,dbg) where{E<:AbstractElement,Nxder,Nuder} 
    X     = NTuple{Nxder}(xs.*scale.X for xs∈Xs)  
    U     = NTuple{Nuder}(us.*scale.U for us∈Us)
    A     =       As.*scale.A
    R,minγfac = getresidual(implemented(eleobj)...,eleobj, X,U,A, t,γ,dbg) 
    hasnan(R) && muscadeerror(dbg,"NaN in a residual or its partial derivatives")
    return R.*scale.Λ 
end
