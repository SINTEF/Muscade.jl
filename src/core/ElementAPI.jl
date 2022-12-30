abstract type AbstractElement  end

# MUST be used by elements to unpack X and U.  Today, the various derivatives are packed into tuples.  Would we use Adiff tomorrow, allowing
# correct computation of e.g. Coriolis terms in beam elements?
∂n(Y,n) = n+1≤length(Y) ? Y[n+1] : zeros(eltype(Y[1]),size(Y[1])...)  # this implementation will be slow if zero is to be returned!
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
getndof(E::DataType,class::Tuple) = (getndof(E,c) for c∈class)

####### Lagrangian from residual and residual from Lagrangian
# an assembler that calls "lagrangian" will call the element's own method if implemented, or this one, which then calls the element's residual method
lagrangian(        eleobj::E,δX,X,U,A, t,ε,dbg) where{E<:AbstractElement} = δX ∘₁ residual(        eleobj,X,U,A, t,ε,dbg)
lagrangian(out,key,eleobj::E,δX,X,U,A, t,ε,dbg) where{E<:AbstractElement} = δX ∘₁ residual(out,key,eleobj,X,U,A, t,ε,dbg)
# an assembler that calls "residual" will call the element's own method if implemented, or this one, which then calls the element's lagrangian method
function residual(eleobj::E, X,U,A, t,ε,dbg) where{E<:AbstractElement} 
    P            = constants(∂0(X),∂0(U),A,t)
    Nx           = length(∂0(X))
    δX           = δ{P,Nx,𝕣}()                        
    L            = lagrangian(eleobj,δX,X,U,A, t,ε,dbg)
    return ∂{P,Nx}(L)
end
# if an element implements neither lagrangian nor residual, the above code will flat-spin recursively

####### For testing: get all the gradients. 
function gradient(eleobj::E,Λ,X,U,A, t,ε,dbg) where{E<:AbstractElement}
    P            = constants(Λ,∂0(X),∂0(U),A,t)
    nX,nU,nA     = length(Λ),length(∂0(U)),length(A)
    N            = 2nX+nU+nA
    iΛ,iX,iU,iA  = (1:nX) , (1:nX) .+ nX , (1:nU) .+ 2nX , (1:nA) .+ (2nX+nU)  
    ΔY           = δ{P,N,𝕣}()                        
    L            = lagrangian(eleobj,Λ+ΔY[iΛ],(∂0(X)+ΔY[iX],),(∂0(U)+ΔY[iU],),A+ΔY[iA], t,ε,dbg)
    Ly           = ∂{P,N}(L)
    return (L=value{P}(L), Lλ=Ly[iΛ], Lx=Ly[iX], Lu=Ly[iU], La=Ly[iA])
end