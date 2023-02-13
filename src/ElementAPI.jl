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
# an assembler that calls "lagrangian" will call the element's own method if implemented, or this one, which then calls the element's residual method
lagrangian(        eleobj::E,δX,X,U,A, t,γ,dbg) where{E<:AbstractElement} = δX ∘₁ residual(        eleobj,X,U,A, t,γ,dbg)
lagrangian(out,key,eleobj::E,δX,X,U,A, t,γ,dbg) where{E<:AbstractElement} = δX ∘₁ residual(out,key,eleobj,X,U,A, t,γ,dbg)
# an assembler that calls "residual" will call the element's own method if implemented, or this one, which then calls the element's lagrangian method
function residual(eleobj::E, X,U,A, t,γ,dbg) where{E<:AbstractElement} 
    P            = constants(∂0(X),∂0(U),A,t)
    Nx           = length(∂0(X))
    δX           = δ{P,Nx,𝕣}()   
    L            = lagrangian(eleobj,δX,X,U,A, t,γ,dbg)
    return ∂{P,Nx}(L)
end
# if an element implements neither lagrangian nor residual, the above code will flat-spin recursively

