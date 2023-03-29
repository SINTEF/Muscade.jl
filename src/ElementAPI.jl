abstract type AbstractElement  end

# MUST be used by elements to unpack X and U.  Today, the various derivatives are packed into tuples.  Would we use Adiff tomorrow, allowing
# correct computation of e.g. Coriolis terms in beam elements?
∂n(Y,n) = n+1≤lastindex(Y) ? Y[n+1] : zeros(eltype(Y[1]),size(Y[1])...)  # this implementation will be slow if zero is to be returned!
∂0(y)   = ∂n(y,0)
∂1(y)   = ∂n(y,1)
∂2(y)   = ∂n(y,2)

doflist(     ::Type{E}) where{E<:AbstractElement}  = muscadeerror(@sprintf("method 'Muscade.doflist' must be provided for elements of type '%s'\n",E))
### Not part of element API, not exported by Muscade

getnnod(E::DataType)              = maximum(doflist(E).inod) 
getdoflist(E::DataType)           = doflist(E).inod, doflist(E).class, doflist(E).field
getidof(E::DataType,class)        = findall(doflist(E).class.==class)  
getndof(E::DataType)              = length(doflist(E).inod)
getndof(E::DataType,class)        = length(getidof(E,class))  
getndof(E::DataType,class::Tuple) = ntuple(i->getndof(E,class[i]),length(class))



