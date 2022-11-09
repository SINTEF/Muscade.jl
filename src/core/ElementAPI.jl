abstract type AbstractElement  end

# TODO macros for neater syntax in element definition.  Someting like 
# @Xdofid   Ballast (nod=[1,1,1],typ=[:dx1,:dx2,:dx3])  
# @Udofid   Ballast (nod=[1,1,1],typ=[:dx1,:dx2,:dx3])  
# @Adofid   Ballast (nod=[2,2  ],typ=[:Δseadrag,:Δbuoyancy])  
# @espyable Ballast (X=(3,),)

# MUST be used by elements to unpack X and U.  Today, the various derivatives are packed into tuples.  Would we use Adiff tomorrow, allowing
# correct computation of e.g. Coriolis terms in beam elements?
∂n(Y,n) = n+1≤length(Y) ? Y[n+1] : zeros(eltype(Y[1]),size(Y[1])...)  # this implementation will be slow if zero is to be returned!
∂0(y)  = ∂n(y,0)
∂1(y)  = ∂n(y,1)
∂2(y)  = ∂n(y,2)

draw(axe,key,out, ::E,args...) where{E<:AbstractElement}    = nothing # by default, an element draws nothing

espyable(    ::Type{E}) where{E<:AbstractElement}  = ()
request2draw(::Type{E}) where{E<:AbstractElement}  = ()
doflist(     ::Type{E}) where{E<:AbstractElement}  = (inod=𝕫[],class=Symbol[],field=Symbol[])
### Not part of element API, not exported by Muscade

getndof(E::DataType)              = length(doflist(E).inod)
getnnod(E::DataType)              = maximum(doflist(E).inod) 
getdoflist(E::DataType)           = doflist(E).inod, doflist(E).class, doflist(E).field
getidof(E::DataType,class)        = findall(doflist(E).class.==class)  
getndof(E::DataType,class)        = length(getidof(E,class))  
getndofs(E::DataType)             = getndof(E,:X),getndof(E,:U),getndof(E,:A)

