abstract type AbstractElement  end

# TODO macros for neater syntax in element definition.  Someting like 
# @Xdofid   Ballast (nod=[1,1,1],typ=[:dx1,:dx2,:dx3])  
# @Udofid   Ballast (nod=[1,1,1],typ=[:dx1,:dx2,:dx3])  
# @Adofid   Ballast (nod=[2,2  ],typ=[:Δseadrag,:Δbuoyancy])  
# @espyable Ballast (X=(3,),)

# MUST be used by elements to unpack X and U.  Today, the various derivatives are packed into tuples.  Would we use Adiff tomorrow, allowing
# correct computation of e.g. Coriolis terms in beam elements?
∂(Y::ℝ11,n) = n+1≤length(Y) ? Y[n+1] : zeros(eltype(Y[1]),size(Y[1])...)  # this implementation will be slow if zero is to be returned!
∂0(y)  = ∂(y,0)
∂1(y)  = ∂(y,1)
∂2(y)  = ∂(y,2)

# to be implemented by elements (or not)
function lagrangian(ele::E,δX,X,U,A, t,ε,dbg) where{E<:AbstractElement} 
    TRe   = promote_type(eltype(δX),eltype(X[1]),eltype(U[1]),eltype(A))
    Re    = zeros(TRe,getndof(E,:X)) # TODO this allocates.  Can we allocate at compilation and zero at each call?
    residual(ele,Re,X,U,A, t,ε,dbg)
    return δX ∘₁ Re
end
function residual(ele::E, R,X,U,A, t,ε,dbg) where{E<:AbstractElement} 
    δX           = zeros(getndof(E,:X))                   # TODO this allocates!
    closure(δX)  = lagrangian(ele,δX,X,U,A, t,ε,dbg)
    R           .= ForwardDiff.gradient(closure,δX)
end





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

