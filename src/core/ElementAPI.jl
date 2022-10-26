abstract type AbstractElement  end

# TODO macros for neater syntax in element definition.  Someting like 
# @Xdofid   Ballast (nod=[1,1,1],typ=[:dx1,:dx2,:dx3])  
# @Udofid   Ballast (nod=[1,1,1],typ=[:dx1,:dx2,:dx3])  
# @Adofid   Ballast (nod=[2,2  ],typ=[:Δseadrag,:Δbuoyancy])  
# @espyable Ballast (X=(3,),)

# MUST be used by elements to unpack X and U.  Today, the various derivatives are packed into tuples.  Would we use Adiff tomorrow, allowing
# correct computation of e.g. Coriolis terms in beam elements?
∂(Y,n) = n+1≤length(Y) ? Y[n+1] : zeros(eltype(Y[1]),size(Y[1])...)
∂0(y)  = ∂(y,0)
∂1(y)  = ∂(y,1)
∂2(y)  = ∂(y,2)

# to be implemented by elements (or not)
lagrangian(      ::eltyp,δX,X,U,A, t,ε,dbg) where{eltyp<:AbstractElement} = muscadeerror(@sprintf "no method lagrangian for %s" eltyp )

draw(axe,key,out, ::eltyp,args...) where{eltyp<:AbstractElement}    = nothing # by default, an element draws nothing

espyable(    ::Type{eltyp}) where{eltyp<:AbstractElement}  = ()
request2draw(::Type{eltyp}) where{eltyp<:AbstractElement}  = ()
Xdofid(      ::Type{eltyp}) where{eltyp<:AbstractElement}  = (nod=𝕫[],typ=Symbol[])
Udofid(      ::Type{eltyp}) where{eltyp<:AbstractElement}  = (nod=𝕫[],typ=Symbol[])
Adofid(      ::Type{eltyp}) where{eltyp<:AbstractElement}  = (nod=𝕫[],typ=Symbol[])

# convenience functions based on the above
dofid(      ::Type{eltyp}) where{eltyp<:AbstractElement}   = (X=Xdofid(eltyp),U=Udofid(eltyp),A=Adofid(eltyp))
neldof(     ::Type{eltyp}) where{eltyp<:AbstractElement}   = (X=length(Xdofid(eltyp).nod),U=length(Udofid(eltyp).nod),A=length(Adofid(eltyp).nod))

Xdofid(  ::eltyp) where{eltyp<:AbstractElement}   = Xdofid(eltyp)
Udofid(  ::eltyp) where{eltyp<:AbstractElement}   = Udofid(eltyp)
Adofid(  ::eltyp) where{eltyp<:AbstractElement}   = Adofid(eltyp)
dofid(   ::eltyp) where{eltyp<:AbstractElement}   =  dofid(eltyp)
neldof(  ::eltyp) where{eltyp<:AbstractElement}   = neldof(eltyp)
espyable(::eltyp) where{eltyp<:AbstractElement}   = espyable(eltyp)

