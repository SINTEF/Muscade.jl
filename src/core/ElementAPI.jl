abstract type Element  end

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
function lagrangian(el::eltyp,δX,X,U,A, t,ε,dbg) where{eltyp<:Element} 
    TRe   = promote_type(eltype(δX),eltype(X[1]),eltype(U[1]),eltype(A))
    Re    = zeros(TRe,neldof(el).X) # TODO this allocates.  Can we allocate at compilation and zero at each call?
    residual(el,Re,X,U,A, t,ε,dbg)
    return δX ∘₁ Re
end
residual( ::eltyp, R,X,U,A, t,ε,dbg) where{eltyp<:Element}  = muscadeerror(@sprintf "no method residual for %s" eltyp )

draw(axe,key,out, ::eltyp,args...) where{eltyp<:Element}    = nothing # by default, an element draws nothing

espyable(    ::Type{eltyp}) where{eltyp<:Element}  = ()
request2draw(::Type{eltyp}) where{eltyp<:Element}  = ()
Xdofid(      ::Type{eltyp}) where{eltyp<:Element}  = (nod=𝕫[],typ=Symbol[])
Udofid(      ::Type{eltyp}) where{eltyp<:Element}  = (nod=𝕫[],typ=Symbol[])
Adofid(      ::Type{eltyp}) where{eltyp<:Element}  = (nod=𝕫[],typ=Symbol[])

# convenience functions based on the above
dofid(      ::Type{eltyp}) where{eltyp<:Element}   = (X=Xdofid(eltyp),U=Udofid(eltyp),A=Adofid(eltyp))
neldof(     ::Type{eltyp}) where{eltyp<:Element}   = (X=length(Xdofid(eltyp).nod),U=length(Udofid(eltyp).nod),A=length(Adofid(eltyp).nod))


Xdofid(  ::eltyp) where{eltyp<:Element}   = Xdofid(eltyp)
Udofid(  ::eltyp) where{eltyp<:Element}   = Udofid(eltyp)
Adofid(  ::eltyp) where{eltyp<:Element}   = Adofid(eltyp)
dofid(   ::eltyp) where{eltyp<:Element}   =  dofid(eltyp)
neldof(  ::eltyp) where{eltyp<:Element}   = neldof(eltyp)
espyable(::eltyp) where{eltyp<:Element}   = espyable(eltyp)

