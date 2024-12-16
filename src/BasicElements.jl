"""
    DofCost{Class,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,
        afield,Tcost,Tcostargs} <: AbstractElement

An element to apply costs on combinations of dofs.  

# Named arguments to the constructor
- `xinod::NTuple{Nx,𝕫}=()`       For each X-dof to enter `cost`, its element-node number.
- `xfield::NTuple{Nx,Symbol}=()` For each X-dof to enter `cost`, its field.
- `uinod::NTuple{Nu,𝕫}=()`       For each U-dof to enter `cost`, its element-node number.
- `ufield::NTuple{Nu,Symbol}=()` For each U-dof to enter `cost`, its field.
- `ainod::NTuple{Na,𝕫}=()`       For each A-dof to enter `cost`, its element-node number.
- `afield::NTuple{Na,Symbol}=()` For each A-dof to enter `cost`, its field.
- `class:Symbol`                 `:A` for cost on A-dofs only, `:I` ("instant") otherwise.
- `cost::Function`               if `class==:I`, `cost(X,U,A,t,costargs...)→ℝ`
                                 if `class==:A`, `cost(A,costargs...)→ℝ` 
                                 `X` and `U` are tuples (derivates of dofs...), and `∂0(X)`,`∂1(X)`,`∂2(X)` 
                                 must be used by `cost` to access the value and derivatives of `X` (resp. `U`) 
- `[costargs::NTuple=() or NamedTuple] of additional arguments passed to `cost``


# Requestable internal variables
- `cost`, the value of the cost.

# Example
```
ele1 = addelement!(model,DofCost,[nod1],xinod=(1,),field=(:tx1,),
       class=:I,cost=(X,U,A,t;X0)->(X[1]-X0)^2,costargs=(;X0=0.27)
```

See also: [`SingleDofCost`](@ref), [`ElementCost`](@ref), [`addelement!`](@ref)  
"""
struct DofCost{Class,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,Tcost,Tcostargs} <: AbstractElement
    cost     :: Tcost     
    costargs :: Tcostargs
end
function DofCost(nod::Vector{Node};xinod::NTuple{Nx,𝕫}=(),xfield::NTuple{Nx,Symbol}=(),
                                uinod::NTuple{Nu,𝕫}=(),ufield::NTuple{Nu,Symbol}=(),
                                ainod::NTuple{Na,𝕫}=(),afield::NTuple{Na,Symbol}=(),
                                class::Symbol=:I,cost::Function ,costargs=()) where{Nx,Nu,Na} # :I for "instantaneous" or "integrand" cost.
    (class==:A && (Nx>0||Nu>0)) && muscadeerror("Cost with Class==:A must have zero X-dofs and zero U-dofs") 
    return DofCost{class,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,typeof(cost),typeof(costargs)}(cost,costargs)
end
doflist(::Type{<:DofCost{Class,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield}}) where
                     {Class,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield} = 
   (inod =(xinod...           ,uinod...           ,ainod...           ), 
    class=(ntuple(i->:X,Nx)...,ntuple(i->:U,Nu)...,ntuple(i->:A,Na)...), 
    field=(xfield...          ,ufield...          ,afield...          ) )
@espy function lagrangian(o::DofCost{:I,Nx,Nu,Na},Λ,X,U,A,t,SP,dbg) where{Nx,Nu,Na} 
    ☼cost = o.cost(X,U,A,t,o.costargs...)
    return cost,noFB
end
@espy function lagrangian(o::DofCost{:A,Nx,Nu,Na},Λ,X,U,A,t,SP,dbg) where{Nx,Nu,Na} 
    ☼cost = o.cost(    A  ,o.costargs...)
    return cost,noFB
end

"""
    ElementCost{Teleobj,Treq,Tcost,Tcostargs} <: AbstractElement

An element to apply costs on another "target" element's dofs and element-results.  
The target element must *not* be added separately to the model.  Instead, the 
`ElementType`, and the named arguments to the target element are provided
as input to the `ElementCost` constructor.

# Named arguments to the constructor
- `req`                 A request for element-results to be extracted from the target element, see [`@request`](@ref).
                        The request is formulated as if adressed directly to the target element (no "prefixing" with 
                        `eleres`, see "Requestable internal variables")
- `cost`                a cost function `cost(eleres,X,U,A,t,costargs...)→ℝ`
                        `X` and `U` are tuples (derivates of dofs...), and `∂0(X)`,`∂1(X)`,`∂2(X)` 
                        must be used by `cost` to access the value and derivatives of `X` (resp. `U`).
                        `X`, `U` and `A` are the degrees of freedom of the element `ElementType`.
- `[costargs::NTuple=() or NamedTuple] of additional arguments passed to `cost``
- `ElementType`         The named of the constructor for the relevant element 
- `elementkwargs`       A named tuple containing the named arguments of the `ElementType` constructor.     


# Requestable internal variables

From `ElementConstraint` one can request
- `cost`               The value of the cost

From the target element once can request
- `eleres(...)`        where `...` is the list of requestables from the target element.  It must be "prefixed" by 
                       `eleres` to prevent possible confusion with variables requestable from `ElementConstraint`.

# Example
```
@once cost(eleres,X,U,A,t;Fh0) = (eleres.Fh-Fh0)^2
ele1 = addelement!(model,ElementCost,[nod1];req=@request(Fh),
                   cost=cost,
                   costargs = (;Fh0=0.27)
                   ElementType=AnchorLine,
                   elementkwargs=(Λₘtop=[5.,0,0], xₘbot=[250.,0], L=290., buoyancy=-5e3))
```

See also: [`SingleDofCost`](@ref), [`DofCost`](@ref), [`@request`](@ref) 
"""
struct ElementCost{Teleobj,Treq,Tcost,Tcostargs} <: AbstractElement
    eleobj   :: Teleobj
    req      :: Treq
    cost     :: Tcost     
    costargs :: Tcostargs
end
function ElementCost(nod::Vector{Node};req,cost,costargs=(),ElementType,elementkwargs)
    eleobj   = ElementType(nod;elementkwargs...)
    return ElementCost(eleobj,(eleres=req,),cost,costargs)
end
doflist( ::Type{<:ElementCost{Teleobj}}) where{Teleobj} = doflist(Teleobj)
@espy function lagrangian(o::ElementCost, Λ,X,U,A,t,SP,dbg)
    req          = merge(o.req)
    L,FB,☼eleres = getlagrangian(o.eleobj,Λ,X,U,A,t,SP,(dbg...,via=ElementCost),req.eleres)
    ☼cost        = o.cost(eleres,X,U,A,t,o.costargs...) 
    return L+cost,FB
end    

"""
    SingleDofCost{Derivative,Class,Field,Tcost} <: AbstractElement

An element with a single node, for adding a cost to a given dof.  

# Named arguments to the constructor
- `class::Symbol`, either `:X`, `:U` or `:A`.
- `field::Symbol`.
- `cost::Function`, where 
    - `cost(x::ℝ,t::ℝ[,costargs...]) → ℝ` if `class` is `:X` or `:U`, and 
    - `cost(x::ℝ,    [,costargs...]) → ℝ` if `class` is `:A`.
- `[costargs::NTuple=() or NamedTuple] of additional arguments passed to `cost``
- `derivative::Int=0` 0, 1 or 2 - which time derivative of the dof enters the cost. 	    

# Requestable internal variables
- `cost`, the value of the cost.

# Example
```
using Muscade
model = Model(:TestModel)
node  = addnode!(model,𝕣[0,0])
e     = addelement!(model,SingleDofCost,[node];class=:X,field=:tx,
                    costargs=(3.,),cost=(x,t,three)->(x/three)^2)
```    

See also: [`DofCost`](@ref), [`ElementCost`](@ref)
"""
struct SingleDofCost <: AbstractElement end
function SingleDofCost(nod::Vector{Node};class::Symbol,field::Symbol,cost::Function,derivative=0::𝕫,costargs=()) 
    ∂=∂n(derivative)
    if     class==:X; DofCost(nod;xinod=(1,),xfield=(field,),class=:I,cost=(X,U,A,t,args...)->cost(∂(X)[1],t,args...),costargs)
    elseif class==:U; DofCost(nod;uinod=(1,),ufield=(field,),class=:I,cost=(X,U,A,t,args...)->cost(∂(U)[1],t,args...),costargs)
    elseif class==:A; DofCost(nod;ainod=(1,),afield=(field,),class=:A,cost=(    A,  args...)->cost(A[1]     ,args...),costargs)
    else              muscadeerror("'class' must be :X,:U or :A")
    end
end    

"""
    SingleUdof{XField,Ufield,Tcost} <: AbstractElement

An element that creates a Udof, and associates a cost to its value.
The value of the Udof is applied as a load to a Xdof on the same node.  

# Named arguments to the constructor
- `Xfield::Symbol`.
- `Ufield::Symbol`.
- `cost::Function`, where `cost(u::ℝ,t::ℝ[,costargs...]) → ℝ` 
- `[costargs::NTuple=() or NamedTuple] of additional arguments passed to `cost``

# Requestable internal variables
- `cost`, the value of the cost.

# Example
```
using Muscade
model = Model(:TestModel)
node  = addnode!(model,𝕣[0,0])
e     = addelement!(model,SingleUdof,[node];Xfield=:tx,Ufield=:utx,
                    costargs=(3.,),cost=(x,t,three)->(x/three)^2)
```    

See also: [`DofCost`](@ref), [`ElementCost`](@ref)
"""
struct SingleUdof{Tcost,Tcostargs,Xfield,Ufield} <: AbstractElement
    cost     :: Tcost     
    costargs :: Tcostargs
end
SingleUdof(nod::Vector{Node};Xfield::Symbol,Ufield::Symbol,cost::Tcost,costargs::Tcostargs=()) where{Tcost,Tcostargs} = SingleUdof{Tcost,Tcostargs,Xfield,Ufield}(cost,costargs)
doflist( ::Type{SingleUdof{Tcost,Tcostargs,Xfield,Ufield}}) where{Tcost,Tcostargs,Xfield,Ufield} = (inod=(1,1),class=(:X,:U),field=(Xfield,Ufield))
@espy function lagrangian(o::SingleUdof, Λ,X,U,A,t,SP,dbg)
    λ, u = Λ[1], ∂0(U)[1]
    return o.cost(u,t,o.costargs...)-λ*u,noFB
end    

#-------------------------------------------------

"""
    DofLoad{Tvalue,Field} <: AbstractElement

An element to apply a loading term to a single X-dof.  

# Named arguments to the constructor
- `field::Symbol`.
- `value::Function`, where `value(t::ℝ) → ℝ`.

# Requestable internal variables
- `F`, the value of the load.

# Examples
```
using Muscade
model = Model(:TestModel)
node  = addnode!(model,𝕣[0,0])
e     = addelement!(model,DofLoad,[node];field=:tx,value=t->3t-1)
```    

See also: [`Hold`](@ref), [`DofCost`](@ref)  
"""
struct DofLoad{Field,Tvalue,Targs} <: AbstractElement 
    value      :: Tvalue # Function
    args       :: Targs
end
DofLoad(nod::Vector{Node};field::Symbol,value::Tvalue,args...) where{Tvalue<:Function} = DofLoad{field,Tvalue,typeof(args)}(value,args)
doflist(::Type{<:DofLoad{Field}}) where{Field}=(inod=(1,), class=(:X,), field=(Field,))
@espy function residual(o::DofLoad, X,U,A,t,SP,dbg) 
    ☼F = o.value(t,o.args...)
    return SVector{1}(-F),noFB
end
#-------------------------------------------------

#McCormick(a,b)= α->a*exp(-(α/b)^2)            # provided as input to solvers, used by their Addin

S(λ,g,γ) = g*λ-γ # complementary slackness 

KKT(λ::𝕣        ,g::𝕣         ,γ::𝕣)                 = 0 # A pseudo-potential with strange derivatives
KKT(λ::∂ℝ{P,N,R},g::∂ℝ{P,N,R},γ::𝕣) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(0, λ.x*g.dx + S(λ.x,g.x,γ)*λ.dx)
KKT(λ:: ℝ       ,g::∂ℝ{P,N,R},γ::𝕣) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(0, λ.x*g.dx                    )
KKT(λ:: 𝕣       ,g::∂ℝ{P,N,R},γ::𝕣) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(0, λ  *g.dx                    )
KKT(λ::∂ℝ{P,N,R},g:: ℝ       ,γ::𝕣) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(0,            S(λ.x,g.x,γ)*λ.dx)
function KKT(λ::∂ℝ{Pλ,Nλ,Rλ},g::∂ℝ{Pg,Ng,Rg},γ::𝕣) where{Pλ,Pg,Nλ,Ng,Rλ<:ℝ,Rg<:ℝ}
    if Pλ==Pg
        R = promote_type(Rλ,Rg)
        return ∂ℝ{Pλ,Nλ}(convert(R,KKT(λ.x,g.x,γ)),convert.(R,     λ.x*g.dx + S(λ.x,g.x,γ)*λ.dx))
    elseif Pλ> Pg
        R = promote_type(Rλ,typeof(b))
        return ∂ℝ{Pλ,Nλ}(convert(R,KKT(λ  ,g.x,γ)),convert.(R,     λ.x*g.dx                    ))
    else
        R = promote_type(typeof(a),Rg)
        return ∂ℝ{Pg,Ng}(convert(R,KKT(λ.x,g  ,γ)),convert.(R,                S(λ.x,g.x,γ)*λ.dx))
    end
end

#-------------------------------------------------

"""
    off(t) → :off

A function which for any value `t` returns the symbol `off`.  Usefull for specifying
the keyword argument `mode=off` in adding an element of type ``DofConstraint` to
a `Model`.

See also: [`DofConstraint`](@ref), [`ElementConstraint`](@ref), [`equal`](@ref), [`positive`](@ref)
"""
off(t)     = :off
"""
    equal(t) → :equal

A function which for any value `t` returns the symbol `equal`.  Usefull for specifying
the keyword argument `mode=equal` in adding an element of type ``DofConstraint` to
a `Model`.

See also: [`DofConstraint`](@ref), [`ElementConstraint`](@ref), [`off`](@ref), [`positive`](@ref)
"""
equal(t)   = :equal
"""
    positive(t) → :positive

A function which for any value `t` returns the symbol `positive`.  Usefull for specifying
the keyword argument `mode=positive` in adding an element of type ``DofConstraint` to
a `Model`.

See also: [`DofConstraint`](@ref), [`ElementConstraint`](@ref), [`off`](@ref), [`equal`](@ref)
"""
positive(t) = :positive
"""
    DofConstraint{λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,
        afield,λinod,λfield,Tg,Tmode} <: AbstractElement

An element to apply physical/optimisation equality/inequality constraints on dofs. 

The constraints are holonomic, i.e. they apply to the values, not the time derivatives, of the involved dofs. 
This element is very general but not very user-friendly to construct, factory functions are provided for better useability. 
The sign convention is that the gap `g≥0` and the Lagrange multiplier `λ≥0`.

This element can generate three classes of constraints, depending on the input argument `λclass`.
- `λclass=:X` Physical constraint.  In mechanics, the Lagrange multiplier dof is a 
   generalized force, dual of the gap. The gap function must be of the form `gap(x,t,gargs...)`.
- `λclass=:U` Time varying optimisation constraint. For example: find `A`-parameters so that
   at all times, the response does not exceed a given criteria. The gap function must be of the form   
   `gap(x,u,a,t,gargs...)`.
- `λclass=:A` Time invariant optimisation constraint. For example: find `A`-parameters such that
   `A[1]+A[2]=gargs.somevalue`. The gap function must be of the form `gap(a,gargs...)`.

# Named arguments to the constructor
- `xinod::NTuple{Nx,𝕫}=()`       For each X-dof to be constrained, its element-node number.
- `xfield::NTuple{Nx,Symbol}=()` For each X-dof to be constrained, its field.
- `uinod::NTuple{Nu,𝕫}=()`       For each U-dof to be constrained, its element-node number.
- `ufield::NTuple{Nu,Symbol}=()` For each U-dof to be constrained, its field.
- `ainod::NTuple{Na,𝕫}=()`       For each A-dof to be constrained, its element-node number.
- `afield::NTuple{Na,Symbol}=()` For each A-dof to be constrained, its field.
- `λinod::𝕫`                     The element-node number of the Lagrange multiplier.
- `λclass::Symbol`               The class (`:X`,`:U` or `:A`) of the Lagrange multiplier. 
                                 See the explanation above for classes of constraints
- `λfield::Symbol`               The field of the Lagrange multiplier.
- `gap::Function`                The gap function.
- `gargs::NTuple`                Additional inputs to the gap function.
- `mode::Function`               where `mode(t::ℝ) -> Symbol`, with value `:equal`, 
                                 `:positive` or `:off` at any time. An `:off` constraint 
                                 will set the Lagrange multiplier to zero.

# Example
```jldoctest; output = false
using Muscade
model           = Model(:TestModel)
n1              = addnode!(model,𝕣[0]) 
e1              = addelement!(model,DofConstraint,[n1],xinod=(1,),xfield=(:t1,),
                              λinod=1, λclass=:X, λfield=:λ1,gap=(x,t)->x[1]+.1,
                              mode=positive)
e2              = addelement!(model,QuickFix  ,[n1],inod=(1,),field=(:t1,),
                              res=(x,u,a,t)->0.4x.+.08+.5x.^2)
initialstate    = initialize!(model)
setdof!(initialstate,1.;field=:λ1)
state           = solve(SweepX{0};initialstate,time=[0.],verbose=false) 
X               = state[1].X[1]

# output

2-element Vector{Float64}:
 -0.09999152289496528
  0.04500254313151041
```    

See also: [`Hold`](@ref), [`ElementConstraint`](@ref), [`off`](@ref), [`equal`](@ref), [`positive`](@ref)
"""
struct DofConstraint{λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield,Tg,Tgargs,Tmode} <: AbstractElement
    gap      :: Tg    # Class==:X gap(x,t,gargs...) ,Class==:U  gap(x,u,a,t,gargs...), Class==:A gap(a,gargs...) 
    gargs    :: Tgargs
    mode     :: Tmode # mode(t)->symbol, or Symbol for Aconstraints
end
function DofConstraint(nod::Vector{Node};xinod::NTuple{Nx,𝕫}=(),xfield::NTuple{Nx,Symbol}=(),
                                      uinod::NTuple{Nu,𝕫}=(),ufield::NTuple{Nu,Symbol}=(),
                                      ainod::NTuple{Na,𝕫}=(),afield::NTuple{Na,Symbol}=(),
                                      λinod::𝕫, λclass::Symbol, λfield::Symbol,
                                      gap::Function ,gargs=(),mode::Function) where{Nx,Nu,Na} 
    (λclass==:X && (Nu>0||Na>0)) && muscadeerror("Constraints with λclass=:X must have zero U-dofs and zero A-dofs") 
    (λclass==:A && (Nx>0||Nu>0)) && muscadeerror("Constraints with λclass=:A must have zero X-dofs and zero U-dofs") 
    return DofConstraint{λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield,
                       typeof(gap),typeof(gargs),typeof(mode)}(gap,gargs,mode)
end
doflist(::Type{<:DofConstraint{λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield}}) where
                              {λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield} = 
   (inod =(xinod...           ,uinod...           ,ainod...           ,λinod ), 
    class=(ntuple(i->:X,Nx)...,ntuple(i->:U,Nu)...,ntuple(i->:A,Na)...,λclass), 
    field=(xfield...          ,ufield...          ,afield...          ,λfield)) 
@espy function residual(o::DofConstraint{:X,Nx}, X,U,A,t,SP,dbg) where{Nx}
    γ          = default{:γ}(SP,0.) # γ=SP.γ - default 0
    P          = constants(∂0(X),t)
    m          = o.mode(t)
    x,☼λ       = ∂0(X)[SVector{Nx}(1:Nx)], ∂0(X)[Nx+1]   
    x∂         = variate{P,Nx}(x) 
    ☼gap,g∂x   = value_∂{P,Nx}(o.gap(x∂,t,o.gargs...)) 
    R = if     m==:equal;    SVector{Nx+1}((       -g∂x*λ)...,-gap       ) # - sign: λ interpreted as an external force on generalised dof g∂x
    elseif     m==:positive; SVector{Nx+1}((       -g∂x*λ)...,-S(λ,gap,γ)) 
    elseif     m==:off;      SVector{Nx+1}(ntuple(i->0,Nx)...,-λ         ) 
    end
    return R,(λ=λ,g=gap,mode=m)
end
@espy function lagrangian(o::DofConstraint{:U,Nx,Nu,Na}, Λ,X,U,A,t,SP,dbg) where{Nx,Nu,Na}
    γ          = default{:γ}(SP,0.)
    m          = o.mode(t)
    x,u,a,☼λ   = ∂0(X),∂0(U)[SVector{Nu}(1:Nu)],A,∂0(U)[Nu+1]
    ☼gap       = o.gap(x,u,a,t,o.gargs...)
    L = if     m==:equal;    -gap*λ         
    elseif     m==:positive; -KKT(λ,gap,γ)  
    elseif     m==:off;      -0.5λ^2         
    end
    return L,(λ=λ,g=gap,mode=m)
end
@espy function lagrangian(o::DofConstraint{:A,Nx,Nu,Na}, Λ,X,U,A,t,SP,dbg) where{Nx,Nu,Na}
    γ          = default{:γ}(SP,0.)
    m          = o.mode(t)
    a,☼λ       = A[SVector{Na}(1:Na)],A[    Na+1] 
    ☼gap       = o.gap(a,o.gargs...)
    L = if     m==:equal;    -gap*λ         
    elseif     m==:positive; -KKT(λ,gap,γ)  
    elseif     m==:off;      -0.5λ^2           
    end
    return L,(λ=λ,g=gap,mode=m)
end


#-------------------------------------------------

"""
    Hold <: AbstractElement

An element to set a single X-dof to zero.  

# Named arguments to the constructor
- `field::Symbol`. The field of the X-dof to constraint.
- `λfield::Symbol=Symbol(:λ,field)`. The field of the Lagrange multiplier.

# Example
```
using Muscade
model = Model(:TestModel)
node  = addnode!(model,𝕣[0,0])
e     = addelement!(model,Hold,[node];field=:tx)
```    

See also: [`DofConstraint`](@ref), [`DofLoad`](@ref), [`DofCost`](@ref) 
"""
struct Hold <: AbstractElement end  
function Hold(nod::Vector{Node};field::Symbol,λfield::Symbol=Symbol(:λ,field)) 
    gap(v,t)=v[1]
    return DofConstraint{:X     ,1, 0, 0, (1,),(field,),(),   (),    (),   (),    1,    λfield, typeof(gap),typeof(()),typeof(equal)}(gap,(),equal)
end

#-------------------------------------------------

"""
    ElementConstraint{Teleobj,λinod,λfield,Nu,Treq,Tg,Tgargs,Tmode} <: AbstractElement

An element to apply optimisation equality/inequality constraints on the element-results of 
another "target" element. The target element must *not* be added separately to the model.  Instead, the 
`ElementType`, and the named arguments to the target element are provided as input to the 
`ElementConstraint` constructor.

This element generates a time varying optimisation constraint. For example: find `A`-parameters so that
   at all times, the element-result von-Mises stress does not exceed a given value. 

The Lagrangian multiplier introduced by this optimisation constraint is of class :U   

# Named arguments to the constructor

- `λinod::𝕫`            The element-node number of the Lagrange multiplier.
- `λfield::Symbol`      The field of the Lagrange multiplier.
- `req`                 A request for element-results to be extracted from the target element, see [`@request`](@ref).
                        The request is formulated as if adressed directly to the target element (no "prefixing" with 
                        `eleres`, see "Requestable internal variables")
- `gₛ::𝕣=1.`             A scale for the gap.
- `λₛ::𝕣=1.`             A scale for the Lagrange multiplier.
- `gap`                 a gap function `gap(eleres,X,U,A,t,gargs...)→ℝ`
                        `X` and `U` are tuples (derivates of dofs...), and `∂0(X)`,`∂1(X)`,`∂2(X)` 
                        must be used by `cost` to access the value and derivatives of `X` (resp. `U`).
                        `X`, `U` and `A` are the degrees of freedom of the element `ElementType`.
- `gargs::NTuple`       Additional inputs to the gap function. 

- `mode::Function`      where `mode(t::ℝ) -> Symbol`, with value `:equal`, 
                        `:positive` or `:off` at any time. An `:off` constraint 
                        will set the Lagrange multiplier to zero.
- `ElementType`         The named of the constructor for the relevant element 
- `elementkwargs`       A named tuple containing the named arguments of the `ElementType` constructor.     

# Requestable internal variables

From `ElementConstraint` one can request
- `λ`                   The constraints Lagrange multiplier
- `gap`                 The constraints gap function

From the target element on can request
- `eleres(...)`         where `...` is the list of requestables from the target element.  It must be "prefixed" by 
                        `eleres` to prevent possible confusion with variables requestable from `ElementConstraint`.
δX
# Example

```
@once gap(eleres,X,U,A,t) = eleres.Fh^2
ele1 = addelement!(model,ElementCoonstraint,[nod1];req=@request(Fh),
                   gap,λinod=1,λfield=:λ,mode=equal, 
                   ElementType=AnchorLine,
                   elementkwargs=(Δxₘtop=[5.,0,0], xₘbot=[250.,0],L=290., buoyancy=-5e3))
```

See also: [`Hold`](@ref), [`DofConstraint`](@ref), [`off`](@ref), [`equal`](@ref), [`positive`](@ref), [`@request`](@ref)
"""
struct ElementConstraint{Teleobj,λinod,λfield,Nu,Treq,Tg,Tgargs,Tmode} <: AbstractElement
    eleobj   :: Teleobj
    req      :: Treq
    gap      :: Tg    
    gargs    :: Tgargs
    mode     :: Tmode 
end
function ElementConstraint(nod::Vector{Node};λinod::𝕫, λfield::Symbol,
    req,gap::Function,gargs=(;),mode::Function,ElementType,elementkwargs)
    eleobj   = ElementType(nod;elementkwargs...)
    Nu       = getndof(typeof(eleobj),:U)
    return ElementConstraint{typeof(eleobj),λinod,λfield,Nu,typeof((eleres=req,)),typeof(gap),typeof(gargs),typeof(mode)}(eleobj,(eleres=req,),gap,gargs,mode)
end
doflist( ::Type{<:ElementConstraint{Teleobj,λinod,λfield}}) where{Teleobj,λinod,λfield} =
    (inod =(doflist(Teleobj).inod... ,λinod),
     class=(doflist(Teleobj).class...,:U),
     field=(doflist(Teleobj).field...,λfield))
@espy function lagrangian(o::ElementConstraint{Teleobj,λinod,λfield,Nu}, Λ,X,U,A,t,SP,dbg) where{Teleobj,λinod,λfield,Nu} 
    req        = merge(o.req)
    γ          = default{:γ}(SP,0.)
    m          = o.mode(t)
    u          = getsomedofs(U,SVector{Nu}(1:Nu)) 
    ☼λ         = ∂0(U)[Nu+1]
    L,FB,☼eleres = getlagrangian(o.eleobj,Λ,X,u,A,t,SP,(dbg...,via=ElementConstraint),req.eleres)
    ☼gap       = o.gap(eleres,X,u,A,t,o.gargs...)
    L += if    m==:equal;    -gap*λ   
    elseif     m==:positive; -KKT(λ,gap,γ) 
    elseif     m==:off;      -0.5λ^2 
    end
    return L,(λ=λ,g=gap,mode=m)
end

#-------------------------------------------------

"""
    QuickFix <: AbstractElement

An element for creating simple elements with "one line" of code.  
Elements thus created have several limitations:
- physical elements with only X-dofs.
- only `R` can be espied.
The element is intended for testing.  Muscade-based applications should not include this in their API. 

# Named arguments to the constructor
- `inod::NTuple{Nx,𝕫}`. The element-node numbers of the X-dofs.
- `field::NTuple{Nx,Symbol}`. The fields of the X-dofs.
- `res::Function`, where `res(X::ℝ1,X′::ℝ1,X″::ℝ1,t::ℝ) → ℝ1`, the residual.

# Examples
A one-dimensional linear elastic spring with stiffness 2.
```jldoctest; output = false
using Muscade
model = Model(:TestModel)
node1  = addnode!(model,𝕣[0])
node2  = addnode!(model,𝕣[1])
e = addelement!(model,QuickFix,[node1,node2];inod=(1,2),field=(:tx1,:tx1),
                res=(X,X′,X″,t)->Svector{2}(2*(X[1]-X[2]),2*(X[2]-X[1])) )

# output

Muscade.EleID(1, 1)                       
```    
"""
struct QuickFix{Nx,inod,field,Tres} <: AbstractElement
    res        :: Tres    # R = res(X,X′,X″,t)
end
QuickFix(nod::Vector{Node};inod::NTuple{Nx,𝕫},field::NTuple{Nx,Symbol},res::Function) where{Nx} = QuickFix{Nx,inod,field,typeof(res)}(res)
doflist(::Type{<:QuickFix{Nx,inod,field}}) where{Nx,inod,field} = (inod =inod,class=ntuple(i->:X,Nx),field=(field)) 
@espy function residual(o::QuickFix, X,U,A, t,SP,dbg) 
    ☼R = o.res(∂0(X),∂1(X),∂2(X),t)
    return R,noFB
end

#-------------------------------------------------

"""
    Monitor <: AbstractElement

Debbuging tool: An element for for monitoring inputs to and outputs from
another element, during an analysis.     

Instead of adding the element to be monitored directly into the model,
add this element with the element to be monitored as argument.

Inputs and outs get @show'n to screen

# Named arguments to the constructor

- `ElementType`         The the type of element to be monitored-
- `trigger`            A function that takes `dbg` as an input and returns a boolean 
                        (`true`) to printout.
- `elementkwargs`       a `NamedTuple` containing the named arguments of the `ElementType` constructor.

"""
struct Monitor{Teleobj,Ttrigger} <: AbstractElement
    eleobj   :: Teleobj
    trigger  :: Ttrigger
end
function Monitor(nod::Vector{Node};ElementType,trigger::Function,elementkwargs)
    eleobj = ElementType(nod;elementkwargs...)
    return Monitor(eleobj,trigger)
end
doflist( ::Type{<:Monitor{Teleobj}}) where{Teleobj} = doflist(Teleobj)
@espy function lagrangian(o::Monitor{Teleobj}, Λ,X,U,A,t,SP,dbg)  where{Teleobj}
    L,FB = getlagrangian(o.eleobj,Λ,X,U,A,t,SP,(dbg...,via=Monitor))
    if o.trigger(dbg)
        @show dbg
        @show SP
        @show VALUE(Λ)
        @show VALUE(X[1])
        @show VALUE(U[1])
        @show VALUE(A)
        @show Teleobj
        @show doflist(Teleobj)
        @show L
    end
    return L,FB
end
