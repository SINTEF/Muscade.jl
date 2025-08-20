"""
    DofCost{Class,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,
        afield,Tcost,Tcostargs} <: AbstractElement

An element to apply costs on combinations of dofs. The cost is "per unit of time".
For once-off costs on A-dofs, see [`Acost`](@ref).

# Named arguments to the constructor
- `xinod::NTuple{Nx,𝕫}=()`       For each X-dof to enter `cost`, its element-node number.
- `xfield::NTuple{Nx,Symbol}=()` For each X-dof to enter `cost`, its field.
- `uinod::NTuple{Nu,𝕫}=()`       For each U-dof to enter `cost`, its element-node number.
- `ufield::NTuple{Nu,Symbol}=()` For each U-dof to enter `cost`, its field.
- `ainod::NTuple{Na,𝕫}=()`       For each A-dof to enter `cost`, its element-node number.
- `afield::NTuple{Na,Symbol}=()` For each A-dof to enter `cost`, its field.
- `cost::Function`               `cost(X,U,A,t,costargs...)→ℝ`
                                 `X` and `U` are tuples (derivates of dofs...), and `∂0(X)`,`∂1(X)`,`∂2(X)` 
                                 must be used by `cost` to access the value and derivatives of `X` (resp. `U`) 
- `costargs::NTuple=()` or `NamedTuple` of additional arguments passed to `cost``

# Requestable internal variables
- `cost`, the value of the cost.

# Example
```
ele1 = addelement!(model,DofCost,[nod1],xinod=(1,),xfield=(:tx1,),
       cost=(X,U,A,t;X0)->(X[1]-X0)^2,costargs=(;X0=0.27)
```

See also: [`Acost`](@ref), [`SingleDofCost`](@ref), [`ElementCost`](@ref), [`addelement!`](@ref)  
"""
struct DofCost{Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,Tcost,Tcostargs} <: AbstractElement
    cost     :: Tcost     
    costargs :: Tcostargs
end
function DofCost(nod::Vector{Node};xinod::NTuple{Nx,𝕫}=(),xfield::NTuple{Nx,Symbol}=(),
                                uinod::NTuple{Nu,𝕫}=(),ufield::NTuple{Nu,Symbol}=(),
                                ainod::NTuple{Na,𝕫}=(),afield::NTuple{Na,Symbol}=(),
                                cost::Function ,costargs=()) where{Nx,Nu,Na} 
    return DofCost{Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,typeof(cost),typeof(costargs)}(cost,costargs)
end
doflist(::Type{<:DofCost{Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield}}) where
                        {Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield} = 
   (inod =(xinod...           ,uinod...           ,ainod...           ), 
    class=(ntuple(i->:X,Nx)...,ntuple(i->:U,Nu)...,ntuple(i->:A,Na)...), 
    field=(xfield...          ,ufield...          ,afield...          ) )
@espy function lagrangian(o::DofCost{Nx,Nu,Na},Λ,X,U,A,t,SP,dbg) where{Nx,Nu,Na} 
    ☼cost = o.cost(X,U,A,t,o.costargs...)
    return cost,noFB
end


"""
    Acost{Na,ainod,afield,Tcost,Tcostargs} <: AbstractElement

An element to apply a once-off cost on a combination of A-dofs. 
For costs per unit of time on A-dofs (not recommended), see [`DofCost`](@ref).

# Named arguments to the constructor
- `inod::NTuple{Na,𝕫}=()`       For each A-dof to enter `cost`, its element-node number.
- `field::NTuple{Na,Symbol}=()` For each A-dof to enter `cost`, its field.
- `cost::Function`              `cost(A,costargs...)→ℝ`
- `costargs::NTuple=()` or `NamedTuple` of additional arguments passed to `cost``

# Requestable internal variables
- `cost`, the value of the cost.

# Example
```
ele1 = addelement!(model,Acost,[nod1],inod=(1,),field=(:EI,),
       cost=(A;A0)->(A[1]-A0)^2,costargs=(;A0=0.27)
```

See also:  [`SingleAcost`](@ref), [`DofCost`](@ref), [`SingleDofCost`](@ref), [`ElementCost`](@ref), [`addelement!`](@ref)  
"""
struct Acost{Na,inod,field,Tcost,Tcostargs} <: AbstractElement
    cost     :: Tcost     
    costargs :: Tcostargs
end
Acost(nod::Vector{Node};inod::NTuple{Na,𝕫}=(),field::NTuple{Na,Symbol}=(),cost::Function ,costargs=()) where{Na} = Acost{Na,inod,field,typeof(cost),typeof(costargs)}(cost,costargs)
doflist(::Type{<:Acost{Na,inod,field}}) where{Na,inod,field} = (inod =inod,class=ntuple(i->:A,Na),field=field)
@espy function lagrangian(o::Acost,Λ,X,U,A,t,SP,dbg)  
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
                        The request is formulated as if adressed directly to the target element.
- `cost`                A cost function `cost(eleres,X,U,A,t,costargs...)→ℝ`. 
                        `eleres` is the output of the above-mentionned request to the target element.     
                        `X` and `U` are tuples (derivates of dofs...), and `∂0(X)`,`∂1(X)`,`∂2(X)` 
                        must be used by `cost` to access the value and derivatives of `X` (resp. `U`).
                        `X`, `U` and `A` are the degrees of freedom of the element `ElementType`.
- `costargs::NTuple=()` or `NamedTuple` of additional arguments passed to `cost`.
- `ElementType`         The named of the constructor for the relevant element.
- `elementkwargs`       A named tuple containing the named arguments of the `ElementType` constructor.     

# Requestable internal variables

Not to be confused with the `req` provided as input to `addelement!` when adding an `ElementCost`, one can, after 
the analysis, request results from `ElementConstraint` 
- `cost`               The value of `cost(eleres,X,U,A,t,costargs...)`.
- `eleres(...)`        where `...` is the list of requestables wanted from the target element.  The "prefix"  
                       `eleres` is there to prevent possible confusion with variables requestable from `ElementConstraint`.  
                       For example `@request cost` would extract the value of the `ElementConstraint`'s function `cost`, while
                       `@request eleres(cost)` refers to the value of a variable called `cost` in the target element. 

# Example
```
@once cost(eleres,X,U,A,t;Fh0) = (eleres.Fh-Fh0)^2
ele1 = addelement!(model,ElementCost,[nod1];req=@request(Fh),
                   cost=cost,
                   costargs = (;Fh0=0.27)
                   ElementType=AnchorLine,
                   elementkwargs=(Λₘtop=[5.,0,0], xₘbot=[250.,0], L=290., buoyancy=-5e3))
```

See also: [`SingleDofCost`](@ref), [`DofCost`](@ref), [`@request`](@ref), [`@once`](@ref) 
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
    #@show typeof(X)
    req          = mergerequest(o.req)
    L,FB,☼eleres = getlagrangian(o.eleobj,Λ,X,U,A,t,SP,(dbg...,via=ElementCost),req.eleres)
    ☼cost        = o.cost(eleres,X,U,A,t,o.costargs...) 
    return L+cost,FB
end   
allocate_drawing(axis,eleobj::AbstractVector{Teleobj};kwargs...)                    where{Teleobj<:ElementCost} = allocate_drawing(axis,[eᵢ.eleobj for eᵢ∈eleobj];kwargs...)
update_drawing(  axis,eleobj::AbstractVector{Teleobj},oldmut,opt, Λ,X,U,A,t,SP,dbg) where{Teleobj<:ElementCost} = update_drawing(  axis,[eᵢ.eleobj for eᵢ∈eleobj],oldmut,opt, Λ,X,U,A,t,SP,dbg)  
display_drawing!(axis,::Type{<:ElementCost{Teleobj}},obs,opt)                       where{Teleobj}              = display_drawing!(axis,Teleobj,obs,opt)

"""
    SingleAcost <: AbstractElement

An element with a single node, for adding a once-off cost to a single A-dof.  

# Named arguments to the constructor
- `field::Symbol`.
- `cost::Function`, where 
    - `cost(a::ℝ,[,costargs...]) → ℝ` 
- `costargs::NTuple=()` or `NamedTuple` of additional arguments passed to `cost`

# Requestable internal variables
- `cost`, the value of the cost.

# Example
```
using Muscade
model = Model(:TestModel)
node  = addnode!(model,𝕣[0,0])
e     = addelement!(model,SingleAcost,[node];field=:EI,
                    costargs=(3.,),cost=(a,three)->(a/three)^2)
```    

See also: [`DofCost`](@ref), [`SingleDofCost`](@ref),  [`Acost`](@ref), [`ElementCost`](@ref)
"""
struct SingleAcost <: AbstractElement end
SingleAcost(nod::Vector{Node};field::Symbol,cost::Function,costargs=()) = 
    Acost(nod;inod=(1,),field=(field,),cost=(A,args...)->cost(A[1],args...),costargs)

"""
    SingleDofCost{Derivative,Class,Field,Tcost} <: AbstractElement

An element with a single node, for adding a cost to a given dof.  

# Named arguments to the constructor
- `class::Symbol`, either `:X` or `:U`.
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
    if     class==:X; DofCost(nod;xinod=(1,),xfield=(field,),cost=(X,U,A,t,args...)->cost(∂(X)[1],t,args...),costargs)
    elseif class==:U; DofCost(nod;uinod=(1,),ufield=(field,),cost=(X,U,A,t,args...)->cost(∂(U)[1],t,args...),costargs)
    else              muscadeerror("'class' must be :X or :U")
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
                        The request is formulated as if adressed directly to the target element.
- `gₛ::𝕣=1.`             A scale for the gap.
- `λₛ::𝕣=1.`             A scale for the Lagrange multiplier.
- `gap`                 A gap function `gap(eleres,X,U,A,t,gargs...)→ℝ`.
                        `eleres` is the output of the above-mentionned request to the target element. 
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

Not to be confused with the `req` provided as input to `addelement!` when adding an `ElementConstraint`, one can, after 
the analysis, request results from `ElementConstraint` 
- `λ`                   The constraints Lagrange multiplier
- `gap`                 The constraints gap function
- `eleres(...)`         where `...` is the list of requestables wanted from the target element.  The "prefix"  
                        `eleres` is there to prevent possible confusion with variables requestable from `ElementConstraint`.  
                        For example `@request gap` would extract the value of the `ElementConstraint`'s function `gap`, while
                        `@request eleres(gap)` refers to the value of a variable called `gap` in the target element. 

# Example

```
@once gap(eleres,X,U,A,t) = eleres.Fh^2
ele1 = addelement!(model,ElementCoonstraint,[nod1];req=@request(Fh),
                   gap,λinod=1,λfield=:λ,mode=equal, 
                   ElementType=AnchorLine,
                   elementkwargs=(Δxₘtop=[5.,0,0], xₘbot=[250.,0],L=290., buoyancy=-5e3))
```

See also: [`Hold`](@ref), [`DofConstraint`](@ref), [`off`](@ref), [`equal`](@ref), [`positive`](@ref), [`@request`](@ref), [`@once`](@ref)
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
    tmp = ElementConstraint{typeof(eleobj),λinod,λfield,Nu,typeof((eleres=req,)),typeof(gap),typeof(gargs),typeof(mode)}(eleobj,(eleres=req,),gap,gargs,mode)
    return tmp
end
doflist( ::Type{<:ElementConstraint{Teleobj,λinod,λfield}}) where{Teleobj,λinod,λfield} =
    (inod =(doflist(Teleobj).inod... ,λinod),
     class=(doflist(Teleobj).class...,:U),
     field=(doflist(Teleobj).field...,λfield))
@espy function lagrangian(o::ElementConstraint{Teleobj,λinod,λfield,Nu}, Λ,X,U,A,t,SP,dbg) where{Teleobj,λinod,λfield,Nu} 
    req        = mergerequest(o.req)
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
allocate_drawing(axis,eleobj::AbstractVector{Teleobj};kwargs...)                          where{Teleobj<:ElementConstraint} = allocate_drawing(axis,[eᵢ.eleobj for eᵢ∈eleobj];kwargs...)
update_drawing(  axis,eleobj::AbstractVector{Teleobj},oldmut,opt, Λ,X,U,A,t,SP,dbg)       where{Teleobj<:ElementConstraint} = update_drawing(  axis,[eᵢ.eleobj for eᵢ∈eleobj],oldmut,opt, Λ,X,U,A,t,SP,dbg)  
display_drawing!(axis,::Type{<:ElementConstraint{Teleobj}},obs,opt)                       where{Teleobj}                    = display_drawing!(axis,Teleobj,obs,opt)

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


