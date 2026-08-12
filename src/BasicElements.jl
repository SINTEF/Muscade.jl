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
- `cost::Functor`               `cost(X,U,A,t,costargs...)→ℝ`
                                 `X` and `U` are tuples (derivates of dofs...), and `∂0(X)`,`∂1(X)`,`∂2(X)` 
                                 must be used by `cost` to access the value and derivatives of `X` (resp. `U`) 
- `costargs::NTuple=()` or `NamedTuple` of additional arguments splatted when calling `cost`.

# Requestable internal variables
- `cost`, the value of the cost.

# Example
```
@functor with() xcost(X,U,A,t;X0) = (X[1]-X0)^2
ele1 = addelement!(model,DofCost,[nod1],xinod=(1,),xfield=(:tx1,),
       cost=xcost,costargs=(;X0=0.27)
```
Note that `X0` is within a `NamedTuple` in the call to `addelement!`, but not in the definition
of the `Functor` `xcost`.

See also: [`Acost`](@ref), [`SingleDofCost`](@ref), [`ElementCostAndConstraint`](@ref), [`addelement!`](@ref)  
"""
struct DofCost{Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,Tcost,Tcostargs} <: AbstractElement
    cost     :: Tcost     
    costargs :: Tcostargs
end
function DofCost(nod::Vector{Node};xinod::NTuple{Nx,𝕫}=(),xfield::NTuple{Nx,Symbol}=(),
                                   uinod::NTuple{Nu,𝕫}=(),ufield::NTuple{Nu,Symbol}=(),
                                   ainod::NTuple{Na,𝕫}=(),afield::NTuple{Na,Symbol}=(),
                                   cost::Tcost ,costargs::Tcostargs=()) where{Nx,Nu,Na,Tcost<:Functor,Tcostargs} 
    return DofCost{Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,Tcost,Tcostargs}(cost,costargs)
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
- `cost::Functor`              `cost(A,costargs...)→ℝ`
- `costargs::NTuple=()` or `NamedTuple` of additional arguments splatted when calling `cost`.

# Requestable internal variables
- `cost`, the value of the cost.

# Example
```
@functor with() acost(A;A0)=(A[1]-A0)^2
ele1 = addelement!(model,Acost,[nod1],inod=(1,),field=(:EI,),
       cost=acost,costargs=(;A0=0.27)
```
Note that `A0` is within a `NamedTuple` in the call to `addelement!`, but not in the definition
of the `Functor` `acost`.


See also:  [`SingleAcost`](@ref), [`DofCost`](@ref), [`SingleDofCost`](@ref), [`ElementCostAndConstraint`](@ref), [`addelement!`](@ref)  
"""
struct Acost{Na,inod,field,Tcost,Tcostargs} <: AbstractElement
    cost     :: Tcost     
    costargs :: Tcostargs
end
Acost(nod::Vector{Node};inod::NTuple{Na,𝕫}=(),field::NTuple{Na,Symbol}=(),cost::Tcost,costargs::Tcostargs=()) where{Na,Tcost<:Functor,Tcostargs} = Acost{Na,inod,field,Tcost,Tcostargs}(cost,costargs)
doflist(::Type{<:Acost{Na,inod,field}}) where{Na,inod,field} = (inod =inod,class=ntuple(i->:A,Na),field=field)
@espy function lagrangian(o::Acost,A,SP,dbg)  
    ☼cost = o.cost(    A  ,o.costargs...)
    return cost,noFB
end


"""
    SingleAcost <: AbstractElement

An element with a single node, for adding a once-off cost to a single A-dof.  

# Named arguments to the constructor
- `field::Symbol`.
- `cost::Functor`, where 
    - `cost(a::ℝ,[,costargs...]) → ℝ` 
- `costargs::NTuple=()` or `NamedTuple` of additional arguments splatted when calling `cost`.

# Requestable internal variables
- `cost`, the value of the cost.

# Example
```
using Muscade
model = Model(:TestModel)
node  = addnode!(model,𝕣[0,0])
@functor with() acost(a,three)=(a/three)^2
e     = addelement!(model,SingleAcost,[node];field=:EI,
                    costargs=(3.,),cost=acost)
```    

Note that `3.` is within a `NamedTuple` in the call to `addelement!`, but `three` is not, in the definition
of the `Functor` `acost`.

See also: [`DofCost`](@ref), [`SingleDofCost`](@ref),  [`Acost`](@ref), [`ElementCostAndConstraint`](@ref)
"""
struct SingleAcost <: AbstractElement end
#@noinline (o::Functor{:fA})(A::SVector,args...) = o.captured.cost(A[1],args...)
(o::Functor{:fA})(A,args...) = o.captured.cost(A[1],args...)
SingleAcost(nod::Vector{Node};field::Symbol,cost::Functor,costargs=()) = 
    Acost(nod;inod=(1,),field=(field,),cost=Functor{:fA}(;cost),costargs)

"""
    SingleDofCost{Derivative,Class,Field,Tcost} <: AbstractElement

An element with a single node, for adding a cost to a given dof.  

# Named arguments to the constructor
- `class::Symbol`, either `:X` or `:U`.
- `field::Symbol`.
- `cost::Functor`, where 
    - `cost(x::ℝ,t::ℝ[,costargs...]) → ℝ` if `class` is `:X` or `:U`, and 
    - `cost(x::ℝ,    [,costargs...]) → ℝ` if `class` is `:A`.
- `costargs::NTuple=()` or `NamedTuple` of additional arguments passed to `cost``
- `derivative::Int=0` 0, 1 or 2 - which time derivative of the dof enters the cost. 	    

# Requestable internal variables
- `cost`, the value of the cost.

# Example
```
using Muscade
model = Model(:TestModel)
node  = addnode!(model,𝕣[0,0])
@functor with() xcost(x,t,three)=(x/three)^2
e     = addelement!(model,SingleDofCost,[node];class=:X,field=:tx,
                    costargs=(3.,),cost=xcost)
```    
Note that `3.` is within a `NamedTuple` in the call to `addelement!`, but `three` is not, in the definition
of the `Functor` `xcost`.

See also: [`DofCost`](@ref), [`ElementCostAndConstraint`](@ref)
"""
struct SingleDofCost <: AbstractElement end
# 1) @functor macro is not called.  Instead, we have a method definition and a call to object constructor
# 2) Method defined for the functor is called by DofCost/lagrangian, extracts correct derivative and calls userfunctor(X,t,args...)
(o::Functor{:singleXcost})(X,U,A,t,args...) = o.captured.cost(∂n(o.captured.derivative)(X)[1],t,args...)
(o::Functor{:singleUcost})(X,U,A,t,args...) = o.captured.cost(∂n(o.captured.derivative)(U)[1],t,args...)
function SingleDofCost(nod::Vector{Node};class::Symbol,field::Symbol,cost::Functor,derivative=0::𝕫,costargs=()) 
    if     class==:X;   DofCost(nod;xinod=(1,),xfield=(field,),cost=Functor{:singleXcost}(;cost,derivative),costargs)
    elseif class==:U;   DofCost(nod;uinod=(1,),ufield=(field,),cost=Functor{:singleUcost}(;cost,derivative),costargs)
    else                muscadeerror("'class' must be :X or :U")
    end
end    

"""
    SingleUdof{XField,Ufield,Tcost} <: AbstractElement

An element that creates a Udof, and associates a cost to its value.
The value of the Udof is applied as a load to a Xdof on the same node.  

# Named arguments to the constructor
- `Xfield::Symbol`.
- `Ufield::Symbol`.
- `cost::Functor`, where `cost(u::ℝ,t::ℝ[,costargs...]) → ℝ` 
- `costargs::NTuple=()` or `NamedTuple` of additional arguments passed to `cost`.

# Requestable internal variables
- `cost`, the value of the cost.

# Example
```
using Muscade
model = Model(:TestModel)
node  = addnode!(model,𝕣[0,0])
@functor with() ucost(u,t,three)->(u/three)^2
e     = addelement!(model,SingleUdof,[node];Xfield=:tx,Ufield=:utx,
                    costargs=(3.,),cost=ucost)
```    

Note that `3.` is within a `NamedTuple` in the call to `addelement!`, but `three` is not, in the definition
of the `Functor` `ucost`.

See also: [`DofCost`](@ref), [`ElementCostAndConstraint`](@ref)
"""
struct SingleUdof{Tcost,Tcostargs,Xfield,Ufield} <: AbstractElement
    cost     :: Tcost     
    costargs :: Tcostargs
end
SingleUdof(nod::Vector{Node};Xfield::Symbol,Ufield::Symbol,cost::Tcost,costargs::Tcostargs=()) where{Tcost<:Functor,Tcostargs} = SingleUdof{Tcost,Tcostargs,Xfield,Ufield}(cost,costargs)
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
- `value::Functor`, where `value(t::ℝ) → ℝ`.
- `args::NTuple=()` or `NamedTuple` of additional arguments passed to `value`.

# Requestable internal variables
- `F`, the value of the load.

# Examples
```
using Muscade
model = Model(:TestModel)
node  = addnode!(model,𝕣[0,0])
@functor with(a=3,b=-1) load(t,c)=a*t^c+b
e     = addelement!(model,DofLoad,[node];field=:tx,value=load,args=(2.,))
```   

Note that `2.` is within a `Tuple` in the call to `addelement!`, but `c` is not, in the definition
of the `Functor` `load`.

See also: [`Hold`](@ref), [`DofCost`](@ref)  
"""
struct DofLoad{Field,Tvalue,Targs} <: AbstractElement 
    value      :: Tvalue # Functor
    args       :: Targs
end
DofLoad(nod::Vector{Node};field::Symbol,value::Tvalue,args::Targs=()) where{Tvalue<:Functor,Targs} = DofLoad{field,Tvalue,Targs}(value,args)
doflist(::Type{<:DofLoad{Field}}) where{Field}=(inod=(1,), class=(:X,), field=(Field,))
@espy function residual(o::DofLoad, X,U,A,t,SP,dbg) 
    ☼F = o.value(t,o.args...)
    return SVector{1}(-F),noFB
end
#-------------------------------------------------

#McCormick(a,b)= α->a*exp(-(α/b)^2)            # provided as input to solvers, used by their Addin

∂KKT(λ,g,γ) = λ.*g.-γ # complementary slacknesses 

KKT(λ::𝕣        ,g::𝕣         ,γ::𝕣)                 = 0 # A pseudo-potential with strange derivatives
KKT(λ::∂ℝ{P,N,R},g::∂ℝ{P,N,R},γ::𝕣) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(0, λ.x*g.dx + ∂KKT(λ.x,g.x,γ)*λ.dx)
KKT(λ:: ℝ       ,g::∂ℝ{P,N,R},γ::𝕣) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(0, λ.x*g.dx                       )
KKT(λ:: 𝕣       ,g::∂ℝ{P,N,R},γ::𝕣) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(0, λ  *g.dx                       )
KKT(λ::∂ℝ{P,N,R},g:: ℝ       ,γ::𝕣) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(0,            ∂KKT(λ.x,g.x,γ)*λ.dx)
function KKT(λ::∂ℝ{Pλ,Nλ,Rλ},g::∂ℝ{Pg,Ng,Rg},γ::𝕣) where{Pλ,Pg,Nλ,Ng,Rλ<:ℝ,Rg<:ℝ}
    if Pλ==Pg
        R = promote_type(Rλ,Rg)
        return ∂ℝ{Pλ,Nλ}(convert(R,KKT(λ.x,g.x,γ)),convert.(R,     λ.x*g.dx + ∂KKT(λ.x,g.x,γ)*λ.dx))
    elseif Pλ> Pg
        R = promote_type(Rλ,typeof(b))
        return ∂ℝ{Pλ,Nλ}(convert(R,KKT(λ  ,g.x,γ)),convert.(R,     λ.x*g.dx                       ))
    else
        R = promote_type(typeof(a),Rg)
        return ∂ℝ{Pg,Ng}(convert(R,KKT(λ.x,g  ,γ)),convert.(R,                ∂KKT(λ.x,g.x,γ)*λ.dx))
    end
end
KKT(λ::SVector{Nλ},g::SVector{Nλ},γ::𝕣) where{Nλ} = sum(KKT(λ[i],g[i],γ) for i=1:Nλ) 


#-------------------------------------------------

@functor with() off(t)     = :off
"""
    Muscade.off(t) → :off

A function which for any value `t` returns the symbol `off`.  Useful for specifying
the keyword argument `mode=off` in adding an element of type `DofConstraint` to
a `Model`.

See also: [`DofConstraint`](@ref), [`ElementCostAndConstraint`](@ref), [`equal`](@ref), [`positive`](@ref)
"""
off
@functor with() equal(t)   = :equal
"""
    Muscade.equal(t) → :equal

A function which for any value `t` returns the symbol `equal`.  Useful for specifying
the keyword argument `mode=equal` in adding an element of type `DofConstraint` to
a `Model`.

See also: [`DofConstraint`](@ref), [`ElementCostAndConstraint`](@ref), [`off`](@ref), [`positive`](@ref)
"""
equal
@functor with() positive(t) = :positive
"""
    Muscade.positive(t) → :positive

A function which for any value `t` returns the symbol `positive`.  Useful for specifying
the keyword argument `mode=positive` in adding an element of type `DofConstraint` to
a `Model`.

See also: [`DofConstraint`](@ref), [`ElementCostAndConstraint`](@ref), [`off`](@ref), [`equal`](@ref)
"""
positive
"""
    DofConstraint{λclass,Nλ,Nx,Nu,Na,
                  λinod,λfield, xinod,xfield, uinod,ufield, ainod,afield,
                  Tg,Tmode} <: AbstractElement

An element to apply physical/optimisation equality/inequality constraints on dofs. 

The constraints are holonomic, i.e. they apply to the values, not the time derivatives, of the involved dofs. 
This element is very general but not very user-friendly to construct, factory functions are provided for better useability. 
The sign convention is that each gap `g≥0` and each Lagrange multiplier `λ≥0`.

This element can generate three classes of constraints, depending on the input argument `λclass`.
- `λclass=:X` Physical constraint.  In mechanics, the Lagrange multiplier dof is a 
   generalized force, dual of the gap. The gap `Functor` must be of the form `gap(x,t,gargs...)`.
- `λclass=:U` Time varying optimisation constraint. For example: find `A`-parameters so that
   at all times, the response does not exceed a given criteria. The gap `Functor` must be of the form   
   `gap(x,u,a,t,gargs...)`.
- `λclass=:A` Time invariant optimisation constraint. For example: find `A`-parameters such that
   `A[1]+A[2]=gargs.somevalue`. The gap `Functor` must be of the form `gap(a,gargs...)`.
All constraints in one `DofConstraint` element are of the same `λclass`. The gap `Functor` must
return a `SVector{Nλ,𝕫}`.

# Named arguments to the constructor
- `λclass::Symbol`               The class (`:X`,`:U` or `:A`) of the Lagrange multipliers. 
- `λinod ::NTuple{Nλ,𝕫     }`    The element-nodes number of the Lagrange multipliers.
- `λfield::NTuple{Nλ,Symbol}`    The field of the Lagrange multipliers.
- `xinod ::NTuple{Nx,𝕫     }=()`       For each X-dof to be constrained, its element-node number.
- `xfield::NTuple{Nx,Symbol}=()` For each X-dof to be constrained, its field.
- `uinod ::NTuple{Nu,𝕫     }=()`       For each U-dof to be constrained, its element-node number.
- `ufield::NTuple{Nu,Symbol}=()` For each U-dof to be constrained, its field.
- `ainod ::NTuple{Na,𝕫     }=()`       For each A-dof to be constrained, its element-node number.
- `afield::NTuple{Na,Symbol}=()` For each A-dof to be constrained, its field.
- `λinod::𝕫`                     The element-node number of the Lagrange multiplier.
- `λclass::Symbol`               The class (`:X`,`:U` or `:A`) of the Lagrange multiplier. 
                                 See the explanation above for classes of constraints
- `λfield::Symbol`               The field of the Lagrange multiplier.
- `gap::Functor`                 The gap function.
- `gargs::NTuple=()` or `NamedTuple`  Additional inputs to the gap function.
- `mode::Functor`                where `mode(t::ℝ) -> Symbol`, with value `:equal`, 
                                 `:positive` or `:off` at any time. An `:off` constraint 
                                 will set the Lagrange multiplier to zero. Applies to all `Nλ`
                                 constraints.

# Example    TODO
```jldoctest; output = false
using Muscade,StaticArrays
model           = Model(:TestModel)
n1              = addnode!(model,𝕣[0])
@functor with() gap(x,    t)=SVector(x+.1)
@functor with() res(x,u,a,t)=0.4x.+.08+.5x.^2)
e1              = addelement!(model,DofConstraint,[n1],λclass=:X, 
                              λinod=(1,),λfield=(:λ1,), 
                              xinod=(1,),xfield=(:t1,),
                              gap=gap,
                              mode=Muscade.positive)
e2              = addelement!(model,Muscade.QuickFix  ,[n1],inod=(1,),field=(:t1,),
                              res=res
initialstate    = initialize!(model)
setdof!(initialstate,1.;field=:λ1)
state           = solve(SweepX{0};initialstate,time=[0.],verbose=false) 
X               = state[1].X[2]

# output

2-element Vector{Float64}:
 -0.09999152289496528
  0.04500254313151041
```    

See also: [`Hold`](@ref), [`ElementCostAndConstraint`](@ref), [`off`](@ref), [`equal`](@ref), [`positive`](@ref)
"""
struct DofConstraint{λclass,Nλ,Nx,Nu,Na,λinod,λfield,xinod,xfield,uinod,ufield,ainod,afield,Tg<:Functor,Tgargs,Tmode<:Functor} <: AbstractElement
    gap      :: Tg    # λclass==:X gap(x,t,gargs...) ,λclass==:U  gap(x,u,a,t,gargs...), λclass==:A gap(a,gargs...) 
    gargs    :: Tgargs
    mode     :: Tmode # mode(t)->symbol, or Symbol for Aconstraints
end
function DofConstraint(nod::Vector{Node};λclass::Symbol,
                                         λinod::NTuple{Nλ,𝕫}=(),λfield::NTuple{Nλ,Symbol}=(), # Nλ: number of constraints
                                         xinod::NTuple{Nx,𝕫}=(),xfield::NTuple{Nx,Symbol}=(),
                                         uinod::NTuple{Nu,𝕫}=(),ufield::NTuple{Nu,Symbol}=(),
                                         ainod::NTuple{Na,𝕫}=(),afield::NTuple{Na,Symbol}=(),
                                         gap::Tgap ,gargs::Tgargs=(),mode::Tmode) where{Nλ,Nx,Nu,Na,Tgap<:Functor,Tgargs,Tmode<:Functor} 
    (λclass==:X && (Nu>0||Na>0)) && muscadeerror("Constraints with λclass=:X must have zero U-dofs and zero A-dofs") 
    (λclass==:A && (Nx>0||Nu>0)) && muscadeerror("Constraints with λclass=:A must have zero X-dofs and zero U-dofs") 
    gaptyp = if     λclass==:X typeof(gap(SVector{Nx}(0. for i=1:Nx),0.,gargs...)) 
             elseif λclass==:U typeof(gap(SVector{Nx}(0. for i=1:Nx),
                                          SVector{Nu}(0. for i=1:Nu),
                                          SVector{Na}(0. for i=1:Na),0.,gargs...)) 
             elseif λclass==:A typeof(gap(SVector{Na}(0. for i=1:Na),   gargs...))    
             end 
    gaptyp <: SVector{Nλ} || muscadeerror(@sprintf("gap expected to return a SVector of length Nλ=%i, got %s.",Nλ,gaptyp))
    return DofConstraint{λclass,Nλ,Nx,Nu,Na,λinod,λfield,xinod,xfield,uinod,ufield,ainod,afield,
                       Tgap,Tgargs,Tmode}(gap,gargs,mode)
end
doflist(::Type{<:DofConstraint{λclass,Nλ,Nx,Nu,Na,λinod,λfield,xinod,xfield,uinod,ufield,ainod,afield}}) where
                              {λclass,Nλ,Nx,Nu,Na,λinod,λfield,xinod,xfield,uinod,ufield,ainod,afield} = 
   (inod =(λinod...               ,xinod...           ,uinod...           ,ainod...           ), 
    class=(ntuple(i->λclass,Nλ)...,ntuple(i->:X,Nx)...,ntuple(i->:U,Nu)...,ntuple(i->:A,Na)...), 
    field=(λfield...              ,xfield...          ,ufield...          ,afield...          )) 
@espy function residual(o::DofConstraint{:X,Nλ,Nx}, X,U,A,t,SP,dbg) where{Nλ,Nx}
    iλ         = SVector{Nλ}(1:Nλ)
    ix         = SVector{Nx}(Nλ+1:Nλ+Nx)
    γ          = default{:γ}(SP,0.) # γ=SP.γ - default 0
    m          = o.mode(t)
    ☼λ,x       = ∂0(X)[iλ], ∂0(X)[ix]    
    P,_,x∂     = variate(x,context=(X,U,A,t))
    ☼gap,g∂x   = value_∂{P,Nx}(o.gap(x∂,t,o.gargs...)) 
    R = if     m==:equal;    SVector{Nλ+Nx}(-gap...          ,(       -λ∘₁g∂x)...) # - sign: λ interpreted as an external force on generalised dof g∂x
    elseif     m==:positive; SVector{Nλ+Nx}(-∂KKT(λ,gap,γ)...,(       -λ∘₁g∂x)...) 
    elseif     m==:off;      SVector{Nλ+Nx}(-λ...            , ntuple(i->0,Nx)...) 
    end
    return R,(λ=λ,g=gap,mode=m)
end
@espy function lagrangian(o::DofConstraint{:U,Nλ,Nx,Nu}, Λ,X,U,A,t,SP,dbg) where{Nλ,Nx,Nu}
    iλ         = SVector{Nλ}(1:Nλ)
    iu         = SVector{Nu}(Nλ+1:Nλ+Nu)
    γ          = default{:γ}(SP,0.)
    m          = o.mode(t)
    ☼λ,x,u,a   = ∂0(U)[iλ], ∂0(X), ∂0(U)[iu], A
    ☼gap       = o.gap(x,u,a,t,o.gargs...)
    L = if     m==:equal;    -λ ∘₁ gap         
    elseif     m==:positive; -KKT(λ,gap,γ)  
    elseif     m==:off;      -0.5*(λ ∘₁ λ)         
    end
    return L,(λ=λ,g=gap,mode=m)
end
@espy function lagrangian(o::DofConstraint{:A,Nλ,Nx,Nu,Na}, Λ,X,U,A,t,SP,dbg) where{Nλ,Nx,Nu,Na}
    iλ         = SVector{Nλ}(1:Nλ)
    ia         = SVector{Na}(Nλ+1:Nλ+Na)
    γ          = default{:γ}(SP,0.)
    m          = o.mode(t)
    ☼λ,a       = A[iλ], A[ia] 
    ☼gap       = o.gap(a,o.gargs...)
    L = if     m==:equal;    -λ ∘₁ gap         
    elseif     m==:positive; -KKT(λ,gap,γ)  
    elseif     m==:off;      -0.5*(λ ∘₁ λ)           
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
@functor with() gap_for_hold_functor(v,t)=v[1]
Hold(nod::Vector{Node};field::Symbol,λfield::Symbol=Symbol(:λ,field)) = 
    DofConstraint{:X   ,1,1, 0, 0,  (1,),(λfield,), (1,),(field,), (),(), (),(),    typeof(gap_for_hold_functor),typeof(()),typeof(equal)}(gap_for_hold_functor,(),equal)


"""
    ElementCostAndConstraint{TargetElement,λinod,λfield,Nu,Treq,Tg,Tgargs,Tmode} <: AbstractElement

An element to apply a single costs and multiple physical and optimisation equality and inequality constraints 
on the element-results of another "target" element. The target element must *not* be added separately to the model.  
Instead, the `TargetElement`, and the named arguments to the target element are provided as input to the 
`ElementCostAndConstraint` constructor.

The Lagrangian multipliers introduced by physical and optimisation constraints are of class `:X` and `:U` respectively.   

!!! warning

    For solvers [`SweepX`](@ref) and [`SweepXA`](@ref), this element does not handle target elements with
    large number of dofs.  This will be taken care of in future versions.

# Named arguments to the constructor

Omit irrelevant arguments if no physical or optimisation constraints, or cost, are required.

- `λxinod::NTuple{Nλx,𝕫}`       The element-nodes number of the physical Lagrange multipliers.
- `λxfield::NTuple{Nλx,Symbol}` The fields of the physical Lagrange multipliers. 
- `λuinod::NTuple{Nλu,𝕫}`       The element-nodes number of the optimisation Lagrange multipliers.
- `λufield::NTuple{Nλu,Symbol}` The fields of the optimisation Lagrange multipliers. 
- `req`                         A request for element-results to be extracted from the target element, see [`@request`](@ref).
                                The request is formulated as if adressed directly to the target element.
- `gap`                         A gap function `gap(eleres,t,gargs...)` returning a `Tuple` of `ℝ`, with `Nλx` gaps of the physical constraints
                                followed by `Nλu` gaps of the optimisation constraints.        
                                `eleres` is the output of the above-mentionned request to the target element. 
- `gargs::NTuple=()` or `NamedTuple`    Additional inputs to the gap function. 
- `mode::Functor` where `mode(t::ℝ)`, returning a `Tuple` of `Symbol`s, with `Nλx` `:off`,`:equal` or `:positive` for the physical constraints
                                followed by `Nλu` similar `Symbol`s for the optimisation constraints.
- `cost`                        A cost `Functor` `cost(eleres,t,costargs...)→ℝ`. 
                                `eleres` is the output of the above-mentionned request to the target element.  
- `costargs::NTuple=()`         or `NamedTuple` of additional arguments splatted when calling `cost`.
- `TargetElement`               The named of the constructor for the target element.  This cannot be a `ElementCostAndConstraint`. 
- `elementkwargs`               A named tuple containing the named arguments of the `TargetElement` constructor.     

# Requestable internal variables

These are variables requestable from `ElementCostAndConstraint` to extract results after the analysis. 
This is distinct from the `req` input described above, which relate to element results to be extracted 
from the *target* element, to evaluate cost or gaps.

- `λ`                           The constraints Lagrange multiplier
- `gap`                         The constraints gap `Functor`
- `mode`                        The vector of modes for each constraint
- `eleres(...)`                 where `...` is the list of requestables wanted from the target element.  The "prefix"  
                                `eleres` is there to prevent possible confusion with variables requestable from `ElementCostAndConstraint`.  
                                For example `@request gap` would extract the value of the `ElementCostAndConstraint`'s varable `gap`, while
                                `@request eleres(gap)` refers to the value of a variable called `gap` in the target element. 


See also: [`DofCost`](@ref), [`Hold`](@ref), [`DofConstraint`](@ref), [`off`](@ref), [`equal`](@ref), [`positive`](@ref), [`@request`](@ref), [`@functor`](@ref)
"""
struct ElementCostAndConstraint{NSO,TargetElement,λxinod,λuinod,λxfield,λufield,Treq,Tgap,Tgargs,Tmode,Tcost,Tcostargs} <: AbstractElement
    eleobj   :: TargetElement
    req      :: Treq
    gap      :: Tgap    
    gapargs  :: Tgargs
    mode     :: Tmode 
    cost     :: Tcost     
    costargs :: Tcostargs
end

@inline function ElementCostAndConstraint(nod::Vector{Node};
    λxinod::NTuple{Nλx,𝕫}=(),λxfield::NTuple{Nλx,Symbol}=(),
    λuinod::NTuple{Nλu,𝕫}=(),λufield::NTuple{Nλu,Symbol}=(),
    req,gap::Tgap=nothing,gapargs::Tgargs=(),mode::Tmode=nothing,
    cost::Tcost=nothing,costargs::Tcostargs=(),
    TargetElement,elementkwargs=(;)
    ) where{Tgap<:Union{Functor,Nothing},Nλx,Nλu,Tgargs,Tmode<:Union{Functor,Nothing},Tcost<:Union{Functor,Nothing},Tcostargs}

    eleobj   = TargetElement(nod;elementkwargs...)
    NSO      = no_second_order(TargetElement) 

    @assert ~(TargetElement isa ElementCostAndConstraint) "ElementCostAndConstraint cannot be nested"
    ElementCostAndConstraint{NSO,typeof(eleobj),λxinod,λuinod,λxfield,λufield,typeof((eleres=req,)),Tgap,Tgargs,Tmode,Tcost,Tcostargs}(
                       eleobj,(eleres=req,),gap,gapargs,mode,cost,costargs)
end


function doflist( ::Type{<:ElementCostAndConstraint{NSO,TargetElement,λxinod,λuinod,λxfield,λufield}}) where{NSO,TargetElement,λxinod,λuinod,λxfield,λufield} 
    λxclass = ntuple(i->:X,length(λxinod)) 
    λuclass = ntuple(i->:U,length(λuinod)) 
    (inod =(λxinod... , λuinod... , doflist(TargetElement).inod... ),
     class=(λxclass..., λuclass..., doflist(TargetElement).class...),
     field=(λxfield..., λufield..., doflist(TargetElement).field...))
end
@inline function dofpartition(::Type{<:ElementCostAndConstraint{NSO,TargetElement,λxinod,λuinod}}) where{NSO,TargetElement,λxinod,λuinod}
    Nλx          = length(λxinod)
    Nλu          = length(λuinod)
    Nx           = getndof(TargetElement,:X) # refers to o.eleobj, not o
    Nu           = getndof(TargetElement,:U) # refers to o.eleobj, not o
    Na           = getndof(TargetElement,:A) # refers to o.eleobj, not o
    iλx          = SVector{Nλx,𝕫}(    1:    Nλx)
    iλu          = SVector{Nλu,𝕫}(    1:    Nλu)
    iλa          = SVector{0  ,𝕫}(             )
    ix           = SVector{Nx ,𝕫}(Nλx+1:Nλx+Nx )
    iu           = SVector{Nu ,𝕫}(Nλu+1:Nλu+Nu )
    ia           = SVector{Na ,𝕫}(    1:    Na )
    return iλx,ix,iλu,iu,iλa,ia
end

no_second_order(::Type{<:ElementCostAndConstraint{NSO}}) where{NSO} = NSO 

@espy function lagrangian(o::ElementCostAndConstraint{Val(false),TargetElement,λxinod,λuinod}, Λ,X,U,A,t,SP,dbg) where{TargetElement,λxinod,λuinod}
    iλx,ix,iλu,iu,_,_ = dofpartition(typeof(o))
    Nλx,Nλu           = length(iλx),length(iλu)
    req               = mergerequest(o.req)
    ♢λ                = SVector(∂0(X)[iλx]...,∂0(U)[iλu]...)
    Λe                = Λ[ix] 
    Xe                = getsomedofs(X,ix) 
    Ue                = getsomedofs(U,iu) 
    if Nλx > 0
        Pv,Nv,Xev     = variate(Xe,context=(Λ,U,A,t)) # TODO MAYBE: make it directional derivative in direction Λe ??? Careful, Λe is adiffed itself !!!
        Lv,FB,eleresv = getlagrangian(o.eleobj,Λe,Xev,Ue,A,t,SP,(dbg...,via=ElementCostAndConstraint),req.eleres)
        L             = value{Pv}(Lv)
        eleres        = value{Pv}(eleresv)
    else
        L,FB,eleres   = getlagrangian(o.eleobj,Λe,Xe,Ue,A,t,SP,(dbg...,via=ElementCostAndConstraint),req.eleres)
    end
    local ☼eleres
    if typeof(o.cost) ≠ Nothing    
        cost          = o.cost(eleres,t,o.costargs...) 
        L            += cost
    else
        cost          = 0.
    end
    local ☼cost
    if Nλx > 0 || Nλu > 0   
        γ            = default{:γ}(SP,0.)
        if Nλx >0
            gapv     = o.gap(eleresv,t,o.gapargs... )
            gap,∇ₓₑgap = value_∂{Pv,Nv}(gapv)
        else
            gap      = o.gap(eleres,t,o.gapargs... )
        end
        mode         = o.mode(t)
        for iλ ∈ iλx
            Λᵢ, λᵢ,mᵢ,gapᵢ,∇ₓₑgapᵢ = Λ[iλ],∂0(X)[iλ],mode[iλ],gap[iλ],∇ₓₑgap[iλ,:]
            L -=if   mᵢ==:equal;    Λe ∘₁ ∇ₓₑgapᵢ .* λᵢ  + Λᵢ *        gapᵢ    
            elseif   mᵢ==:positive; Λe ∘₁ ∇ₓₑgapᵢ .* λᵢ  + Λᵢ *∂KKT(λᵢ,gapᵢ,γ) 
            elseif   mᵢ==:off;                             Λᵢ.*     λᵢ 
            end
        end
        for iλ ∈ iλu.+Nλx
            λᵢ,mᵢ,gapᵢ = ∂0(U)[iλ],mode[iλ],gap[iλ]
            L -=if   mᵢ==:equal;    gapᵢ * λᵢ     
            elseif   mᵢ==:positive; KKT(λᵢ,gapᵢ,γ) 
            elseif   mᵢ==:off;      0.5*(λᵢ * λᵢ) 
            end
        end
        return L,noFB  #(λ=λ,g=gap,mode=mode)
    end
    #TODO Does not set FB, and DirectXUA/solve has no line search...     
    return L,noFB
end

allocate_drawing(       axis,o::AbstractVector{Teleobj};kwargs...) where{Teleobj<:ElementCostAndConstraint} = allocate_drawing(axis,[eᵢ.eleobj for eᵢ∈o];kwargs...)
function update_drawing(axis,o::AbstractVector{Teleobj},oldmut,opt, Λ,X,U,A,t,SP,dbg) where{Teleobj<:ElementCostAndConstraint{NSO,TargetElement,λxinod,λuinod}} where{NSO,TargetElement,λxinod,λuinod} 
    _,ix,_,iu,_,_ = dofpartition(eltype(o))
    Λe            = view(Λ,ix,:) 
    Xe            = getsomedofs(X,ix) 
    Ue            = getsomedofs(U,iu) 
    mut           = update_drawing(  axis,[oᵢ.eleobj for oᵢ∈o],oldmut,opt, Λe,Xe,Ue,A,t,SP,dbg)  
    return mut
end
display_drawing!(axis,::Type{<:ElementCostAndConstraint{NSO,TargetElement}},obs,opt) where{NSO,TargetElement} = display_drawing!(axis,TargetElement,obs,opt)

"""
    Muscade.QuickFix <: AbstractElement

An element for creating simple elements with "one line" of code.  
Elements thus created have several limitations:
- physical elements with only X-dofs.
- only `R` can be espied.
The element is intended for testing.  Muscade-based applications should not include this in their API. 

# Named arguments to the constructor
- `inod::NTuple{Nx,𝕫}`. The element-node numbers of the X-dofs.
- `field::NTuple{Nx,Symbol}`. The fields of the X-dofs.
- `res::Functor`, where `res(X::ℝ1,X′::ℝ1,X″::ℝ1,t::ℝ) → ℝ1`, the residual.

# Examples
A one-dimensional linear elastic spring with stiffness 2.
```jldoctest; output = false
using Muscade
model = Model(:TestModel)
node1  = addnode!(model,𝕣[0])
node2  = addnode!(model,𝕣[1])
@functor with() res(x,u,a,t)=0.4x.+.08+.5x.^2) 
e = addelement!(model,Muscade.QuickFix,[node1,node2];inod=(1,2),field=(:tx1,:tx1),res=res)

# output

Muscade.EleID(1, 1)                       
```    
"""
struct QuickFix{Nx,inod,field,Tres} <: AbstractElement
    res        :: Tres    # R = res(X,X′,X″,t)
end
QuickFix(nod::Vector{Node};inod::NTuple{Nx,𝕫},field::NTuple{Nx,Symbol},res::Functor) where{Nx} = QuickFix{Nx,inod,field,typeof(res)}(res)
doflist(::Type{<:QuickFix{Nx,inod,field}}) where{Nx,inod,field} = (inod =inod,class=ntuple(i->:X,Nx),field=(field)) 
@espy function residual(o::QuickFix, X,U,A, t,SP,dbg) 
    ☼R = o.res(∂0(X),∂1(X),∂2(X),t)
    return R,noFB
end


