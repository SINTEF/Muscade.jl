struct DofCost{Class,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,Tcost,Tcostargs} <: AbstractElement
    cost     :: Tcost    # Class==:instant cost(X,U,A,t,costargs...), Class==:A cost(A,costargs...) X and U are tuples (derivates of dofs...) 
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


@espy function lagrangian(o::DofCost{:I,Nx,Nu,Na},δX,X,U,A,t,χ,χcv,SP,dbg) where{Nx,Nu,Na} 
    ☼cost = o.cost(X,U,A,t,o.costargs...)
    return cost,noχ,noFB
end


@espy function lagrangian(o::DofCost{:A,Nx,Nu,Na},δX,X,U,A,t,χ,χcv,SP,dbg) where{Nx,Nu,Na} 
    ☼cost = o.cost(    A  ,o.costargs...)
    return cost,noχ,noFB
end

struct ElementCost{Teleobj,Treq,Tcost,Tcostargs}
    eleobj   :: Teleobj
    req      :: Treq
    cost     :: Tcost     
    costargs :: Tcostargs
end
function ElementCost(nod::Vector{Node};req,cost,costargs=(;),ElementType,elementkwargs...)
    eleobj   = ElementType(nod;elementkwargs...)
    return ElementCost(eleobj,req,cost,costargs)
end
doflist( ::Type{<:ElementCost{Teleobj}}) where{Teleobj} = doflist(Teleobj)
@espy function lagrangian(o::ElementCost, δX,X,U,A,t,χ,χcv,SP,dbg)
    @show typeof(o.eleobj)
    L,χ,FB,eleres  = ☼lagrangian(o.eleobj,δX,X,U,A,t,χ,χcv,SP,(dbg...,via=ElementCost),o.req)
    ☼cost          = o.cost(eleres,X,U,A,t,o.costargs...) 
    return L+cost,χ,FB
end    

#-------------------------------------------------
"""
`DofCost{Derivative,Class,Field,Tcost} <: AbstractElement`

An element with a single node, for adding a cost to a given dof.  

# Named arguments to the constructor
- `class::Symbol`, either `:X`, `:U` or `:A`.
- `field::Symbol`.
- `cost::Function`, where `cost(x::ℝ,t::ℝ[,costargs...]) → ℝ`.

# Requestable internal variables
- `cost`, the value of the cost.

# Examples
```jldoctest; output = false
using Muscade
model = Model(:TestModel)
node  = addnode!(model,𝕣[0,0])
e     = addelement!(model,DofCost,[node];class=:X,field=:tx,costargs=(3.,),cost=(x,t,three)->(x/three)^2)

# output

EleID(1, 1)
```    
See also: [`Hold`](@ref), [`DofLoad`](@ref)
"""
struct SingleDofCost <: AbstractElement end
function SingleDofCost(nod::Vector{Node};class::Symbol,field::Symbol,cost::Function,derivative=0::𝕫,costargs=()) 
    if     class==:X; DofCost(nod;xinod=(1,),xfield=(field,),class=:I,cost=(X,U,A,t,args...)->cost(∂n(X,derivative)[1],t,args...),costargs)
    elseif class==:U; DofCost(nod;uinod=(1,),ufield=(field,),class=:I,cost=(X,U,A,t,args...)->cost(∂n(U,derivative)[1],t,args...),costargs)
    elseif class==:A; DofCost(nod;ainod=(1,),afield=(field,),class=:A,cost=(    A,  args...)->cost(A[1]                 ,args...),costargs)
    else              muscadeerror("'class' must be :X,:U or :A")
    end
end    

#-------------------------------------------------

"""
`DofLoad{Tvalue,Field} <: AbstractElement`

An element to apply a loading term to a single X-dof.  

# Named arguments to the constructor
- `field::Symbol`.
- `value::Function`, where `value(t::ℝ) → ℝ`.

# Requestable internal variables
- `F`, the value of the load.

# Examples
```jldoctest; output = false
using Muscade
model = Model(:TestModel)
node  = addnode!(model,𝕣[0,0])
e     = addelement!(model,DofLoad,[node];field=:tx,value=t->3t-1)

# output

EleID(1, 1)
```    

See also: [`Hold`](@ref), [`DofCost`](@ref)  
"""
struct DofLoad{Tvalue,Field} <: AbstractElement
    value      :: Tvalue # Function
end
DofLoad(nod::Vector{Node};field::Symbol,value::Tvalue) where{Tvalue<:Function} = DofLoad{Tvalue,field}(value)
doflist(::Type{DofLoad{Tvalue,Field}}) where{Tvalue,Field}=(inod=(1,), class=(:X,), field=(Field,))
@espy function residual(o::DofLoad, X,U,A,t,χ,χcv,SP,dbg) 
    ☼F = o.value(t)
    return SVector{1}(-F),noχ,noFB
end

#-------------------------------------------------

#McCormick(a,b)= α->a*exp(-(α/b)^2)            # provided as input to solvers, used by their Addin
decided(λ,g,γ)  = abs(VALUE(λ)-VALUE(g))/γ    # used by constraint elements

S(λ,g,γ) = (g+λ-hypot(g-λ,2γ))/2 # Modified interior point method's take on KKT's-complementary slackness 

KKT(λ::𝕣        ,g::𝕣         ,γ::𝕣,λₛ,gₛ)                 = 0 # A pseudo-potential with strange derivatives
KKT(λ::∂ℝ{P,N,R},g::∂ℝ{P,N,R},γ::𝕣,λₛ,gₛ) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(0, λ.x*g.dx + gₛ*S(λ.x/λₛ,g.x/gₛ,γ)*λ.dx)
KKT(λ:: ℝ       ,g::∂ℝ{P,N,R},γ::𝕣,λₛ,gₛ) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(0, λ.x*g.dx                           )
KKT(λ::∂ℝ{P,N,R},g:: ℝ       ,γ::𝕣,λₛ,gₛ) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(0,            gₛ*S(λ.x/λₛ,g.x/gₛ,γ)*λ.dx)
function KKT(λ::∂ℝ{Pλ,Nλ,Rλ},g::∂ℝ{Pg,Ng,Rg},γ::𝕣,λₛ,gₛ) where{Pλ,Pg,Nλ,Ng,Rλ<:ℝ,Rg<:ℝ}
    if Pλ==Pg
        R = promote_type(Rλ,Rg)
        return ∂ℝ{Pλ,Nλ}(convert(R,KKT(λ.x,g.x,γ,λₛ,gₛ)),convert.(R,     λ.x*g.dx + gₛ*S(λ.x/λₛ,g.x/gₛ,γ)*λ.dx))
    elseif Pλ> Pg
        R = promote_type(Rλ,typeof(b))
        return ∂ℝ{Pλ,Nλ}(convert(R,KKT(λ  ,g.x,γ,λₛ,gₛ)),convert.(R,     λ.x*g.dx                            ))
    else
        R = promote_type(typeof(a),Rg)
        return ∂ℝ{Pg,Ng}(convert(R,KKT(λ.x,g  ,γ,λₛ,gₛ)),convert.(R,                gₛ*S(λ.x/λₛ,g.x/gₛ,γ)*λ.dx))
    end
end

#-------------------------------------------------

"""
`off(t) → :off`

See also: [`Constraint`](@ref), [`equal`](@ref), [`inequal`](@ref)
"""
off(t)     = :off
"""
`equal(t) → :equal`

See also: [`Constraint`](@ref), [`off`](@ref), [`inequal`](@ref)
"""
equal(t)   = :equal
"""
`inequal(t) → :inequal`

See also: [`Constraint`](@ref), [`off`](@ref), [`equal`](@ref)
"""
inequal(t) = :inequal
# length of comment                                           stop here|
"""
`Constraint{λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield,Tg,Tmode} <: AbstractElement`

An element to apply physical/optimisation equality/inequality constraints on dofs. 

The constraints are holonomic, i.e. they apply to the values, not the time derivatives, of the involved dofs. 
This element is very general but not very user-friendly to construct, factory functions are provided for better useability. 
The sign convention is that the gap `g≥0` and the Lagrange multiplier `λ≥0`.

# Named arguments to the constructor
- `xinod::NTuple{Nx,𝕫}=()` For each X-dof to be constrained, its element-node number.
- `xfield::NTuple{Nx,Symbol}=()` For each X-dof to be constrained, its field.
- `uinod::NTuple{Nu,𝕫}=()` For each U-dof to be constrained, its element-node number.
- `ufield::NTuple{Nu,Symbol}=()` For each U-dof to be constrained, its field.
- `ainod::NTuple{Na,𝕫}=()` For each A-dof to be constrained, its element-node number.
- `afield::NTuple{Na,Symbol}=()` For each A-dof to be constrained, its field.
- `λinod::𝕫` The element-node number of the Lagrange multiplier.
- `λclass::Symbol` The class of the Lagrange multiplier. `:X` for physical constraints, `:U` for optimisation constraints. `:A` is experimental.
- `λfield::Symbol` The field of the Lagrange multiplier.
- `gₛ::𝕣=1.` A scale for the gap.
- `λₛ::𝕣=1.` A scale for the Lagrange multiplier.
- `g::Function` For physical constraints: `g(X::ℝ1,t::ℝ) -> ℝ`, for physical constraints and `g(X::ℝ1,U::ℝ1,A::ℝ1,t::ℝ) -> ℝ`, for optimisation constraints.
- `mode::Function`, where `mode(t::ℝ) -> Symbol`, with value `:equal`, `:inequal` or `:off` at any time. An `:off` constraint will set the Lagrange multiplier to zero.

# Examples
```jldoctest
using Muscade
model           = Model(:TestModel)
n1              = addnode!(model,𝕣[0]) 
e1              = addelement!(model,Constraint,[n1],xinod=(1,),xfield=(:t1,),
                              λinod=1, λclass=:X, λfield=:λ1,g=(x,t)->x[1]+.1,mode=inequal)
e2              = addelement!(model,QuickFix  ,[n1],inod=(1,),field=(:t1,),res=(x,u,a,t)->0.4x.+.08+.5x.^2)
state           = solve(staticX;model,time=[0.],verbose=false) 
X               = state[1].X[1]

# output

2-element Vector{Float64}:
 -0.09999867546403915
  0.045000397353771225
```    

See also: [`Hold`,`off`,`equal`,`inequal`](@ref)
"""
struct Constraint{λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield,Tg,Tgargs,Tmode} <: AbstractElement
    g        :: Tg    # Class==:X g(x,t,gargs...) ,Class==:U  g(x,u,a,t,gargs...), Class==:A g(a,gargs...) 
    gargs    :: Tgargs
    mode     :: Tmode # mode(t)->symbol, or Symbol for Aconstraints
    gₛ        :: 𝕣
    λₛ        :: 𝕣  
end
function Constraint(nod::Vector{Node};xinod::NTuple{Nx,𝕫}=(),xfield::NTuple{Nx,Symbol}=(),
                                      uinod::NTuple{Nu,𝕫}=(),ufield::NTuple{Nu,Symbol}=(),
                                      ainod::NTuple{Na,𝕫}=(),afield::NTuple{Na,Symbol}=(),
                                      λinod::𝕫, λclass::Symbol, λfield::Symbol,
                                      gₛ::𝕣=1.,λₛ::𝕣=1.,
                                      g::Function ,gargs=(),mode::Function) where{Nx,Nu,Na} 
    (λclass==:X && (Nu>0||Na>0)) && muscadeerror("Constraints with λclass=:X must have zero U-dofs and zero A-dofs") 
    (λclass==:A && (Nu>0||Na>0)) && muscadeerror("Constraints with λclass=:A must have zero X-dofs and zero U-dofs") 
    return Constraint{λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield,
                       typeof(g),typeof(gargs),typeof(mode)}(g,gargs,mode,gₛ,λₛ)
end
doflist(::Type{<:Constraint{λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield}}) where
                            {λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield} = 
   (inod =(xinod...           ,uinod...           ,ainod...           ,λinod ), 
    class=(ntuple(i->:X,Nx)...,ntuple(i->:U,Nu)...,ntuple(i->:A,Na)...,λclass    ), 
    field=(xfield...          ,ufield...          ,afield...          ,λfield)) 
@espy function residual(o::Constraint{:X,Nx}, X,U,A,t,χ,χcv,SP,dbg) where{Nx}
    γ          = default{:γ}(SP,0.)
    P,gₛ,λₛ     = constants(∂0(X)),o.gₛ,o.λₛ
    x,☼λ       = ∂0(X)[SVector{Nx}(1:Nx)], ∂0(X)[Nx+1]   
    x∂         = variate{P,Nx}(x) 
    ☼g,g∂x     = value_∂{P,Nx}(o.g(x∂,t,o.gargs...)) 
    return if  o.mode(t)==:equal;   SVector{Nx+1}((       -g∂x*λ)...,-g              ) ,noχ,(α=∞                  ,)
    elseif     o.mode(t)==:inequal; SVector{Nx+1}((       -g∂x*λ)...,-gₛ*S(λ/λₛ,g/gₛ,γ)) ,noχ,(α=decided(λ/λₛ,g/gₛ,γ),)
    elseif     o.mode(t)==:off;     SVector{Nx+1}(ntuple(i->0,Nx)...,-gₛ/λₛ*λ         ) ,noχ,(α=∞                  ,)
    end
end
@espy function lagrangian(o::Constraint{:U,Nx,Nu,Na}, δX,X,U,A,t,χ,χcv,SP,dbg) where{Nx,Nu,Na}
    γ          = default{:γ}(SP,0.)
    x,u,a,☼λ   = ∂0(X),∂0(U)[SVector{Nu}(1:Nu)],A,∂0(U)[Nu+1]
    ☼g         = o.g(x,u,a,t,o.gargs...)
    return if  o.mode(t)==:equal;   -g*λ                  ,noχ,(α=∞                      ,)
    elseif     o.mode(t)==:inequal; -KKT(λ,g,γ,o.λₛ,o.gₛ)  ,noχ,(α=decided(λ/o.λₛ,g/o.gₛ,γ),)
    elseif     o.mode(t)==:off;     -o.gₛ/(2o.λₛ)*λ^2      ,noχ,(α=∞                      ,)  
    end
end
@espy function lagrangian(o::Constraint{:A,Nx,Nu,Na}, δX,X,U,A,t,χ,χcv,SP,dbg) where{Nx,Nu,Na}
    γ          = default{:γ}(SP,0.)
    a,☼λ       = A[SVector{Na}(1:Na)],A[    Na+1] 
    ☼g         = o.g(a,o.gargs...)
    return if  o.mode(t)==:equal;   -g*λ                  ,noχ,(α=∞                      ,) 
    elseif     o.mode(t)==:inequal; -KKT(λ,g,γ,o.λₛ,o.gₛ)  ,noχ,(α=decided(λ/o.λₛ,g/o.gₛ,γ),)
    elseif     o.mode(t)==:off;     -o.gₛ/(2o.λₛ)*λ^2      ,noχ,(α=∞                      ,)   
    end
end

#-------------------------------------------------

"""
`Hold <: AbstractElement`

An element to set a single X-dof to zero.  

# Named arguments to the constructor
- `field::Symbol`. The field of the X-dof to constraint.
- `λfield::Symbol=Symbol(:λ,field)`. The field of the Lagrange multiplier.

# Examples
```jldoctest; output = false
using Muscade
model = Model(:TestModel)
node  = addnode!(model,𝕣[0,0])
e     = addelement!(model,Hold,[node];field=:tx)

# output

EleID(1, 1)
```    

See also: [`Constraint`](@ref), [`DofLoad`](@ref), [`DofCost`](@ref) 
"""
struct Hold <: AbstractElement end  
function Hold(nod::Vector{Node};field::Symbol,λfield::Symbol=Symbol(:λ,field)) 
    g(v,t)=v[1]
    return Constraint{:X     ,1, 0, 0, (1,),(field,),(),   (),    (),   (),    1,    λfield, typeof(g),typeof(()),typeof(equal)}(g,(),equal,1.,1.)
    #      Xconstraint{λclass,Nx,Nu,Na,xinod,xfield, uinod,ufield,ainod,afield,λinod,λfield}
end

#-------------------------------------------------

"""
`QuickFix <: AbstractElement`

An element for creating simple elements with "one line" of code.  
Elements thus created have several limitations:
- no internal state.
- no initialisation.
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

EleID(1, 1)                       
```    
"""
struct QuickFix{Nx,inod,field,Tres} <: AbstractElement
    res        :: Tres    # R = res(X,X′,X″,t)
end
QuickFix(nod::Vector{Node};inod::NTuple{Nx,𝕫},field::NTuple{Nx,Symbol},res::Function) where{Nx} = QuickFix{Nx,inod,field,typeof(res)}(res)
doflist(::Type{<:QuickFix{Nx,inod,field}}) where{Nx,inod,field} = (inod =inod,class=ntuple(i->:X,Nx),field=(field)) 
@espy function residual(o::QuickFix, X,U,A, t,χ,χcv,SP,dbg) 
    ☼R = o.res(∂0(X),∂1(X),∂2(X),t)
    return R,noχ,noFB
end

#-------------------------------------------------
