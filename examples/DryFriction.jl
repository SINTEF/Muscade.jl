using Muscade,StaticArrays

# # DryFriction
#    
# Besides providing a general example of how to implement an element in `Muscade`, this element illustrates how to implement hysteretic behaviour
# (here: dry friction, but this would also apply to plasticity) without internal variables, since these make problems with `U`dofs intractable.
# This is handled by making what would otherwise have been an internal variable into an additional `X`dof.  In this case, the element's second degree
# of freedom is the friction force.
# ## Type
# The struct contains the values provided (indirectly) by the user. Note `Fx` and `Ff` which are type parameters: these will be `Symbol`s that
# represent the field of the `X`dof on which to apply the friction, and a `X`dof to represent the friction force
struct DryFriction{Fx,Ff} <: AbstractElement
    fric    :: 𝕣
    x′scale :: 𝕣  
    k⁻¹     :: 𝕣   # ∈ [0,∞[, so k ∈ ]0,∞]
end
# ## Constructor
# We provide a constructor, which will be called by `AddElement!`. The keyword arguments can, or must be given by the user when calling  `AddElement!`, and are passed on
# to the constructor. Note that the constructor is type unstable: it gets `fields` and `fieldf` as values and uses them as type parameters. This is not deemed to be a problem for
# the constructor (type instability in `residual` would be another matter)
"""
    DryFriction <: AbstractElement

Add a single-node "dry-friction" resistance to a single X-dof.  Because `Muscade`does not allow internal variables,
the element has a second dof which is the friction force.

# Named arguments to the constructor
- `fieldx::Symbol`. The field of the dof to which to apply the dry friction.
- `fieldf::Symbol = :f`. The field of the friction force dof.
- `fric::𝕣`. The absolute value of the friction force.
- `Δx::𝕣=0`. The width over which the friction force builds up.
- `x′scale::𝕣=1.`. A typical order of magnitude of the velocity of the dof to which dry friction is applied.

"""
DryFriction(nod::Vector{Node};fieldx::Symbol,fieldf::Symbol=:f,
           friction::𝕣,Δx::𝕣=0.,x′scale::𝕣=1.) =
           DryFriction{fieldx,fieldf}(friction,x′scale,Δx/friction)
# ## `residual`
# The `residual` function is prepended by `@espy` to facilitate the extraction of element-results (see [`getresults`(@ref)).
# The full name `Muscade.residual` must be used, because we are adding a method to a function defined in the `module` `Muscade`.
@espy function Muscade.residual(o::DryFriction, X,U,A, t,SP,dbg) 
    x,x′,f,f′ = ∂0(X)[1],∂1(X)[1], ∂0(X)[2], ∂1(X)[2]       
    conds     = (stick = (x′-o.k⁻¹*f′)/o.x′scale,           
                 slip  =  abs(f)/o.fric -1      )                      
    ☼old      = argmin(map(abs,conds))                      
    if        old==:stick && abs(f)>o.fric   ☼new = :slip   
    elseif    old==:slip  && f*x′<0          ☼new = :stick  
    else                                     ☼new =  old    
    end        
    return SVector(f,conds[new]), noFB
end
# In the above `f` (a force) uses the "nod-on-el" convention (force exterted by the element's node on the element), so the sign is unusual.
#
# `conds = ...`: Was the system in `stick` of `slip` at the previous iteration? Each condition is matched if the 
# expression (`(x′-o.k⁻¹*f′)/o.x′scale` or `abs(f)/o.fric -1`) evaluates to a near zero value. `conds` is a `NamedTuple`.  

# The expression `argmin(map(abs,conds))` applies `abs` to each term of the tuple, and returns the index (`:stick` or `:slip`) of the
# smallest argument: it identifies the `old` state of the element.
#
# Variables `old` and `new` are prepended by `☼` (`\sun`), to tell `@espy` this is an element-result. 
#
# The `if` construct can be read as follows:
# - If we were in stick but now `|f|` exceeds `o.fric`, we now slip.
# - If we were in slip but now the force is in the wrong direction, we now stick.
# - Otherwise, no change.
#
# The function returns a 2-vector of residuals (corresponding to the two `X`dofs), and we have no "feedback" to the solver (as opposed to constraint elements).
#
# ## `doflist`
# Another function that must be overloaded, in order to tell `Muscade` what dofs the element provides. Note that this is a function of the element *type*, not
# of the element *variable*: elements of the same concrete type must have the same dofs.
Muscade.doflist( ::Type{DryFriction{Fx,Ff}}) where{Fx,Ff} =
    (inod =(1 ,1 ), class=(:X,:X), field=(Fx,Ff)) 

