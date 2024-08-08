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
# # DryFriction
# See [`Muscade.DryFriction`](@ref) for reference manual.
#    
# The struct contains the values provided (indirectly) by the user. Note `Fx` and `Ff` which are type parameters: these will be `Symbol`s that
# represent the field of the Xdof on which to apply the friction, and a Xdof to represent the friction force
struct DryFriction{Fx,Ff} <: AbstractElement
    fric    :: 𝕣
    x′scale :: 𝕣  
    k⁻¹     :: 𝕣   # ∈ [0,∞[, so k ∈ ]0,∞]
end
# We provide a constructor, which will be clled by `AddElement!`. The keyword arguments can, or must be given by the user when calling  `AddElement!`, and are passed on
# to the constructor. Note that the constructor is type unstable: it gets `fields` and `fieldf` as values and uses them as type parameters. This is not deemed to be a problem for
# the constructor (type instability in `residual` would be another matter)
DryFriction(nod::Vector{Node};fieldx::Symbol,fieldf::Symbol=:f,friction::𝕣,Δx::𝕣=0.,x′scale::𝕣=1.) = DryFriction{fieldx,fieldf}(friction,x′scale,Δx/friction)
# The `residual` function is prepended by `@espy` to facilitate the extraction of element-results (see [`getresults`(@ref)).
# The full name `Muscade.residual` must be used, because we are adding a method to a function defined in the `module` `Muscade`.
@espy function Muscade.residual(o::DryFriction, X,U,A, t,SP,dbg) 
    x,x′,f,f′ = ∂0(X)[1],∂1(X)[1], ∂0(X)[2], ∂1(X)[2]       # f: nod-on-el convention, the sign is unusual.
    conds     = (stick = (x′-o.k⁻¹*f′)/o.x′scale,           # Was the system in stick of slip at the previous iteration?
                 slip  =  abs(f)/o.fric -1      )           #  - each condition is matched if expression evals to 0.   
# Variables `old` and `new` are prepended by `☼` (`\sun`), to tell `@espy` this is an element-result.            
    ☼old      = argmin(map(abs,conds))                      # Symbol-index of the "most matched" condition
    if        old==:stick && abs(f)>o.fric   ☼new = :slip   # if we were in stick but now |f| exceeds o.fric, we now slip
    elseif    old==:slip  && f*x′<0          ☼new = :stick  # if we were in slip but now the force is is the wrong direction, we now stick
    else                                     ☼new =  old    # otherwise, no change
    end        
# We return a 2-vector or residuals (corresponding to the two Xdofs), and we have no "feedback" to the solver (as opposed to constraint elements).
    return SVector(f,conds[new]), noFB
end
# Another function that must be overloaded, in order to tell `Muscade` what dofs the element provides. Note that this is a function of the element *type*, not
# of the element *variable*: elements of the same concrete type must have the same dofs.
Muscade.doflist( ::Type{DryFriction{Fx,Ff}}) where{Fx,Ff} = (inod =(1 ,1 ), class=(:X,:X), field=(Fx,Ff)) 

