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
struct DryFriction{Fx,Ff} <: AbstractElement
    fric    :: 𝕣
    x′scale :: 𝕣  
    k⁻¹     :: 𝕣   # ∈ [0,∞[, so k ∈ ]0,∞]
end
DryFriction(nod::Vector{Node};fieldx::Symbol,fieldf::Symbol=:f,friction::𝕣,Δx::𝕣=0.,x′scale::𝕣=1.) = DryFriction{fieldx,fieldf}(friction,x′scale,Δx/friction)
@espy function Muscade.residual(o::DryFriction, X,U,A, t,SP,dbg) 
    x,x′,f,f′ = ∂0(X)[1],∂1(X)[1], ∂0(X)[2], ∂1(X)[2]       # f: nod-on-el convention, the sign is unusual.
    conds     = (stick = (x′-o.k⁻¹*f′)/o.x′scale,           # Was the system in stick of slip at the previous iteration?
                 slip  =  abs(f)/o.fric -1      )           #  - each condition is matched if expression evals to 0.       
    ☼old      = argmin(map(abs,conds))                      # Symbol-index of the "most matched" condition
    if        old==:stick && abs(f)>o.fric   ☼new = :slip   # if we were in stick but now |f| exceeds o.fric, we now slip
    elseif    old==:slip  && f*x′<0          ☼new = :stick  # if we were in slip but now the force is is the wrong direction, we now stick
    else                                     ☼new =  old    # otherwise, no change
    end                  
    return SVector(f,conds[new]), noFB
end
Muscade.doflist( ::Type{DryFriction{Fx,Ff}}) where{Fx,Ff} = (inod =(1 ,1 ), class=(:X,:X), field=(Fx,Ff)) 

