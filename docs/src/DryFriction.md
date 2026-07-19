```@meta
EditURL = "../../examples/DryFriction.jl"
```

````@example DryFriction
using Muscade,StaticArrays
````

# DryFriction

Besides providing a general example of how to implement an element in `Muscade`, this element illustrates how to implement hysteretic behaviour
(here: dry friction, but this would also apply to plasticity) without internal variables, since these make problems with `U`dofs intractable.
This is handled by making what would otherwise have been an internal variable into an additional `X`dof.  In this case, the element's second degree
of freedom is the friction force.
## Type
The struct contains the values provided (indirectly) by the user. Note `Fx` and `Ff` which are type parameters: these will be `Symbol`s that
represent the field of the `X`dof on which to apply the friction, and a `X`dof to represent the friction force

````@example DryFriction
struct DryFriction{Fx,Ff} <: AbstractElement
    fric    :: 𝕣
    x′scale :: 𝕣
    k⁻¹     :: 𝕣   # ∈ [0,∞[, so k ∈ ]0,∞]
end
````

## Constructor
We provide a constructor, which will be called by `AddElement!`. The keyword arguments can, or must be given by the user when calling  `AddElement!`, and are passed on
to the constructor. Note that the constructor is type unstable: it gets `fields` and `fieldf` as values and uses them as type parameters. This is not deemed to be a problem for
the constructor (type instability in `residual` would be another matter)

````@example DryFriction
DryFriction(nod::Vector{Node};fieldx::Symbol,fieldf::Symbol=:f,
           friction::𝕣,Δx::𝕣=0.,x′scale::𝕣=1.) =
           DryFriction{fieldx,fieldf}(friction,x′scale,Δx/friction)
````

## `residual`
The `residual` function is prepended by `@espy` to facilitate the extraction of element-results .
Variables `old` and `new` are prepended by `☼` (`\sun`), to tell `@espy` the values of these variables can
be requested using [`getresult`](@ref).

The full name `Muscade.residual` must be used, because we are adding a method to a function defined in the `module` `Muscade`.

````@example DryFriction
@espy function Muscade.residual(o::DryFriction, X,U,A, t,SP,dbg)
    x,x′,f,f′ = ∂0(X)[1],∂1(X)[1], ∂0(X)[2], ∂1(X)[2]
    stick = (x′-o.k⁻¹*f′)/o.x′scale
    slip  = abs(f)/o.fric -1
    ☼old  = abs(slip)<abs(stick) ? :slip : :stick
    if        old==:stick && abs(f)>o.fric   ☼new = :slip
    elseif    old==:slip  && f*x′<0          ☼new = :stick
    else                                     ☼new = old
    end
    return (new==:slip ? SVector(f,slip) : SVector(f,stick)), noFB
end;
nothing #hide
````

In the above `f` (a force) uses the "nod-on-el" convention (force exterted by the element's node on the element), so the sign is unusual.

If the element was in stick the previous iteration, the variable `stick = (x′-o.k⁻¹*f′)/o.x′scale` will be very close to zero.
Similarly, if the element was in slip the previous iteration, the variable `slip  = abs(f)/o.fric -1` will be very close to zero.
`old` is then either `:stick` or `:slip` depending on which of the two above is smallest.

The `if` construct can be read as follows:
- If we were in stick at previous iteration but `|f|` from previous iteration exceeds `o.fric`, we slip in this iteration.
- If we were in slip at previous iteration but force from previous iteration is in the wrong direction, we stick in this iteration.
- Otherwise, no change.

The function returns a 2-vector of residuals (corresponding to the two `X`dofs).  The first residual `f` is the friction force applied to the dof `fieldx`.

The second residual corresponds to dof `fieldf`.  Importantly, this later dof must not be shared with any other element, so the that solver will set
the second residual to zero in this iteration: depending on `new`, this will enforce `slip==0` or `stick==0`. The relevant condition will
be enforced exactly (to rounding errors) because both conditions are linear in `X`.

## `doflist`
Another function that must be overloaded, in order to tell `Muscade` what dofs the element provides. Note that this is a function of the element *type*, not
of the element *variable*: elements of the same concrete type must have the same dofs.

````@example DryFriction
Muscade.doflist( ::Type{DryFriction{Fx,Ff}}) where{Fx,Ff} =
    (inod =(1 ,1 ), class=(:X,:X), field=(Fx,Ff))
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

