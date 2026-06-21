
# [Automatic differentiation](@id adiff_appendix)

## Theory
This section aims to provide an understanding of a variety of usefull techniques in automatic differentiation.

`Muscade` has its own implementation of forward automatic differentiation for historical reasons: Prototypes automatic differentiation of `Muscade` where developed in parallel with [`ForwardDiff.jl`](https://juliadiff.org/ForwardDiff.jl/stable/).  The inner workings of [`ForwardDiff.jl`](https://juliadiff.org/ForwardDiff.jl/stable/) and `Muscade`'s version are very similar. 

The APIs, on the other hand, are quite different: `ForwardDiff.jl` is designed for safety, `Muscade`'s automatic differentiation for flexibility. When using the techniques described in the following, it is possible to introduce some confusing bugs if one is not quite careful.

 `Muscade` also offers some innovative functionality (see the Section on [chain rule](@ref chain_rule)).

### Type definition
`∂ℝ` is defined as

```julia
struct ∂ℝ{P,N,R} <:ℝ where{R<:ℝ}  
    x  :: R
    dx :: SVector{N,R}
end
```

where `N` is the number of partial derivatives, `P` is the *precedence* (to prevent [perturbation confusion](https://arxiv.org/pdf/1211.4892v3), discussed in the following), and `R` is the underlying datatype (either `Float64`, noted `𝕣` in Muscade, or another `∂ℝ` when computing higher order derivatives).

### Techniques and tools

**Gradient evaluation**

```julia
using StaticArrays
using Muscade: variate,value,∂
f(x) = x.*x .+ sum(x)
x    = SVector(1.,2.,3.)
P,N  = 1,length(x)
x1   = variate{P,N}(x)
y1   = f(x1)
y    = value{P}(y1)
yₓ   = ∂{P,N}(y1)    
```

[`variate`](@ref) builds a `SVector` of `∂{P,N,R}`, where `R == eltype(x)`. `x` *must* be a `SVector{N}`. Once the calculations to be differentiated are completed, results (here:`y`) are `StaticArrays` of `∂{P,N,R}`.  The partials are "seeded", so that `∂{P,N}(x1)[i]` is a `SVector{N}` of `zero(R)`, except for the `i`-th term that is `one(R)`. 

[`value`](@ref)`{P}`, where `y1` is a `StaticArray{S,∂{P,N,R}}`, retrieves the value `y`, of type `StaticArray{S,R}`

[`∂`](@ref)`{P,N}`, where `y1` is a `StaticArray{S,∂{P,N,R}}`, retrieves the gradient `yₓ`, of type `StaticArray{Syₓ,R}`: `yₓ` has one more dimension than `y`, of size `N`, corresponding to the pratial derivatives with respect to the `SVector` `x`.

**Joint derivatives**
One may wish to compute partial derivatives with respect to various variables at the same time. An example of the technique is

```julia
using StaticArrays
using Muscade: variate,value,∂
res(x,u) = x.*x .+ sum(u) # SVector-valued 
x        = SVector(1.,2.)
u        = SVector(3.,4.,5.)
P,N,R    = 1,5,𝕣
ix       = SVector(1,2)
iu       = SVector(3,4,5)
δs       = δ{P,N,R}()
x1       = x + δs[ix]
u1       = u + δs[iu]
r1       = res(x1,u1)
r        = value{P}(r1)
rₓ       = ∂{P,N}(r1)[:,ix]    
rᵤ       = ∂{P,N}(r1)[:,iu]    
```

`δ` is similar to `variate`, except that all the `value`s of `δs` are `zero(R)`.  In the above the elements of `x1`, `u2` and `r1` all have 5 partials - the first two corresponding to `x` and the rest to `u`. 

The alove code is quite clumsy, and the function [`revariate`](@ref) allows to do the same more elegantly:

XXXXXX revariate for simpler code XXXX  

**Cross derivatives** provide the best way to understand the mesting of `∂ℝ` variables can be nested.

```julia
using StaticArrays
using Muscade: variate,value,∂
res(x,u) = SVector(x[1]*x[1]+u[2],x[2]*u[3])
x        = SVector(1.,2.)
u        = SVector(3.,4.,5.)
Px,Nx    = 1,2
Pu,Nu    = 2,3 # Pu > Px
x1       = variate{Px,Nx}(x)
u2       = variate{Pu,Nu}(u)
r2       = res(x1,u2)
r        = value{Px}(value{Pu}(r2))
rₓ       = ∂{Px,Nx }(value{Pu}(r2))    
rᵤ       = value{Px}(∂{Pu,Nu }(r2))    
rᵤₓ      = ∂{Px,Nx }(∂{Pu,Nu }(r2))  
```

In the above, `typeof(x1) == ∂ℝ{1,2,𝕣}`, and `typeof(u2) == ∂ℝ{2,3,𝕣}`.  When these two variables "touch" in a binary operation in `res`, this results 
in `typeof(r2) == ∂ℝ{2,3,∂ℝ{1,2,𝕣}}`. 

Note the specific order in which `r2` is unpacked to obtain the derivatives: first precedence `2` is unpacked (because `Muscade.precedence(r2) == 2`), then precedence `1`.  Unpacking `r2` directly, with a precedence lower than `2` causes an error, as it amounts to "peeling the onion from the inside".

Unpacking `r2` with any precedence higher than `2` will cause `r2` to be treated as a constant: `value{3}(r2) == r2`, and `∂{3,6}(r2)` returns an array of zeros. 

**Higher order derivatives**

**Nested derivatives**

**Directional derivatives**

**[Chain rule](@id chain_rule)**

## Tools for differentiation
motion
revariate, McLaurin, apply

-------------------------------------

**Automatic differentiation** generaly does not affect how `Muscade.residual` and `Muscade.lagrangian` are written.  There are two performance-related exceptions to this:

1. If a complicated sub-function in `Muscade.residual` and `Muscade.lagrangian` (typicaly a material model or other closure) operates on an array (for example, the strain) that is smaller than the number of degrees of freedom of the system, computing time can be saved by computing the derivative of the output (in the example, the stress) with respect to the input to the subfunction, and then chainrule the derivatives.
2. Iterative precedures are sometimes used within `Muscade.residual` and `Muscade.lagrangian`, a typical example being in plastic material formulations.  There is no need to propagate automatic differentiation through all the iterations - doing so with the result of the iteration provides the same result.
3. Elements with corotated reference system (e.g. [beam elements](StaticBeamAnalysis.md)) can use automatic differentiation to transform the residual back to the global reference system.

See the page on [automatic differentiation](Adiff.md).

## Automatic differentiation within element code

Some advanced elements (in particular, elements with co-rotated element systems) can be implemented elegantly by using automatic differentiation within `residual` or `lagrangian`.  These are advanced techniques, requiring a good understanding of [`automatic differentiation`](Adiff.md).  Example of usage can be found in [`toolbox/BeamElement.jl`](StaticBeamAnalysis.md).

Helper functions [`motion`](@ref) and [`motion⁻¹`](@ref) allow to transform a `tuple` of `SVectors`, like the input `X` given to `residual` and `lagrangian`, into a an automatic differentiation structure, so that functions of `∂0(X)` only can be differentiated with respect to time. 

It is sometimes possible to improve performance by identifying a part of `residual` or `lagrangian` which takes a single, `SVector` as an input: A vector shorter than the list of dofs differentiated by the solver allow to accelerate computations, by using [`apply`](@ref), or for more advanced usage, [`revariate`](@ref) in combination with [`chainrule`](@ref). 

In [`toolbox/BeamElement.jl`](StaticBeamAnalysis.md), in function `kinematics`, [`apply`](@ref) is applied to accelerate a process of differentiation to the 2nd order.  In `residual`, [`revariate`](@ref) and [`chainrule`](@ref) in order to differentiate `kinematics` and accelerate computations by exploiting the fact that `kinematic` is a function of `∂0(X)` only.




## Usage

`Muscade`s automatic differentiation is used as follows:





 Where automatic differentiation is nested, the extractions must be carried out in reverse order:

```julia
x   = SVector(1.,2.,3.)
N   = length(x)
x1  = variate{1,N}(x)
x2  = variate{2,N}(x1)
y2  = f(x2)
y   = value{1}(value{2}(y2))
yₓ  = ∂{1,N}(value{2}(y2))  
yₓ  = value{1}(∂{2,N}(y2))
yₓₓ = ∂{1,N}(∂{2,N}(y2))
```

The first type parameter given to [`variate`](@ref), [`value`](@ref) and [`∂`](@ref) is the *precedence* `P`.  When automatic differentiation is nested, `variate` must be called with a value of `P` that is higher than the precedence of an variable that will influence the output. Inside a function, this precedence can vary depending on the type of arguments provided to the function. To this end, the function [`constants`](@ref) is provided:

```julia
using SVector
...
N   = length(x)
P   = constants(x,a,b,c)
x1  = variate{P,N}(x)
y1  = f(x1,a,b,c)
y   = value{P}(y1)
yₓ  = ∂{P,N}(y2)
```

## Taylor expansions

One way to accelerate automatic differentiation of complicated functions can be to chainrule the differentiation of simpler functions, in particular function with smaller inputs.  Use [`apply`](@ref) for basic applications, and [`revariate`](@ref) and [`chainrule`](@ref) for advanced usage. 
