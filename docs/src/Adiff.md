
# [Automatic differentiation](@id adiff_appendix)

## Introduction
This section aims to provide an understanding of a variety of useful techniques in automatic differentiation. In `Muscade`, element-developers generaly need not concern themselves with automatic differentiation, except on two points:

First, do not specify the type of inputs `Λ`, `X`, `U`, `A` and `t` to `Muscade.residual` and `Muscade.lagrangian`.  For example

```julia
@espy function Muscade.lagrangian(o::MyElement,Λ::Vector{𝕣},X,U,A,t,SP,dbg)
#                                               |bad_idea|
```

(where `𝕣` is a `Muscade` alias for `Float64`) would cause a `MethodError`, because `Muscade` will attempt to call `lagrangian` with a `SVector{L,∂ℝ...}`, and the method only accepts a `Vector{𝕣}`.

Second, beware of branches (`if` and `cond ? opt1 : opt2`).  In the following code

```julia
if x>0
    y=x  # y is of type ∂ℝ
else
    y=0. # y is of type 𝕣
end    
```

the type of `y` depends on the *value* of `x` and thus cannot be determined at compilation.  This is [type instability](TypeStable.md) and results in *very* slow execution.

However, for advanced element implementation, `residual` or `lagrangian` may themselves use automatic differentiation, and the following provides the necesary guidance.

!!! warning
    Read and understand the section on "perturbation confusion" before using any of the techniques described below.

`Muscade` has its own implementation of forward automatic differentiation for historical reasons: `Muscade`'s developers started dabling in automatic differentiation while [`ForwardDiff.jl`](https://juliadiff.org/ForwardDiff.jl/stable/) was still very much in the works.  

The inner workings of [`ForwardDiff.jl`](https://juliadiff.org/ForwardDiff.jl/stable/) and `Muscade`'s version are very similar. The APIs, on the other hand, are quite different: [`ForwardDiff.jl`](https://juliadiff.org/ForwardDiff.jl/stable/) is designed for safety, `Muscade`'s automatic differentiation for flexibility. When using the techniques described in the following, it is possible to introduce some confusing bugs if one is not quite careful. `Muscade` also offers some innovative differentiation functionality (see the Section on chain rule).

## Type definition

`Muscade`'s data-type for automatic differentiation, `∂ℝ`, is defined as

```julia
struct ∂ℝ{P,N,R} <:ℝ where{R<:ℝ}  
    x  :: R             # value
    dx :: SVector{N,R}  # partials
end
```

where `ℝ` is a `Muscade` alias for `Real`.  The type parameters of `∂ℝ` is

- `P` is the *precedence*, or order of differentiation (more on this later).
- `N` is the number of partial derivatives.
- `R` is the underlying datatype (either `𝕣`, or another `∂ℝ` when computing higher order derivatives).

A variable `v` of precedence 1 is printed as `∂₁⟨v.x|v.dx⟩`, for example `∂₁⟨1.23|0,0.6,1.2⟩` (here with three partials).

This structure is nested to represent higher derivatives, so `R` can itself be a `∂ℝ`. The rule is that in any `∂ℝ{P,N,R}`, we must have `P==1+predence(R)`, where the precedence of `𝕣`  is 0. A variable `v` of precedence 2 is printed as `∂₂⟨∂₁⟨v.x.x|v.x.dx⟩|∂₁⟨x.dxᵢ.x|x.dxᵢ.dx⟩ ∀i⟩`. For example, 

```julia
v = SVector(2.,3.);
P,N,x2=variate{2}(v);
v2
```

(`variate` will be discussed in the following) gives

```julia
2-element SVector{2, Muscade.∂ℝ{2, 2, Muscade.∂ℝ{1, 2, Float64}}} with indices SOneTo(2):
 ∂₂⟨∂₁⟨2|1,0⟩|∂₁⟨1|0,0⟩,∂₁⟨0|0,0⟩⟩
 ∂₂⟨∂₁⟨3|0,1⟩|∂₁⟨0|0,0⟩,∂₁⟨1|0,0⟩⟩
```

In the following, the informal notations `∂₂⟨∂₁⟨v.x.x|v.x.dx⟩|∂₁⟨v.dx.x|v.dx.dx⟩⟩` or `∂₂⟨∂₁⟨v|vₓ⟩|∂₁⟨vₓ|vₓₓ⟩⟩` may be used.  Both notations gloss over the fact that `dx` is a vector, but they provide a useful mental picture for when unpacking such variables. These notations are not executable: they can not be used to create a `∂ℝ` variable.

## Gradient evaluation

Gradient evaluation is the simplest form for automatic differentiation.

```julia
using StaticArrays
using Muscade: variate,value,∂

f(x)   = x.*x .+ sum(x)
x      = SVector(1.,2.,3.)
P,N,x1 = variate(x) # P=1,N=length(x)
y1     = f(x1)
y      = value{P}(y1)
yₓ     = ∂{P,N}(y1)    
```

[`Muscade.variate`](@ref) here builds a `SVector` of `∂ℝ{P,N,R}`, where `R == eltype(x) == 𝕣`, `P == Muscade.precedence(x)+1 == 1` (the precedence of an `𝕣` is 0) and `N == 3` is the length of `x`.
. `Muscade` developers speak of `x1` being the "variated" form of `x`". The partials are "seeded", so that `∂{P,N}(x1)[i]` is a `SVector{N}` of `zero(R)`, except for the `i`-th term that is `one(R)`.  

```julia
julia> x1
3-element SVector{3, Muscade.∂ℝ{1, 3, Float64}} with indices SOneTo(3):
 ∂₁⟨1|1,0,0⟩ 
 ∂₁⟨2|0,1,0⟩ 
 ∂₁⟨3|0,0,1⟩ 
```

Here we can note the identify matrix that appears: the Jacobian of the vector `x` with respect to itself.

```julia
julia> y1
3-element SVector{3, Muscade.∂ℝ{1, 3, Float64}} with indices SOneTo(3):
 ∂₁⟨ 7|3,1,1⟩ 
 ∂₁⟨10|1,5,1⟩ 
 ∂₁⟨15|1,1,7⟩ 
``` 

[`value`](@ref)`{P}(y1)`, where `y1` is a `StaticArray{S,∂{P,N,R}}`, retrieves the value `y`, of type `StaticArray{S,R}`. 

[`∂`](@ref)`{P,N}(y1)`, where `y1` is a `StaticArray{S,∂{P,N,R}}`, retrieves the gradient `yₓ`, of type `StaticArray{Syₓ,R}`. `size(yₓ) = (size(y)...,N)`
where the last index refers to the partial derivatives with respect to the `SVector` `x`.

```julia
julia> y
3-element SVector{3, Float64} with indices SOneTo(3):
  7.0
 10.0
 15.0

julia> yₓ
3×3 SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 3.0  1.0  1.0
 1.0  5.0  1.0
 1.0  1.0  7.0 
```

## Joint derivatives

Joint derivatives are partial derivatives with respect to several variables, computed at the same time. An example of the technique is

```julia
using StaticArrays
using Muscade: variate,variate_indices,value,∂,δ

res(x,u)    = x.*x .+ sum(u)         # SVector-valued 
x           = SVector(1.,2.)
u           = SVector(3.,4.,5.)
P,N,(x1,u1) = variate((x,u))         # note the parentheses
(ix,iu)     = variate_indices((x,u))
r1          = res(x1,u1)
r           = value{P}(r1)
rₓ          = ∂{P,N}(r1)[:,ix]       # the last index refers to  
rᵤ          = ∂{P,N}(r1)[:,iu]       # the partial derivatives   
```

yields

```julia
julia> x1
2-element SVector{2, Muscade.∂ℝ{1, 5, Float64}} with indices SOneTo(2):
 ∂₁⟨1|1,0,0,0,0⟩ 
 ∂₁⟨2|0,1,0,0,0⟩ 

julia> u1
3-element SVector{3, Muscade.∂ℝ{1, 5, Float64}} with indices SOneTo(3):
 ∂₁⟨3|0,0,1,0,0⟩ 
 ∂₁⟨4|0,0,0,1,0⟩ 
 ∂₁⟨5|0,0,0,0,1⟩ 

julia> rₓ
2×2 SMatrix{2, 2, Float64, 4} with indices SOneTo(2)×SOneTo(2):
 2.0  0.0
 0.0  4.0

julia> rᵤ
2×3 SMatrix{2, 3, Float64, 6} with indices SOneTo(2)×SOneTo(3):
 1.0  1.0  1.0
 1.0  1.0  1.0 
``` 

## Directional derivatives

Derivatives with respect to a given linear combination `dir` of the input to a function or expression, can be computed as follows

```julia
using StaticArrays
using Muscade: value,∂,variate0

f(x)   = x.*x .+ sum(x)
x      = SVector(1.,2.,3.)
dir    = SVector(1.,0.,2.)
P,N,δs = variate0(0.)[1]
x1     = x + dir*δs
y1     = f(x1)
y      = value{P}(y1)
yₓ     = ∂{P,N}(y1)    
```

giving

```julia
julia> δs    
 ∂₁⟨0|1⟩ 

julia> x1    
3-element SVector{3, Muscade.∂ℝ{1, 1, Float64}} with indices SOneTo(3):
 ∂₁⟨1|1⟩ 
 ∂₁⟨2|0⟩ 
 ∂₁⟨3|2⟩ 
```

## Higher order derivatives

Higher order derivatives are obtained by using `variate{O}(...)` where `O` (defaulting to 1, as seen in previous examples) is the order of differentiation. 

```julia
using StaticArrays
using Muscade: variate,value,∂

lag(x,u) = x[1]*x[1]+u[2]*u[3]+x[2]*u[3]
x        = SVector(1.,2.)
u        = SVector(3.,4.,5.)
P,N,x2   = variate{2}(x)
ℓ2       = lag(x2,u)
ℓ        = value{P-1}(value{P}(ℓ2))
ℓₓ       = ∂{P-1,N  }(value{P}(ℓ2))    
ℓₓ       = value{P-1}(∂{P,N  }(ℓ2))    
ℓₓₓ      = ∂{P-1,N  }(∂{P,N  }(ℓ2))  
```

produces

```julia
julia> x2
2-element SVector{2, Muscade.∂ℝ{2, 2, Muscade.∂ℝ{1, 2, Float64}}} with indices SOneTo(2):
 ∂₂⟨∂₁⟨1|1,0⟩|∂₁⟨1|0,0⟩,∂₁⟨0|0,0⟩⟩
 ∂₂⟨∂₁⟨2|0,1⟩|∂₁⟨0|0,0⟩,∂₁⟨1|0,0⟩⟩

julia> ℓ2
∂₂⟨∂₁⟨31|2,5⟩|∂₁⟨2|2,0⟩,∂₁⟨5|0,0⟩⟩
```

As mentionned above `variate{2}` produces a nested structure, that can be loosely represented as `x2==∂₂⟨∂₁⟨x|xₓ⟩|∂₁⟨xₓ|xₓₓ⟩⟩`, and where `xₓ` is the Jacobian of `x` with respect to itself, and `xₓₓ` is, for each element of `x`, a matrix of zeros. `x2` has precedence `P==2`, as the variable name used here suggests.

The computations result in `ℓ2=∂₂⟨∂₁⟨ℓ|ℓₓ⟩|∂₁⟨ℓₓ|ℓₓₓ⟩⟩`, of precedence 2 "inherited" from `x2`. 

To unpack this nested structure, we first use `value{2}(ℓ2) → ∂₁⟨ℓ|ℓₓ⟩` and `∂{2,N}(ℓ2) → ∂₁⟨ℓₓ|ℓₓₓ⟩`.  On these, we then use: 

- `ℓ = value{1}(∂₁⟨ℓ|ℓₓ⟩)` 
- `ℓₓ = ∂{1,N}(∂₁⟨ℓ|ℓₓ⟩)` 
- `ℓₓ = value{1}(∂₁⟨ℓₓ|ℓₓₓ⟩)`
- `ℓₓₓ = ∂{1,N}(∂₁⟨ℓₓ|ℓₓₓ⟩)` 

`ℓₓ` is stored twice in `ℓ2` and so can be extracted by two different paths.  

A caution about performance: higher derivatives means that variated numbers carry a significant amount of information, and all operations on them become slower.  Also, because the compiler rolls out all the loops in operations with variated numbers, compilation can become very heavy. 

## Nested derivatives

Nested derivatives typicaly occur when a function like `Muscade.residual`, which is called by `Muscade` with automatic differentiation, itself uses automatic differentiation internaly. An example in `Muscade.Toolbox` is the implementation of `EulerBeam3D`. This beam element uses a co-rotated reference system, which poses the problem of transforming, for example, bending moments `mᵢ` at Gauss quadrature points, back into nodal loads `R`.  Mathematicaly, their contribution to the residual is `mᵢ ∘₁ κ∂X₀`, where `κ∂X₀` is the Jacobian of the curvatures `κ` with respect to the 0-th order time derivative `X₀` of the element's `X`-dofs.

Here is a simple (if contrived) example:

```julia
using StaticArrays
using Muscade: variate,value,∂,𝕣

function f(x1) 
    a1     = x1.^2
    Pa,Na,a2 = variate(a1)         # P=2
    g2     = sum(a2)*sin(sum(a2))
    gₐ     = ∂{Pa,Na}(g2) 
    return x.*x + sum(u).*gₐ 
end

x      = SVector(1.,2.,3.)
Px,Nx,x1 = variate(x)               # P=1
y1     = f(x1)
y      = value{Px}(y1)
yₓ     = ∂{Px,Nx}(y1)    
```

The function `f` calls `variate`, on `x1` which has itself been obtained by calling `variate`.  `x`, `x1` and `a2` have precedence 0, 1 and 2, respectively. The unpacking of the variables is as for higher order derivatives, but here, each segment of code takes care to unpack what it has variated.

`a2` has structure `a2==∂₂⟨∂₁⟨a|aₓ⟩|∂₁⟨aₐ|aₓₐ⟩⟩`.  Note that `a2` now stores two distinct first order derivatives.

## Understanding and avoiding perturbation confusion

Nested differentiation in all its forms carries the danger of "perturbation confusion". 

Here is an example of how things go wrong. Variable names include the precedence. This is only done here to help readability. When implementing an element in `Muscade` the precedence of input variables is not known to the developer, because `residual` and `lagrangian` are called by various `Muscade`-solvers with input variables `Λ`, `X`, `U`, `A` and `t` that may or may not be variated, and that to different orders.

```julia
using StaticArrays
using Muscade: variate,value,∂,𝕣

function residual(x1,u) 
    Pu,Nu,u1 = variate(u)           # Pu=Px !!!
    g1       = sum(x1)*sin(sum(u1)) # crash
    gᵤ       = ∂{Pu,Nu}(g1) 
    return sum(x1.*x1) + sum(u.*gᵤ) 
end

x        = SVector(1.,2.,3.)
u        = SVector(4.,5.)
Px,Nx,x1 = variate(x)               # Px=1
r1       = residual(x1,u)
r        = value{Px}(r1)
rₓ       = ∂{Px,Nx}(r1)    
```

The above, *luckily*, crashes. In `residual`, we have `u1[1]==∂₁⟨4|1,0⟩` and `x1[1]==∂₁⟨1|1,0,0⟩`. The problem is that `u1` contains partial derivatives with respect to `u` while `x1` contains partial derivatives with respect to `x`.  However `Muscade`'s implementation of automatic differentiation assumes that since both variables have the same precedence 1 (as marked by the subscript in `∂₁`), they both refer to partial derivatives with respect to the same variable. Mayhem follows when attempting to compute linear combinations of the two vectors of partial derivatives.  A similar counter example, with `u` of length 3 leads to an even more undesirable situation, since the "perturbation confusion" then can go undetected.

To avoid these problems, the call to `variate` in the definition of `residual` must be rewriten:

```julia
using StaticArrays
using Muscade: variate,value,∂,𝕣

function residual(x1,u) 
    Pu,Nu,u2 = variate(u,context=x1) # Pu=Px+1 
    g2       = sum(x1)*sin(sum(u2))  # (u1 renamed u2)
    g1ᵤ      = ∂{Pu,Nu}(g2) 
    return sum(x1.*x1) + sum(u.*g1ᵤ) 
end

x        = SVector(1.,2.,3.)
u        = SVector(4.,5.)
Px,Nx,x1 = variate(x)                # Px==1
r1       = residual(x1,u)
r        = value{Px}(r1)
rₓ       = ∂{Px,Nx}(r1)    
```

Thanks to the `context` input, `variate` can now detect that although `u` has precedence 0, `x1` already "occupies" precedence 1.  Variate thus returns `u2`, of precedence `u==2`. `u2` has the structure `∂₂⟨∂₁⟨u|uₓ⟩|∂₁⟨uᵤ|uₓᵤ⟩⟩`, where `uₓ` and `uₓᵤ` are zero (`u` does not depend on `x`).  `x1` and `u2` operate correctly with each other, in spite of the different precedences.

When implementing an element, the safest is to set `context=(Λ,X,U,A,t)` (`lagrangian`) or `context=(X,U,A,t)` (`residual`). 

But if, for example, a solver was calling `residual` with only `A` variated, and a `∂ᵢ` issued by a call to `variate` within `residual` never touches `A` or any variable influenced by `A`, then omitting `A` from the `context` allows `variate` to issue a variable of lower precedence, enabling faster execution. _Be very deliberate about this_:

!!! important
    Never allow two `∂ᵢ` issuing from different calls to [`variate`](@ref), [`variate0`](@ref), [`revariate`](@ref) and [`motion`](@ref) to mix.

Consider using `let` blocks to create a local scope to contain "dangerous variables":

```julia
function f(x) 
    gₓ = let P,N,x_ = variate(x)
        g_ = sum(x_)*sin(sum(x_))
        ∂{P,N}(g_) 
    end
    x.*x + sum(u).*gₓ 
end
```

## Time derivatives

Dynamic solvers call `Muscade.residual` and `Muscade.lagrangian` with `X` and `U` that are `NTuple`s of `SVector`s.  The elements of the `NTuple`s represent the time derivatives of the element's degrees of freedom.

Finite elements solve differential equations by introducing interpolations schemes to represent fields within the element as a function of the degrees of freedom.  For example, this could be a relation of the type `ε₀ = kinematics(X₀)`, where `X₀` and `ε₀` are the nodal displacements and the strain at integration points.  In some elements, `kinematics` can be non-linear and complicated to differentiate. [`Toolbox.EulerBeam3D`](@ref) provide a good example of this. In some material models, the time derivative of `ε₁` can be required. The function `motion` uses directional derivatives to transform the `NTuple` `X` into nested `∂ℝ` `X_`. `X_` is used to call `kinematics`, and `motion⁻¹` used to recover a tuple containing the time derivatives of `ε`. This works whether `X` is already variated or not.

```julia
kinematics(X₀) = sum(X₀.*X₀)
X₀,X₁,X₂       = (SVector(1.),SVector(2.),SVector(3.))
X              = (X₀,X₁,X₂) # argument to residual
# in residual:
(P,ND,X_)      = motion(X) 
ε_             = kinematics(X_) 
ε              = motion⁻¹{P,ND}(ε_)
ε₀, ε₁         = ∂0(ε), ∂1(ε)
```

Like `variate`, `motion` has an optional keyword argument `context` and the same rules of caution apply. On the other hand, the above example is quite typical: `kinematics` is a function of `X_` only, and there is no need to extend the context.

If `x` (lower case) is the i-th degree of freedom, then `X_[i]` has the structure  `∂₂⟨∂₁⟨x₀|x₁⟩|∂₁⟨x₁|x₂⟩⟩`, where `xᵢ` is the i-th time derivative.  Hence the number of partials is 1 for both `∂₂` and `∂₁`: while higher order derivatives in general are expensive to compute, `motion` tends to be affordable. 

If `X₀`,`X₁`,`X₂` or `context` contains variables of precedence for example 1, then `motion` returns `∂₃⟨∂₂⟨x₀|x₁⟩|∂₂⟨x₁|x₂⟩⟩`.

## Chain rule

Applying the chain rule *can*, under the right circumstances accelerate computations. We are differentiating a function `f` with respect to a long vector `x`. Within `f`, a computationaly expensive subfunction `z=g(y)` is only function of a significantly shorter vector `y`.

Instead of traversing the function `g` with a large number of partials with respect to `x`, one can instead evaluate the derivatives of `z` with respect to `y`, and, having the derivatives of `y` with respect to `x`, recover the derivatives of `z` with respect to `x`.  

When differentiating to the first order, this is the chain rule
``\frac{\partial z_i}{\partial x_k} = \sum_j \frac{\partial z_i}{\partial y_j}\cdot\frac{\partial y_j}{\partial x_i}``. When diffrentiating to an arbitrary order, the computation `z=g(y)` yields a variated number `z` that contains value and multiple derivatives. This data is used to evaluate a multi-dimensional Taylor expansion of `g`, with respect to `y-VALUE(y)`.

```julia
using StaticArrays
using Muscade: variate,apply

g(y)  = sin.(y).*y              # computationaly expensive
x     = SVector(1.,2.,3.)       # long x
P,N,X = variate(x)              # N is large
Y     = 3. .+ X[1] .+ X[2]*X[3] # short y
# Z = g(Y)                      # too slow
Z     = apply{:chainrule}(g,Y)  # faster (well, sometimes)
```

Looking under the hood, the code of `apply` is

```julia
function apply{:chainrule}(f,x) 
    x_ = revariate(x)           # with the code above, N=1
    y_ = f(x_)                  # this is accelerated, as wanted
    y  = chainrule(y_,x)        # evaluating the Taylor expansion can be slow
    return y
end
apply{:direct}(f,x)=f(x)
```

Try `apply{:chainrule}` vs. `apply{:direct}` to see which is faster.
