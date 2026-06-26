
# [Automatic differentiation](@id adiff_appendix)

## Introduction
This section aims to provide an understanding of a variety of usefull techniques in automatic differentiation. In Muscade, element developper generaly do not need to be concerned about automatic differentiation, except on two points:

First, do not specify the type of inputs `Λ`, `X`, `U`, `A` and `t`.  For example

```julia
@espy function Muscade.lagrangian(o::MyElement,Λ::Vector{Float64},X,U,A,t,SP,dbg)
#                                               |____bad_idea___|
```

would cause a `MethodError`, because `Muscade` will attempt to call `lagrangian` with a `SVector{L,∂ℝ...}`, and the method only accepts a`Vector{𝕣}`.

Second, beware of branches (`if` and `cond ? opt1 : opt2`).  In the following code

```julia
if x>0
    y=x  # y is of type ∂ℝ
else
    y=0. # y is of type 𝕣
end    
```

the type of `y` depends on the *value* of `x` and can thus not be determined at compilation.  This is [type instability](TypeStable.md) and results in *very* slow execution.

What follows is a description of a range of techniques only needed in some demanding situations.

`Muscade` has its own implementation of forward automatic differentiation for historical reasons: Prototypes automatic differentiation of `Muscade` where developed in parallel with [`ForwardDiff.jl`](https://juliadiff.org/ForwardDiff.jl/stable/).  The inner workings of [`ForwardDiff.jl`](https://juliadiff.org/ForwardDiff.jl/stable/) and `Muscade`'s version are very similar. 

The APIs, on the other hand, are quite different: [`ForwardDiff.jl`](https://juliadiff.org/ForwardDiff.jl/stable/) is designed for safety, `Muscade`'s automatic differentiation for flexibility. When using the techniques described in the following, it is possible to introduce some confusing bugs if one is not quite careful.

`Muscade` also offers some innovative differentiation functionality (see the Section on chain rule).

## Type definition
`Muscade`'s data-type for automatic differentiation, `∂ℝ`, is defined as

```julia
struct ∂ℝ{P,N,R} <:ℝ where{R<:ℝ}  
    x  :: R
    dx :: SVector{N,R}
end
```

where `ℝ` is `Muscade` shorthand for `Real`.  The type parameters of `∂ℝ` is

- `P` is the *precedence* (to prevent [perturbation confusion](https://arxiv.org/pdf/1211.4892v3), discussed in the following).
- `N` is the number of partial derivatives.
- `R` is the underlying datatype (either `Float64`, noted `𝕣` in Muscade, or another `∂ℝ` when computing higher order derivatives).

## Techniques and tools

**Gradient evaluation** is the simplest form for automatic differentiation.

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

[`variate`](@ref) builds a `SVector` of `∂{P,N,R}`, where `R == eltype(x)`. `x` *must* be a `SVector{N}`. `Muscade` developers speak of `x1` being the "adiffed" form of `x`". Once the calculations to be differentiated are completed, results (here:`y`) are `StaticArrays` of `∂{P,N,R}`.  The partials are "seeded", so that `∂{P,N}(x1)[i]` is a `SVector{N}` of `zero(R)`, except for the `i`-th term that is `one(R)`. 

```julia
julia> x1
3-element SVector{3, Muscade.∂ℝ{1, 3, Float64}} with indices SOneTo(3):
 1+∂₁⟨1,0,0⟩ 
 2+∂₁⟨0,1,0⟩ 
 3+∂₁⟨0,0,1⟩ 
```

where `∂₁` can be thought of as a number following all the rules of algebra, but with `∂₁⋅∂₁ = 0`, and `⟨...⟩` delimits a vector of partials.

```julia
julia> y1
3-element SVector{3, Muscade.∂ℝ{1, 3, Float64}} with indices SOneTo(3):
  7+∂₁⟨3,1,1⟩ 
 10+∂₁⟨1,5,1⟩ 
 15+∂₁⟨1,1,7⟩ 
 ``` 

[`value`](@ref)`{P}`, where `y1` is a `StaticArray{S,∂{P,N,R}}`, retrieves the value `y`, of type `StaticArray{S,R}`. [`∂`](@ref)`{P,N}`, where `y1` is a `StaticArray{S,∂{P,N,R}}`, retrieves the gradient `yₓ`, of type `StaticArray{Syₓ,R}`: `yₓ` has one more dimension than `y`, of size `N`, corresponding to the pratial derivatives with respect to the `SVector` `x`.

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

**Joint derivatives** are partial derivatives with respect to various variables, computed at the same time. An example of the technique is

```julia
using StaticArrays
using Muscade: variate,value,∂,δ,𝕣

res(x,u) = x.*x .+ sum(u) # SVector-valued 
x        = SVector(1.,2.)
u        = SVector(3.,4.,5.)
P,N,R    = 1,length(x)+length(u),𝕣
ix       = SVector(1,2)
iu       = SVector(3,4,5)
δ1       = Muscade.Muscade.δ_{P,N,R}()
x1       = x + δ1[ix]
u1       = u + δ1[iu]
r1       = res(x1,u1)
r        = value{P}(r1)
rₓ       = ∂{P,N}(r1)[:,ix]    
rᵤ       = ∂{P,N}(r1)[:,iu]    
```

yields

```julia
julia> δ1
5-element SVector{5, Muscade.∂ℝ{1, 5, Float64}} with indices SOneTo(5):
 0+∂₁⟨1,0,0,0,0⟩ 
 0+∂₁⟨0,1,0,0,0⟩ 
 0+∂₁⟨0,0,1,0,0⟩ 
 0+∂₁⟨0,0,0,1,0⟩ 
 0+∂₁⟨0,0,0,0,1⟩ 

julia> x1
2-element SVector{2, Muscade.∂ℝ{1, 5, Float64}} with indices SOneTo(2):
 1+∂₁⟨1,0,0,0,0⟩ 
 2+∂₁⟨0,1,0,0,0⟩ 

julia> u1
3-element SVector{3, Muscade.∂ℝ{1, 5, Float64}} with indices SOneTo(3):
 3+∂₁⟨0,0,1,0,0⟩ 
 4+∂₁⟨0,0,0,1,0⟩ 
 5+∂₁⟨0,0,0,0,1⟩ 
``` 

`Muscade.δ_` is similar to `variate`, except that all the `value`s of `δs` are `zero(R)`.  In the above the elements of `x1`, `u2` and `r1` all have 5 partials - the first two corresponding to `x` and the rest to `u`. 

The above code is quite clumsy, and the function [`revariate`](@ref) allows to do the same more elegantly.

**Directional derivatives**, that is, derivatives with respect to a given linear combination of the input to a function or expression, can be computed as follows

```julia
using StaticArrays
using Muscade: variate,value,∂,Muscade.δ_,𝕣

f(x)  = x.*x .+ sum(x)
x     = SVector(1.,2.,3.)
dir   = SVector(1.,0.,2.)
P,N,R = 1,1,𝕣
δs    = δ_{P,N,R}()[1]
x1    = x + dir*δs
y1    = f(x1)
y     = value{P}(y1)
yₓ    = ∂{P,N}(y1)    
```

```julia
julia> δs    = Muscade.δ_{P,N,R}()[1]
0+∂₁⟨1⟩ 

julia> x1    = x + dir*δs
3-element SVector{3, Muscade.∂ℝ{1, 1, Float64}} with indices SOneTo(3):
 1+∂₁⟨1⟩ 
 2+∂₁⟨0⟩ 
 3+∂₁⟨2⟩ 
```

**Cross derivatives** provide the best way to understand the mesting of `∂ℝ` variables can be nested.

```julia
using StaticArrays
using Muscade: variate,value,∂

lag(x,u) = x[1]*x[1]+u[2]*u[3]+x[2]*u[3]
x        = SVector(1.,2.)
u        = SVector(3.,4.,5.)
Px,Nx    = 1,length(x)
Pu,Nu    = 2,length(u) # NB: Pu > Px
x1       = variate{Px,Nx}(x)
u2       = variate{Pu,Nu}(u)
ℓ2       = lag(x1,u2) # \ell
ℓ        = value{Px}(value{Pu}(ℓ2))
ℓₓ       = ∂{Px,Nx }(value{Pu}(ℓ2))    
ℓᵤ       = value{Px}(∂{Pu,Nu }(ℓ2))    
ℓᵤₓ      = ∂{Px,Nx }(∂{Pu,Nu }(ℓ2))  
```

We now use variations of different precedence.  The structure of `ℓ2` is `ℓ+∂₁⟨ℓₓ⟩+∂₂⟨ℓᵤ+∂₁⟨ℓᵤₓ⟩⟩`, and its type `∂ℝ{Pu,Nu,∂ℝ{Px,Nx,𝕣}}`.
 Note the presence of both `∂₁` and `∂₂`, and how how `∂₂` wraps outside of `∂₁`. The previously-mentionned rule generalises to ``∀i, ∂ᵢ⋅∂ᵢ = 0`` (and only that: `∀(i,j), i≠j ⇒ ∂ᵢ⋅∂ⱼ ≠ 0`). 

```julia
julia> x1
2-element SVector{2, Muscade.∂ℝ{1, 2, Float64}} with indices SOneTo(2):
 1+∂₁⟨1,0⟩ 
 2+∂₁⟨0,1⟩ 

julia> u2
3-element SVector{3, Muscade.∂ℝ{2, 3, Float64}} with indices SOneTo(3):
 3+∂₂⟨1,0,0⟩ 
 4+∂₂⟨0,1,0⟩ 
 5+∂₂⟨0,0,1⟩  

julia> ℓ2
31+∂₁⟨2,5⟩+∂₂⟨0+∂₁⟨0,0⟩,5+∂₁⟨0,0⟩,6+∂₁⟨0,1⟩⟩ 

julia> ℓₓ 
2-element SVector{2, Float64} with indices SOneTo(2):
 2.0
 5.0

julia> ℓᵤ
3-element SVector{3, Float64} with indices SOneTo(3):
 0.0
 5.0
 6.0

 julia> ℓᵤₓ 
3×2 SMatrix{3, 2, Float64, 6} with indices SOneTo(3)×SOneTo(2):
 0.0  0.0
 0.0  0.0
 0.0  1.0
```

Note the specific order in which `ℓ2` is unpacked to obtain the derivatives: first precedence `2` is unpacked (because `Muscade.precedence(ℓ2) == 2`), then precedence `1`.  Unpacking `ℓ2` directly, with a precedence lower than `2` causes an error, as it amounts to "peeling the onion from the inside".

Unpacking `ℓ2` with any precedence higher than `2` will cause `ℓ2` to be treated as a constant: `value{3}(ℓ2) == ℓ2`, and `∂{3,6}(ℓ2)` returns an array of zeros. Similarly, `value{Pu}(x1) == x1` and `x1ᵤ₁ == ∂{Pu,Nu }(x1)` is zeros.

**Higher order derivatives** are a special case of cross derivatives (the cross derivative of a variable with itself), but this requires nested calls to `variate`:

```julia
using StaticArrays
using Muscade: variate,value,∂

lag(x,u) = x[1]*x[1]+u[2]*u[3]+x[2]*u[3]
x        = SVector(1.,2.)
u        = SVector(3.,4.,5.)
P1,P2,N  = 1,2,length(x)
x1       = variate{P1,N}(x )
x2       = variate{P2,N}(x1)
ℓ2       = lag(x2,u)
ℓ        = value{P1}(value{P2}(ℓ2))
ℓₓ       = ∂{P1,N  }(value{P2}(ℓ2))    
ℓₓ       = value{P1}(∂{P2,N  }(ℓ2))    
ℓₓₓ      = ∂{P1,N  }(∂{P2,N  }(ℓ2))  
```

Note in particular the nested calls

```julia
x1       = variate{P1,N}(x )
x2       = variate{P2,N}(x1)
```

in which the lower precedence `P1` must applied before the higher precedence `P2`, and the fact that `ℓₓ` can be obtained in two different ways.  Indeed, the values are stored twice in `ℓ2` which has structure `ℓ+∂₁⟨ℓₓ⟩+∂₂⟨ℓₓ+∂₁⟨ℓₓₓ⟩⟩`.

```julia
julia> x2
2-element SVector{2, Muscade.∂ℝ{2, 2, Muscade.∂ℝ{1, 2, Float64}}} with indices SOneTo(2):
 1+∂₁⟨1,0⟩+∂₂⟨1+∂₁⟨0,0⟩,0+∂₁⟨0,0⟩⟩ 
 2+∂₁⟨0,1⟩+∂₂⟨0+∂₁⟨0,0⟩,1+∂₁⟨0,0⟩⟩ 

julia> ℓ2
31+∂₁⟨2,5⟩+∂₂⟨2+∂₁⟨2,0⟩,5+∂₁⟨0,0⟩⟩ 
```

**Nested derivatives** occur when a function, which uses automatic differentiation internaly is itself called with automatic differentiation. An example in `Muscade.Toolbox` is the implementation of `EulerBeam3D`'s method for `Muscade.residual`. This beam element uses a co-rotated reference system, which pauses the problem of how to transform, for example, bending moments `mᵢ` at Gauss quadrature points, back into nodal loads `R`.  Mathematicaly, their contribution to the residual is `mᵢ ∘₁ κ∂X₀`, where `κ∂X₀` is the Jacobian of the curvatures `κ` with respect to the 0-th order time derivative `X₀` of the element's `X`-dofs.

Here is a simple example:

```julia
using StaticArrays
using Muscade: variate,value,∂,δ_,𝕣,npartial,constants

function f(x) 
    P  = constants(x)
    N  = length(x)
    xP = variate{P,N}(x)
    gP = sum(xP)*sin(sum(xP))
    gₓ = ∂{P,N}(gP) 
    x.*x + sum(u).*gₓ 
end

x    = SVector(1.,2.,3.)
P,N  = 1,length(x)
x1   = variate{P,N}(x)
y1   = f(x1)
y    = value{P}(y1)
yₓ   = ∂{P,N}(y1)    
```

In a context where the programmer of `f` cannot know to what order the input variable `x` is already adiffed, [`constants`](@ref) provides a precedance that is higher than the precedence of `x` (`constants` can handle multiple variables), to prevent [perturbation confusion](https://arxiv.org/pdf/1211.4892v3).

!!! warning
     Any pair og variables of identical precedence (and this applies recursively, where higher order differentiation is used), but issued from different factories (including [`variate`](@ref), [`Muscade.δ_`](@ref), [`revariate`](@ref) and [`motion`](@ref)) must never be allowed to be the two arguments of a binary operation.  Failure to enforce this may cause silent failure, producing erroneous results.  

Consider using `let` blocks to contain "dangerous variables".

```julia
function f(x) 
    P  = constants(x)
    N  = length(x)
    gₓ = let xP = variate{P,N}(x)
        gP = sum(xP)*sin(sum(xP))
        ∂{P,N}(gP) 
    end
    x.*x + sum(u).*gₓ 
end
```

A caution about performance: nesting derivatives (and any form of differentiation to higher order) means that addif'ed numbers carry a significant amount of information, and all operations on them become slower (linearly with the amount fo information). 

**Time derivatives** appear in dynamic solvers. `Muscade.residual` and `Muscade.lagrangian` receive inputs `X` and `U` that are `ntuples` of `SVectors`.  The elements of the `ntuples` represent the time derivatives of the degrees of freedom.

Finite elements solve differential equations by introducing interpolations schemes to represent fields within the element as a function of the degrees of freedom.  For example, this could be a relation of the type `ε₀ = kinematics(X₀)`, where `X₀` and `ε₀` are the nodal displacements and the strain.  In some cases, `kinematics` can be non-linear and complicated to differentiate. [`Toolbox.EulerBeam3D`](@ref) provide a good example of this. In some cases, the time derivative of `ε₁` can be required. The function `motion` uses directional derivatives to transform the `ntuple` `C` into nested `∂ℝ`. These are used to call `kinematics`, and `motion⁻¹` used to recover a tuple containing the time derivatives of `ε`. This works wether `X` is already adiffed or not.

```julia
    ...
    P  = constants(X,U,A,t)
    N  = length(X)
    ε  = let X_  = motion{P,N}(X)
        ε_  = kinematics(X_) 
        motion⁻¹(ε_)
    end   
# or: ε  = motion⁻¹(kinematics(motion{P,N}(X))) 
    ε₀, ε₁ = ∂0(ε), ∂1(ε)
```

For example

```julia
julia> X=(SVector(1.),SVector(2.),SVector(3.));

julia> Muscade.motion{1}(X)
1-element SVector{1, ∂ℝ{2, 1, ∂ℝ{1, 1, Float64}}} with indices SOneTo(1):
 1+∂₁⟨2⟩+∂₂⟨2+∂₁⟨3⟩⟩
```

**Chain rule** is a technique that *can*, under the right circumstances accelerate computations. We are differentiating a function `f` with respect to a long vector `x`. A computationaly expensive subfunction `z=g(y)` of `f` is only function of a significantly shorter vector `y`.

Instead of traversing the function `g` with a large number of partials with respect to `x`, one can instead evaluate the derivatives of `z` with respect to `y`, and, having the derivatives of `y` with respect to `x`, recover the derivatives of `z` with respect to `x`.  When differentiating to the first order, this is the chain rule
``\frac{\partial z_i}{\partial x_k} = \sum_j \frac{\partial z_i}{\partial y_j}\cdot\frac{\partial y_j}{\partial x_i}``.

When differentiating to a higher order, the computation `z=g(y)` yields a data structure `z` that contains value and multiple derivatives. This data can be used to evaluate a Taylor expansion of `g`, with respect to `y-VALUE(y)`.

In this example `x` is not actualy long, to make outputs more readable.

```julia
using StaticArrays
using Muscade: variate,VALUE,McLaurin

g(y) = sin.(y).*y                   # computationaly expensive

x = variate{3,1}(SVector(1.))       # long x
y = 3. .+ x .+ x.^2                 # short y

# z = g(y)                          # too slow

z = let y_ = variate{3,1}(VALUE(y)) # fewer partials
    z_ = g(y_)                      # faster execution
    McLaurin(z_,y-VALUE(y))         # ah! Not cheap
end
```
