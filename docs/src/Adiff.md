# Automatic differentiation

## History

`Muscade` has its own implementation of forward automatic differentiation for historical reasons: Prototypes automatic differentiation of `Muscade` where developed in parallel with [`ForwardDiff.jl`](https://juliadiff.org/ForwardDiff.jl/stable/).  While the inner workings of [`ForwardDiff.jl`](https://juliadiff.org/ForwardDiff.jl/stable/) and `Muscade`'s version are similar ( with `ForwardDiff.jl` probably having better performance), the API are quite different.

`Muscade` evaluates second derivative of the Lagrangian, using nested *forward* differentiation, which is far from optimal.  An ambition is to make use of reverse differentiation (using [`Zygote.jl`](https://fluxml.ai/Zygote.jl/latest/), [`Enzyme.jl`](https://docs.sciml.ai/Enzyme/stable/) or similar).

## Usage

`Muscade`s automatic differentiation is used as follows:

```julia
using SVector
x   = SVector(1.,2.,3.)
N   = length(x)
x1  = variate{1,N}(x)
y1  = f(x1)
y   = value{1}(y1)
yₓ  = ∂{1,N}(y2)    
```

`x` *must* be a `SVector`. In `yₓ`, the index over `x` is post-pended to the indices of `y`.

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

One way to accelerate automatic differentiation of complicated functions can be to chainrule the differentiation of simpler functions, in particular function with smaller inputs.  Use [`fast`](@ref) for basic applications, and [`revariate`](@ref) and [`chainrule`](@ref) for advanced usage. 
