# Creating an element

## Introducton 

New element types, providing new element formulations or adressing new domains in physics are implemented by implementing a datatype and methods for function defined by `Muscade`: 

- A **constructor** which is called when the user adds an element to the model.
- **`Muscade.doflist`** specifies the degrees of freedom (dofs) of the element.
- **`Muscade.residual`** takes element dofs as input and returns the element's additive contribution to a non-linear system of equations to be solved (a residual aka. a right hand side),
- **`Muscade.lagrangian`** takes element dofs as input and returns the element's additive contribution to a function to be extremized (Lagrangian aka. cost function, surprisal...)

Each element must implement a constructor following a format specified by `Muscade`, as well as `Muscade.doflist`.  Each element must also implement *either* `Muscade.lagrangian` *or* `Muscade.residual`, depending on what is more natural.

## DataType
For a new element type `MyELement`, the datatype is defined as

```julia
struct MyElement <: AbstractElement
    ...
end
```

`MyElement` *must* be a subtype of `AbstractElement`.


## Constructor

The element must provide a constructor of the form

```julia
function MyElement(nod::Vector{Node};kwargs...)
    ...
    return eleobj
end
```

which will then call the default constructor provided by Julia.  

`nod` can be used to access the coordinates of the nodes:

```julia
x = nod[inod].coord[icoord]
```

where `inod` is the element-node number and `icoord` the index into a vector of coordinates. `coord` is 
provided by the user when adding a `Node` to the `Model`. `Muscade` has no opinion about, and provides no
check of, the length of `coord` provided by the user. In this way elements can define what coordinate system
(how many coordinates, and their interpretation) is to be used.  Coordinate systems can even differ from one
node to the next. See also [`coord`](@ref).

`kwargs...` is any number of named arguments, typicaly defining the material properties of the element.

The user does not call the above-defined element directly.  Instead, an element is added to the model by
a call of the form

```julia
e1 = addelement!(model,MyElement,nodid,kwargs...)
```

See [`addelement!`](@ref).

## `Muscade.doflist`

The element must provide a method of the form

```julia
function Muscade.doflist(::Type{MyElement})
    return (inod =(...),
            class=(...),
            field=(...))
end
```

The syntax `::Type{MyElement}` is because `doflist` will be called by `Muscade` with a `DataType` (the type `MyElement`). The function name must begin
with `Muscade.` to make it possible to overload a function defined in the module `Muscade`. 

The return value of the function is a `NamedTuple` with the fields `inod`, `class` and `field`. 

- `inod`  is a `NTuple` of `Int64`: for each dof, its element-node number. 
- `class` is a `NTuple` of `Symbol`: for each dof, its class (must be `:X`, `:U` or `:A`).
- `field` is a `NTuple` of `Symbol`: for each dof, its field.

Importantly, `doflist` does not mention dofs of class `:Λ`: if the element implements `lagrangian`, there is automaticaly a one-to-one correspondance between Λ-dofs and X-dofs.

For example (using Julia's syntax for one-liner functions):

```julia
Muscade.doflist( ::Type{Turbine}) = (inod =(1   ,1   ,2        ,2        ),
                                     class=(:X  ,:X  ,:A       ,:A       ),
                                     field=(:tx1,:tx2,:Δseadrag,:Δskydrag))
```

See [`doflist`](@ref).

## `Muscade.lagrangian`

Elements that implement a cost on the degrees of freedom must implement a method of the form

See [`lagrangian`](@ref) for the list of arguments and outputs.

### Automatic differentiation

The gradients and Hessians of `L` do not need to be implemented, because `Muscade` uses automatic 
differentiation. Because of this, it is important not to over-specify the inputs.  For example, 
implementing a function header with

```julia
@espy function Muscade.lagrangian(o::MyElement,Λ::Vector{Float64},X,U,A,t,χ,χcv,SP,dbg)
#                                               ____bad_idea____
```

would cause a `MethodError`, because `Muscade` will attempt to call with a `SVector` instead of `Vector`, and a special
datatype supporting automatic diffeentiation instead of `Float64`.

### Extraction of intermediate element results

The function definition must be anotated with the macro call `@espy`.  Variables within the body of `lagrange`, which
the user may want to obtain must be anotated with `☼` (by typing `\sun` they pressing `TAB`) at the place where they
are calculated. An example would be 

```julia
    ☼σ = E*ε
```

The macro will generate two versions of `lagrange`.  One in which the anotations `☼` are taken away, which is used to
solve the numerical problem.  Another with additional input and output variables, and code inserted into the body
of the function to extract results wanted by the user.

See [`@espy`](@ref) for a complete guide on code anotations. 


## `Muscade.residual`

Elements that implement "physics" will typicaly implement `residual` (they could implement the same using `lagrangian`, but the resulting code would be less performant).

The interface is mostly the same as for `lagrangian` with the differences that

- `residual` returns a vector `R`
- there is no argument `Λ`

```julia
@espy function Muscade.residual(o::MyElement,X,U,A,t,χ,χcv,SP,dbg) 
    ...
    return R,noχ,noFB
end
```

### Immutables and Gauss quadrature

`residual` is called many times and it is critical to obtain high performance. Thus, allocating the vector `R` on the heap 
within the function must be avoided. A design option would be to have `Muscade` pass a preallocated vector `R` and have 
residual mutate its argument.  For "forward" automatic differentiation, it can be difficult to predict the element type of `R`,
and other techniques of automatic differentiaton do not accomodate mutations.

For this reason, `residual` and `lagrange` must be writen in a functional style, using only immutable variables, and in particular
immutable arrays.  This can be done using `StaticArrays.jl` (tested) or `Tensorial.jl` (not tested with `Muscade`), and generaly
results in very readable code that directly expresses concepts of linear algebra.

One difficulty arises with Gauss quadrature.  Typical implementations would rely on setting `R` to zero, then adding the 
contributions from quadrature points to `R` within a `for` loop over the Gauss points.  The pseudo code (not valid in `Muscade`):

```julia
R .= 0
for igp = 1:ngp
    F = ...
    Σ = ...
    R += F ∘ Σ ∘ ∇N * dV
end
```

shows that `R` is mutated. A pseudocode in immutable style would be

```julia
@espy function residual(x,χ)
    t = ntuple(ngp) do igp
        ☼F = ...
        ☼Σ = ...
        r = F ∘ Σ ∘ ∇N * dV
        @named(χ,r)
    end
    χ = ntuple(igp->t[igp].χ,ngp)
    R = sum(   igp->t[igp].r,ngp)
    return R,χ
end
```

`r` are the contributions to `R` at each quadrature point.  The operation `t = ntuple ...` returns a datastructure `t` such that `t[igp].χ` are the memory
variable and ``t[igp].r` the contribution to the residual from the `igp`-th quadrature point. This is because

```julia
    t = ntuple(ngp) do igp
        expr(igp)
    end    
```

is equivalent to 

```julia
    t = ntuple(expr for igp=1:ngp)
```

which returns

```julia
    t = (expr(1),expr(2),...,expr(ngp))
```

where the value of `expr` is that of its last line `@named(χ,r)` which is a macro provided by `Muscade` that inserts the code `(χ=χ,r=r)`.

The code 

```julia
    χ = ntuple(igp->t[igp].χ,ngp)
    R = sum(   igp->t[igp].r,ngp)
```

gathers the memories of all quadrature points into a `Tuple` and adds together the contributions `r` into the residual `R`.

See [`residual`](@ref).

## Help functions

`Muscade` provides functions and constants to make it easier to comply with the API:

- Element constructors can use function [`coord`](@ref) to extract the coordinates fron the `Vector{Node}` they get as first argument.
- `residual` and `lagrangian` **must** use [`∂0`](@ref), [`∂1`](@ref) and [`∂2`](@ref) when extracting the zeroth, first and second time derivatives from arguments `X` and `U`.
- Constants [`noχ`](@ref) and [`noFB`](@ref) (which have value `nothing`) can be used by elements that do not have memory or no feedback to the solving procedure.

## Performance

For a given element formulation, the performance of `residual` and `lagrangian` can vary with a factor up to 100 between a good and a bad implementation.

**Type stable code** allows the compiler to know the type of every variable in a function given the type of its parameters. Code that is type unstable is significantly slower. See the page on [type stability](TypeStable.md).

**Allocation**, and the corresponding deallocation of memory *on the heap* takes time. By contrast, allocation and deallocation *on the stack* is fast.  In Julia, only immutable variables can be allocated on the stack. See the page on [memory management](Memory.md)

**Automatic differentiation** generaly does not affect how `residual` and `lagrangian` are writen.  There are two performance-related exceptions to this:

1. If a complex sub function in `residual` and `lagrangian` (typicaly a material model or other closure) operates on an array (for example, the strain) that is smaller than the number of degrees of freedom of the system, computing time can be saved by computing the derivative of the output (in the example, the stress) with respect to the input to the subfunction, and then compose the derivatives.
2. Iterative precedures are sometimes used within `residual` and `lagrangian`, a typical example being in plastic material formulations.  There is no need to propagate automatic differentiation through all the iterations - doing so with the result of the iteration provides the same result.

See the page on [automatic differentiation](Adiff.md)





