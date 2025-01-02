# [Implementing new elements](@id implementelement)

## Introduction

The implementation of a element requires 

- A `struct` defining the element. 
- A **constructor** which is called when the user adds an element to the model, and constructs the above `struct`.
- **[`Muscade.doflist`](@ref)** specifies the degrees of freedom (dofs) of the element.
- **[`Muscade.residual`](@ref)** (either this of [`Muscade.lagrangian`](@ref)) takes element dofs as input and returns the element's additive contribution to the residual of a non-linear system of equations,
- **[`Muscade.lagrangian`](@ref)** (either this of [`Muscade.residual`](@ref)) takes element dofs as input and returns the element's additive contribution to a target function,
- **[`Muscade.draw`](@ref)** (optional) which draws all the elements of the same element type.

Each element must implement *either* [`Muscade.lagrangian`](@ref) *or* [`Muscade.residual`](@ref), depending on what is more natural: a beam element will implement [`Muscade.residual`](@ref) (element reaction forces as a function of nodal displacements), while an element representing a strain sensor will implement [`Muscade.lagrangian`](@ref) (log-of the probability density of the strain, given an uncertain measurement).

## No internal variables

In classical finite element formulations, plastic strain is implemented by letting an element have plastic strain and a hardening parameter as an *internal variable* at each quadrature point of the element.  Internal variables are not degrees of freedom.  Instead, they are a memory of the converged state of the element at the previous load or time step, used to affect the residual computed at the present time step.  Internal variables are also used when modeling friction and damage processes.

`Muscade` does not allow elements to have internal variables. The reason is that problem involing dofs of class `U` are not causal: our estimation of the state of a system a step `i` also depends on on measurements taken at steps `j` with `j>i`.  Hence "sweep" procedures, that is, procedures that solve a problem one load or time step at a time, are not applicable to such problems.  Solvers must hence solve for all dofs and steps at once.  When using Newton-Raphson iterations to solve problems of this class, internal variables make the Hessian matrix full: the value of a stress at a given stress depends (through the internal variable) on the strain at all preceeding steps. Or more formaly: Internal variables transforms the problem from *differential* to *integral*. A full matrix quickly leads to impossibly heavy computations.

To model phenomena usualy treated using internal variable, it is necessary in `Muscade` to make the "internal" variable into a degree of freedom, and describe the equation of evolution of this degree of freedom. See [`examples/DecayAnalysis.jl`](DecayAnalysis.md) for an example of implementation.  
!!! warning
    Because the equations of evolution involves first order time derivative, one can not use a static solver in combination with such elements.

## DataType

For a new element type `MyELement`, the datatype is defined as

```julia
struct MyElement <: AbstractElement
    ...
end
```

`MyElement` *must* be declared a subtype of [`AbstractElement`](@ref).

## Constructor

The element must provide a constructor of the form

```julia
function MyElement(nod::Vector{Node};kwargs...)
    ...
    return eleobj
end
```

which will then call the default constructor provided by Julia.  

`nod` can be used to access the coordinates of the nodes directly:

```julia
x = nod[inod].coord[icoord]
```

where `inod` is the element-node number and `icoord` the index into a vector of coordinates. `coord` is 
provided by the user when adding a `Node` to the `Model`. `Muscade` has no opinion about, and provides no
check of, the length of `coord` provided by the user. In this way elements can define what coordinate system
(how many coordinates, and their interpretation) is to be used.  Coordinate systems can even differ from one
node to the next. See also the helper function [`coord`](@ref) to get all node coordinates.

`kwargs...` is any number of named arguments, typicaly defining the material properties of the element.

The user does not call the above-defined constructor directly.  Instead, an element is added to the model by
a call of the form

```julia
e1 = addelement!(model,MyElement,nodid,kwargs...)
```

See [`addelement!`](@ref).

## Method for `Muscade.doflist`

The element must provide a method of the form

```julia
function Muscade.doflist(::Type{MyElement})
    return (inod =(...),
            class=(...),
            field=(...))
end
```

The syntax `::Type{MyElement}` is because `Muscade.doflist` will be called by `Muscade` with *a `DataType`* 
(the type `MyElement`), not with an object *of type* `MyELement` . The function name must begin
with `Muscade.` to make it possible to overload a function defined in the module `Muscade`. 

The return value of the function is a `NamedTuple` with the fields `inod`, `class` and `field`. 

- `inod`  is a `NTuple` of `Int64`: for each dof, its element-node number. 
- `class` is a `NTuple` of `Symbol`: for each dof, its class (must be `:X`, `:U` or `:A`).
- `field` is a `NTuple` of `Symbol`: for each dof, its field.

Importantly, `Muscade.doflist` does not mention dofs of class `:Λ`: if the element implements `Muscade.lagrangian`, there is automaticaly a one-to-one correspondance between Λ-dofs and X-dofs.

For example (using Julia's syntax for one-liner functions):

```julia
Muscade.doflist( ::Type{Turbine}) = (inod =(1   ,1   ,2        ,2        ),
                                     class=(:X  ,:X  ,:A       ,:A       ),
                                     field=(:tx1,:tx2,:Δseadrag,:Δskydrag))
```

See [`Muscade.doflist`](@ref).

## `Muscade.lagrangian` or `Muscade.residual`

An element must implement at least one of `Muscade.lagrangian` or `Muscade.residual`.

### Method for `Muscade.lagrangian`

Elements that implement a contribution to a target function must implement `Muscade.lagrangian`.

```julia
@espy function Muscade.lagrangian(o::MyElement,Λ,X,U,A,t,SP,dbg) 
    ...
    return L,noFB
end
```

See [`Muscade.lagrangian`](@ref) for the list of arguments and outputs.

### Method for `Muscade.residual`

Elements that implement "physics" will typicaly implement `Muscade.residual` (they could implement the same using `lagrangian`, but the resulting code would be less performant).

The interface is mostly the same as for `Muscade.lagrangian` with the differences that

- `Muscade.residual` returns a vector `R`
- there is no argument `Λ`

```julia
@espy function Muscade.residual(o::MyElement,X,U,A,t,SP,dbg) 
    ...
    return R,noFB
end
```

See [`Muscade.residual`](@ref) for the list of arguments and outputs.

### Automatic differentiation

The gradients and Hessians of `R` or `L` do not need to be implemented, because `Muscade` uses automatic 
differentiation. Because of this, it is important not to over-specify the inputs.  For example, 
implementing a function header with

```julia
@espy function Muscade.lagrangian(o::MyElement,Λ::Vector{Float64},X,U,A,t,SP,dbg)
#                                                |___bad_idea___|
```

would cause a `MethodError`, because `Muscade` will attempt to call with a `SVector` instead of `Vector`, and a special
datatype supporting automatic differentiation instead of `Float64`.

### Extraction of element-results

The function definitions of `Muscade.lagrangian` and `Muscade.residual` must be anotated with the macro call `@espy`.  
Variables within the body of `Muscade.lagrangian` and `Muscade.residual`, which
the user may want to obtain must be anotated with `☼` (by typing `\sun` they pressing `TAB`) at the place where they
are calculated. An example would be 

```julia
    ☼σ = E*ε
```

The macro will generate two versions of `Muscade.lagrangian` and/or `Muscade.residual`.  One in which the anotations 
`☼` are taken away, which is used to
solve the numerical problem.  Another with additional input and output variables, and code inserted into the body
of the function to extract results wanted by the user.

See [`Muscade.@espy`](@ref) for a complete guide on code anotations. 

### Immutables and Gauss quadrature

`Muscade.residual` and `Muscade.lagrangian` must be written in a specific style in order maximize performance
and to facilitate automatic differentiation. 

For performance, no allocation on the heap must occur.  This implies in particular that no `Array`s, and only 
(stack allocated) `StaticArray`s must be used.  For example. the code `a=zeros(n)`, creates an `Array` 
(allocated on the heap), and should be replaced with `a = SVector{N}(0. for i=1:N)` where `N` must be known 
at compile time to ensure type stability.

To facilitate automatic differentiation, no mutation must occur. `StaticArray`s are anyway not mutable.

One difficulty arises with Gauss quadrature.  Typical implementations would rely on setting `R` to zero, then adding the 
contributions from quadrature points to `R` within a `for` loop over the Gauss points, which is a mutation.  
The pseudo code:

```julia
R .= 0
for igp = 1:ngp
    F = ...
    Σ = ...
    R += F ∘ Σ ∘ ∇N * dV
end
```

shows that `R` is mutated. For the extraction of results from a loop, the following patterm **must** be used:

```julia
@espy function Muscade.residual(x,χ)
    t = ntuple(ngp) do igp
        ☼F = ...
        ☼Σ = ...
        r = F ∘ Σ ∘ ∇N * dV
        @named(r)
    end
    R = sum(   igp->t[igp].r,ngp)
    return R,...
end
```

`r` are the contributions to `R` at each quadrature point.  The operation `t = ntuple ...` returns a datastructure `t` such that  ``t[igp].r` are the contribution to the residual from the `igp`-th quadrature point. This is because

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

where the value of `expr` is that of its last line `@named(r,a,b,c)` which is a macro provided by `Muscade` that inserts the code `(r=r,a=a,b=b,c=c)`.

The code 

```julia
    a = ntuple(igp->t[igp].a,ngp)
    R = sum(   igp->t[igp].r,ngp)
```

gathers the hypothetic `a` of all quadrature points into a `Tuple` and adds together the contributions `r` into the residual `R`.  Variables behaving like `a` *might* come into play if solver feedback is provided from each Gauss point.

The macro `@named` is peculiar in that neither Julia nor `Muscade` defines a `macro named`.  Instead, it is a syntactic token identified and transformed by `@espy`.  This has one important implication: if the loop over the Gauss points only accumulates `R` (or `L`, in `Muscade.lagrangian`), it would be tempting to use a simpler pattern:

```julia
@espy function Muscade.residual(x,χ)
    R = sum(1:ngp) do igp
        ☼F = ...
        ☼Σ = ...
        F ∘ Σ ∘ ∇N * dV
    end
    return R,...
end
```

However, result extraction from inside the loop will not work: `@espy` only supports "`ntuple...do...@named` for the purpose.

See [`Muscade.residual`](@ref), [`Muscade.lagrangian`](@ref).

### Performance

For a given element formulation, the performance of `Muscade.residual` and `Muscade.lagrangian` can vary with a factor up to 100 between a good and a bad implementation.

**Type stable code** allows the compiler to know the type of every variable in a function given the type of its parameters. Code that is type unstable is significantly slower. See the page on [type stability](TypeStable.md).

**Allocation**, and the corresponding deallocation of memory *on the heap* takes time. By contrast, allocation and deallocation *on the stack* is fast.  In Julia, only immutable variables can be allocated on the stack. See the page on [memory management](Memory.md)

**Automatic differentiation** generaly does not affect how `Muscade.residual` and `Muscade.lagrangian` are written.  There are two performance-related exceptions to this:

1. If a complicated sub-function in `Muscade.residual` and `Muscade.lagrangian` (typicaly a material model or other closure) operates on an array (for example, the strain) that is smaller than the number of degrees of freedom of the system, computing time can be saved by computing the derivative of the output (in the example, the stress) with respect to the input to the subfunction, and then compose the derivatives.
2. Iterative precedures are sometimes used within `Muscade.residual` and `Muscade.lagrangian`, a typical example being in plastic material formulations.  There is no need to propagate automatic differentiation through all the iterations - doing so with the result of the iteration provides the same result.
3. Elements with corotated reference system (e.g. [beam elements](StaticBeamAnalysis.md)) can use automatic differentiation to transform the residual back to the global reference system.

See the page on [automatic differentiation](Adiff.md).

## Method for `Muscade.draw`

### Vectorization

Elements *can* implement a [`Muscade.draw`](@ref) method. If no method is implemented, the element will be invisible if the user requests a drawing of the element.

While the API may remind that of [`Muscade.lagrangian`](@ref), there is one significant difference: 
because it is more efficient to create few graphical object (few calls to `lines!`, `scatter!`) etc., the element's method for `draw` 
will be called once to draw several elements of the same type. Multiple lines can be drawn in one call to `lines!` by using `NaN`s to "lift the pen".

### Keyword arguments

When requesting a drawing of all or part of the model, the user can provide specifications (line thickness, line colors, what quantity to visualise as colored patches and so forth).  
The user can for example require

```julia
draw(model;linewidth=2)
```

The element's `draw` method *must* accept an arbitrary list of keyword arguments.  Keywords arguments not used by the method are automaticaly ignored.  In order not to fail if a *used* keyword argument is not provided by the user, the following syntax can be used in the element's `draw` method.   

```julia
function Muscade.draw(...)
    ...
    linewith = default{:linewidth}(kwargs,2.)
    ...
end
```

which can be read: if `kwargs.linewidth` exists, the set `linewidth` to its value, otherwise, set it to `2.`.

The user has facilities to draw only selected element types or selected elements, so the element's `draw` method does not need to implement a switch on *wether* to draw.

See [`examples/BeamElements.jl`](StaticBeamAnalysis.md) for an example of implementation.

### Getting element results

In many cases, drawing provides a graphical representation of element-results (see [Extracting results](@ref)).  A pattern is that `draw` creates a request and calls `residual` or `lagrangian` (which ever the element implements), with an additional last input argument `req` (the request created using [`@request`](@ref)), and an additional last output argument `out` (containing the element-results).

## Help functions

`Muscade` provides functions and constants to make it easier to comply with the API:

Element constructors can use function [`coord`](@ref) to extract the coordinates fron the `Vector{Node}` they get as first argument.

`Muscade.residual` and `Muscade.lagrangian` **must** use [`∂0`](@ref), [`∂1`](@ref) and [`∂2`](@ref) when extracting the zeroth, first and second time derivatives from arguments `X` and `U`.

Constant [`noFB`](@ref) (which have value `nothing`) can be used by elements that do not have feedback to the solving procedure.

It is sometimes convenient to handle time derivatives using automatic differentiation: elements with corotated reference systems can thus handle a moving corotated system, and thus centripetal and Coriolis forces.  See [`examples/BeamElements.jl`](StaticBeamAnalysis.md) for an example. Helper functions [`Muscade.motion`](@ref), [`Muscade.position`](@ref), [`Muscade.velocity`](@ref) and [`Muscade.acceleration`](@ref) are provided. These helper functions are not exported by `Muscade`, so their invocation must be qualified with `Muscade.`.

See [`Muscade.test_static_element`](@ref) (not exported) to compute the gradient of a lagrangian. 

For those prefering to think in terms of Cartesian tensor algebra, rather than matrix algebra, operators [`⊗`](@ref), [`∘₁`](@ref) and [`∘₂`](@ref) provide the exterior product, the single dot product and the double dot product respectively.

Elements with a corotated reference system, can make use of [`examples/Rotations.jl`](StaticBeamAnalysis.md) that provides functionality to handle rotations in ℝ³.  See [`examples/BeamElements.jl`](StaticBeamAnalysis.md) for an example.



