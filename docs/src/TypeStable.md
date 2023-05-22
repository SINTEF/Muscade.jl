# Type-stability

## Introduction

This text presents *type stability*, which is one of the important concepts that one needs to understand in order to write high-performance Julia code.

This text is aimed at Julia users that are familiar with composite types, abstract types, functions, methods and multiple dispatch. At the same time, as little advanced Julia syntax as possible is used, to make the text accessible.

## To type, or not to type

The developers of Julia wanted to solve the two-language problem.  They have achieved this and produced a language that "walks like Python and runs like C".  Julia "walks like Python", because it is not necessary to systematically define the type of every variable that appears in the code.  It "runs like C" because it is a compiled language, and produces (or rather, *can* produce) highly efficient machine code.

Python and MATLAB are examples of interpreted language.  In a pure interpreted language, the type of the variables is computed at run time, at the same time as the value of the variables.  As long as the values of the inputs to the code are known at the top level (in the REPL or the top script), the interpretation infers, step by step the type of the variables, all the way down the call stack. This allows to write functions without specifying types, and this in turn allows to write generic code (for example an iterative solver that works just as well with `Float64` and `Float32` variables).  The disadvantage is that inferring the type of variables on the fly introduces significant overhead at run time.

At the other end of the scale C and Fortran are examples of strictly typed compiled languages. Because the source code specifies the type of every variable in the function (both variables in the function interface, and local variables), the compiler can create efficient machine code for each function, just by considering the code of that function alone. The disadvantage is that type declaration takes time to write and clutters the source code, and (unless the language offers "templates", as C++ does), it may be necessary to write several methods, identical in all but types of variables, to make an algorithm available to various data types.

## Julia's approach to type specifications

Julia takes the sweet spot in between, not requiring to specify the type of each variable, yet producing fast machine code. The trick is as follows: every time a method is called (so, at run time), with a combination of *concrete types* of arguments that has not yet been encountered for this method, the compiler quicks in.  A "concrete type" is the information returned by `typeof()` when called on a variable.  One example is `Float64`.  This is as opposed to an abstract type, like `Real`, which is a set of concrete types, and includes `Float64` and `Float32`. *In the rest of this text "type" will refer to "concrete type"*.

The compiler now has the source code of the method, and the types of all the arguments. The compiler will produce a *method instance* (or instance, for short), which is machine code for this combination. One interesting implication is that writing strictly typed method interfaces in Julia does not provide any improvement of machine code performance: the compiler takes the type of the arguments from the calling context anyway. A strictly typed interface has the disadvantage of offering no flexibility. A method that only accepts a `Vector` will not accept other vector-like things like a `SubArray` (an array view), a `Adjoint` (a transposed array), a `SparseMatrix` or a `StaticArray`, even thought the method probably implements an algorithm that would compile perfectly well for all of these.

However, providing partial specification of the type of the arguments of a method serves important purposes in Julia:

1. If a function has several methods, it allows to specify which method should be executed (multiple dispatch). This is where abstract types like `Real`,  `AbstractVector` and `AbstractVector{<:Real}` come into their own.
2. It improves code readability, stating for example "this method expects some vector of some real numbers - but not a string".
3. It provides more graceful failures: "function `foo` has no method that takes in a string" is more informative that some esoteric failure down the line when attempting to add two strings.

## What is type stability?

If the source code of the method is well written, the source code and the concrete type of all arguments is enough information for the compiler to infer the concrete type of every variable and expression within the method.  The method is then said to be "typestable", and the Julia compiler will produce efficient code.

If, for a variety of reasons that will be studied in the following, the type of a local variable cannot be inferred from the types of the arguments, the compiler will produce machine code full of "if"s, covering all options of what the type of each variable could be. The loss in performance is often significant, easily by a factor of 10.

**If you are yourself able to infer the type of every local variable, and every expression in a method (or script) from the types (*not* the values) of the arguments or from constants in the code, the function will be typestable.**  Actually, as will be seen below, this inference of types is also allowed access to `struct` declarations, and to the types of the return values of functions called by the function you are studying.

The rest of this text will examine a variety of situations, ranging from obvious to more tricky tricky, in which it is not possible to infer the types of local variables from the types of the arguments, resulting in type instability.

For this purpose, it will be useful to write down the information available to the compiler.  So for example, if the method

```julia
function add(a::Number,b::Number)
    c = a+b
    return c
end
```

is called with `a` of type `Float64` and `b` of type `Int32`, then we will write the information available to the compiler to create an instance as

```julia
instance add(a::Float64,b::Int32)
    c = a+b
    return c
end
```
`instance` is not Julia syntax, it is just a notation introduced in this text to describe an instance.  In such `instance` description, a *concrete* type must be associated with every argument.

## If, then

Consider the following method instance

```julia
instance largest(a::Float64,b::Int64)
    if a > b
        c = a
    else
        c = b
    end
    return c
end
```

The variable `c` will be set to either `a` or `b`. `c` will take the value *and the type* of either *a* or *b*.  The type of `c` depends on an operation `a > b` on the *values* of `a` and `b`: the type of `c` cannot be inferred from the type of arguments alone, and this code is not typestable.

Several approaches might be relevant to prevent type instability.  The simplest is to code `largest` so that it only accepts two arguments of the same type.

```julia
function largest(a::R,b::R) where{R<:Real}
    if a > b
        c = a
    else
        c = b
    end
    return c
end
```

The method is general, it can result in the generation of method instances like `instance largest(a::Float64,b::Float64)`, `instance largest(a::Int64,b::Int64)` and many others. It cannot result in the generation of machine code for `instance largest(a::Float64,b::Int64)` (because `R` cannot be both `Int64` and `Float64`). If we need to be able to handle variables of different types, yet want type stability, a solution is to use *promotion* to ensure that `c` is always of the same type.

```julia
function largest(a,b)
    pa,pb = promote(a,b)
    if a > b
        c = pa
    else
        c = pb
    end
    return c
end
```

`promote` is defined so that `pa` and `pb` have the same type, and this type is inferred from the types of `a` and `b`. For example, for a call `instance largest(a::Float64,b::Int64)`, the types of `pa`, `pb` and `c` will be `Float64`, to which one can convert a `Int64` variable without loss of information (well, mostly).

**Do not allow an if-then construct to return a variable which type depends on the branch taken.**

## Method return value

A method `foo` that would call the above first, not typestable, version of the method instance `largest` would receive as output a variable of a type that is value dependent: `foo` itself would not be typestable.  The workaround here is to create typestable methods for `largest`, as suggested above.

One example is the method `Base.findfirst(A)`, which given a `Vector{Boolean}` returns the index of the first `true` element of the vector.  The catch is that if all the vector's elements are `false`, the method returns `nothing`. `nothing` is of type `Nothing`, while the index is of type `Int64`.  Using this method will make the calling method not typestable.

**Avoid methods that return variables of value-dependant types.**

## Array of abstract element type

Consider the following code

```julia
v = [3.,1,"Hello world!"]
function showall(v)
    for e ∈ v
        @show e
    end
end
showall(v)
```

The above call `showall(v)` generates a method instance

```julia
instance showall(v::Array{Any,1})
    for e ∈ v
        @show e
    end
end
```

The concrete type of `e` cannot be inferred from `Array{Any,1}`, because `Any` is not a concrete type. More specifically, the type of `e` changes from one iteration to the next: the code is not typestable. If `v` is of type `Array{Any,1}`, even if `V` has elements that are all of the same type, this does not help:

```julia
v = Vector{Any}(undef,3)
v[1] = 3.
v[2] = 1.
v[3] = 3.14
showall(v)
```

`e` may have the same type at each iteration, but this type still cannot be inferred from the type `Array{Any,1}` of the argument.

If we define `w = randn(3)`, `w` has type `Array{Float64,1}`.  This is much more informative: every element of `w` is known to have the same concrete type `Float64`. Hence the call `showall(w)` generates a method instance

```julia
instance showall(v::Array{Float64,1})
    for e ∈ v
        @show e
    end
end
```

and the compiler can infer that `e` is a `Float64`.

**Wherever possible use arrays with a concrete element type.**

Sometimes, the use of array with abstract element type is deliberate.  One may really wish to iterate over a heterogeneous collection of elements and apply various methods of the same function to them: we design for dynamic dispatch, and must accept that the process of deciding which method to call takes time.  Two techniques can be used to limit the performance penalty.

The first is the use of a "function barrier": The loop over the heterogenous array should contain as little code as possible, ideally only the access to the arrays element, and the call to a method.

```julia
for e ∈ v
    foo(e)
end
```

If `v` contains elements of different type, the loop is not typestable and hence slow. Yet each value of `e` at each iteration has its unique concrete type, for which an instance of `foo` will be generated: `foo` can be made typestable and fast.

The second, a further improvement of the first, is to group elements by concrete type, for example, using a heterogenous arrays of homogeneous arrays.

```julia
vv = [[1.,2.,3.],[1,2]]
for v ∈ vv  # outerloop
    innerloop(v)
end
function innerloop(v)
    for e ∈ v
        foo(e)
    end
end
```

Here `vv` is an `Array{Any,1}`, containing two vectors of different types. `vv[1]` is a `Array{Float64,1}` and `vv[2]` is a `Array{Int64,1}`.
Function `innerloop` is called twice and two instances are generated

```julia
instance innerloop(v::Array{Float64,1})
    for e ∈ v  # e is Float64
        foo(e)
    end
end
instance innerloop(v::Array{Int64,1})
    for e ∈ v  # e is Int64
        foo(e)
    end
end
```

and in both instances, the type of `e` is clearly defined: the instances are typestable.

The with this second approach is that the loop `for v ∈ vv` has few iterations (if the number of types is small compared to the number of elements in each types).

## Structure of abstract field type

A similar loss of type stability arises when reading data from structures that have a field of abstract type:

```julia
struct SlowType
    a
end
struct JustAsBad
    a::Real
end
struct MuchBetter
    a::Float64
end
function show_a(s)
    @show s.a
end
show_a(SlowType(3.))
show_a(JustAsBad(3.))
show_a(MuchBetter(3.))
```

The first call to `show_a` generates

```julia
instance show_a(s::SlowType)
    @show s.a # The concrete type of field a of type SlowType cannot be
              # inferred from the definition of SlowType
end
```

The second call to `show_a` has the same problem.  The third call generates a typestable instance

```julia
instance show_a(s::Better)
    @show s.a # That's a Float64
end
```

It is often interesting to create structures with fields that can have various types. A classic example is Julia's `Complex` type, which can have real and imaginary components which are either both `Float64`, both `Float32` or other more exotic choices. This can be done without losing type stability by using parametric types:

```julia
struct FlexibleAndFast{R}
    a::R
end
show_a(FlexibleAndFast(3.))
show_a(FlexibleAndFast(3 ))
```

The above calls generate two typestable instances of `show_a`

```julia
instance show_a(s::FlexibleAndFast{Float64})
    @show s.a # That's a Float64
end
instance show_a(s::FlexibleAndFast{Int64})
    @show s.a # That's an Int64
end
```

**Always use `struct` with fields of concrete types.  Use parametric structure where necessary**.

## A note on constructors for parametric types

Consider a `struct` definition without inner constructor:

```julia
struct MyType{A,B}
    a::A
    b::B
end
```

Julia will automatically generate a constructor *method* with signature

```julia
MyType{A,B}(a::A,b::B)
```

Julia will also produce another method with signature

```julia
MyType(a::A,b::B)
```

because for `MyType`, it is possible to infer all type parameters from the types of the inputs to the constructor. Other constructors like

```julia
MyType{A}(a::A,b::B)
```

have to be defined explicitly (how should the compiler decide whether to interpret a single type-parameter input as `A` or `B`...).

Consider another example:

```julia
struct MyType{A,B,C}
    a::A
    b::B
end
```

Julia will automatically generate a constructor method with signature

```julia
MyType{A,B,C}(a::A,b::B)
```

but will not generate other methods.  A method like

```julia
MyType{C}(a::A,b::B)
```

would have to be defined explicitly.

## StaticArrays

Julia `Array`s are an example of parametric type, where the parameters are the type of elements, and the dimension (the number of indices). Importantly, the *size* of the array is not part of the *type*, it is a part of the *value* of the array.

The package `StaticArrays.jl` provides the type `StaticArray`, useful for avoiding another performance problem: garbage collection that follows the allocation of `Array`s on the heap. This is because `StaticArray` are allocated on the stack, simplifying runtime memory management.

```julia
using StaticArrays
SA = SVector{3,Float64}([1.,2.,3.])
SA = SVector(1.,2.,3.)
SA = SVector([1.,2.,3.])
```

The first call to `SVector` is typestable: all the information needed to infer the type of `SA` is provided in curly braces. The second call is typestable too, because the compiler can deduce the same information from the type and number of inputs. The third call is problematic: while the type of the elements of `SA` can be inferred by the compiler, the length of `[1.,2.,3.]` is part of this array's value, not type. The type of `SA` has a parameter that depends on the value (the size) of the argument passed to the constructor.  Not only does this generate an instance of the constructor that is not type stable, but the non-inferable type of `SA` "contaminates" the calling code with type instability.

## Val

What if we want to write a function that takes an `Vector` as an input, processes it (for example just keeps it as it is), and returns a `SVector` of the same shape. Of course we want this function to be general and not be limited to a given array size *and* we want this function to be typestable, for good performance.

First attempt:

```julia
function static(v::Vector)
    return SVector{length(v),eltype(v)}(v)
end
```

This function is not typestable. It constructs a variable of type `StaticArray{(3,),Float64}`, where `3` is obtained as the `length` of `v`, and the length is part of the value of an `Array`.  Value-to-type alarm!

One possible solution is to use `Val`. Let us say that `static` is called by a function `foo` within which the length of `v` can be inferred at compile time.  We could create the following code

```julia
function static(v,::Val{L}) where{L}
    return SVector{L,Float64}(v)
end
function foo()
    Val3 = Val(3)
    Val4 = Val(4)
    @show static([1.,2.,3.]   ,Val3)
    @show static([1.,2.,3.,4.],Val4)
end
```

The call `Val(3)` generates a variable, of type `Val{3}`. Clearly, `Val` as a function is not typestable, since it creates a variable of a type depending on the value of its argument.

However, function `foo` is typestable.  This may come as a surprise, but two things conspire to allow this:

1. The source code of `foo` explicitly mentions the *constants* `3` and `4`, and the compiler has access to it.
2. The compiler is greedy - it evaluates at compile time whenever possible.  Hence the call `Val(3)` is evaluated during compilation, and `Val3` is known to the compiler to be a a value-empty variable of type `Val{3}`.

In `foo`, the method `static` is called twice, leading to the generation of two typestable instances

```julia
instance static(v,::Val{3})
    return SVector{3,Float64}(v)
end
instance static(v,::Val{4})
    return SVector{4,Float64}(v)
end
```

What if the length of the vectors is not defined as a constant in `foo`?  If this length is the result of some computation, the call to `Val` with not be typestable. If `foo` is high enough in the call hierarchy, and outside any time-critical loop, this is not an issue: only `foo` will not be typestable, but functions that it calls can still be typestable (cf. the function barrier pattern).

**`Val` allows to move type instability up the call hierarchy, or eliminate it altogether.**

## Functions

The type `Function` is an *abstract* datatype, and every function in Julia has its own type.  Here we refer not to the type of the variables returned by the function, but to the function being a variable in itself.

The implication is that if we have a scalar-valued function `energy` that takes a function `signal` as an input, and computes the energy of the signal over some interval, then a new instance of `energy` will be compiled every time it is called with an new argument `signal`.

This also has implications on how to store functions in a `struct`. This is not typestable

```julia
struct MyType
    foo::Function
end
```

but this is

```julia
struct MyType{Tfoo}
    foo::Tfoo
end
```

## @code_warntype

One important tool to check that an instance is typestable is the macro `@code_warntype`. For example

```julia
v = randn(3)
@code_warntype Val(length(v))
@code_warntype static(v,Val(length(v)))
```

The first invocation of `@code_warntype` outputs a semi-compiled code, and highlights some of the types in red: the call `Val(3)` is not typestable. The second invocation of `@code_warntype` produces an output in which all types are highlighted in blue: the call to `static` is typestable.  Note that `@code_warntype` only analyses the compilation of the outermost function `static` - given the arguments `v` and `Val(length(v))`.

## Profile.jl

`Profile.jl` and `ProfileView.jl` together provide a "flame graph", a graphical representation of where processor time goes, in which code that is not typestable is highlighted. Output from the profiler often shows how type instability propagates: a single variable that is not typestable makes "anything it touches" type unstable.

Particularly useful, one can click on a function, and then type `warntype_last()` in the REPL to get to see a `@code_warntype` output for that function.