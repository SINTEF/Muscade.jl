# Memory management

## Introduction

Mutable structs, arguments passes by value to functions, garbage collection, variable names as tags. A variety of concept in Julia become easier to understand by considering how Julia manages memory, and how this affect performance.

In line with programming languages of the past decades, Julia uses two approaches to store variables, namely the *heap*, and the *stack*. Julia's designer made a few important choices on how they use heap and stack, and these strongly shape how we use the language, in particular when writing high-performance code.

## The heap

The heap is a swath of memory made available to a program.  When a variable is *allocated*, a preferably contiguous amount of available memory of the right size is found on the heap, and its adress (a pointer) is returned.  Values within that segment of memory can be written to or read, and the variable is *deallocated* when it is no longer in use.  

This is very flexible, but this flexibilty comes at a cost: To allocate a variable one must find available memory, thus there must be a heap-ledger describing what memory is available. Browsing the heap-ledger takes time.

Early in the execution of a program, the heap is unused, and finding space for new variables is easy.  But as variables are allocated and deallocated in arbitrary order, large contiguous area of free memory become rarer: the heap is fragmented, and there might be a need to swap things around (defragmentation).

Further, language designers must make a choice.  Option one is to let the programmer explicitely deallocate a variable when it is no longer needed.  Unfortunately, bugs in which a variable is not deallocate easily occur. If this happens in some loop, memory is allocated but not deallocated (a memory leak) and these bugs are nasty to track down.

The safer option 2, adopted by Julia, is to automaticaly deallocate memory on the heap when nothing anymore points to it.  That implies that there must be a pointer-ledger of all pointers into the heap, an the language must periodicaly go through the pointer-ledger to find orphaned heap-memory, update the heap-ledger, and possibly defragment the heap (move the variables in order to create large contiguous unallocated memory).  This process is called garbage collection.  While invisible to the programmer, the user sees it: it takes time.

When a profiler reports the number of allocations, it actualy refers to allocation on the heap, not counting variables created on the stack.

## The stack

The stack is a limited amount of memory managed on a "last in first out" basis: as a picture, the stack is vertical and when the stack start empty, the stack pointer points tot he base of the stack. If a new variable needs to be stored, it is added to the top of the stack and the pointer is updated. When the variable no longer needs to be stored, the pointer is just updated back to its previous value (pop the stack).

The mechanism for creating or destroying a variable is simple, and thus extremely fast.  The drawback of course is that variables are destroyed in reverse order of their creation, which is not always convenient.  However this works perfectly for: 

1. input arguments to a function, in languages where these are passed by values, and
2. local variables for the function.

All these variables are created when entering the function, and destroyed when leaving it.  The (memory) stack thus fills as the code goes deeper into the call stack.  Recursions gone beserk typicaly result in a (memory) stack overflow.

That, at least, is how things looked like in the late pre-internet age.  Since then, CPUs have been equipped with a nested system of fast access caches. An outer cache, smaller and fast than RAM, a inner cache, yet smaller and faster, and the CPU's registers are, in a sense, the innermost caches.  It makes sense to store the data currently being processed into one of these inner caches.  But this comes at the cost of complexity to the language and compiler designers: *moving* data to the cache means to free its original location (which is awkward in the context of a stack). *Copying* data introduces the classic problem of keeping all copies of the data up to date.

Thus the designers of the Julia language made one important choice: *once a variable on the stack has been assigned a value, that value can never be changed*.  The variable is said to be *immutable*. This allows Julia copy any part of the stack to the CPU caches to optimize performance, without worrying about out-of-date copies. "Immutables are easier to reason about" says the doc: this may not be true when learning Julia, but applies when writing the compiler. 

Now we see why variable names in Julia are best viualised as tags on a value (you copy the tag if you copy the value, and there can be multiple tags for the same value) as opposed to the parable of a variable as a box in which you can store different value (valid in other languages).

Julia's designers made another important choice: the size of all variables stored the stack must be known at compile time. This drasticaly simplifies the process of creating or destroying space on the stack when entering a function: just increment or decrement to stack pointer by a value determined for the function (actualy the method instance) at compile time.  Accessing a local variable from inside the function is likewise very fast: add a compile-time constant to the stack pointer, and that's were your data is.

As we will see, the fact that variables must be 

1. Immutable
2. Of size known at compile time

significantly affects how one programs "on the stack" in Julia.

## An `Array`

In Julia, an `Array{2,Float64}` is "copied", for example by passing it as an argument to a function. What happens?

The array comes in several parts

1. Heap-memory enough to store all the values in the array, in column major order. 
2. A smaller amount of heap-memory containing.
    a.  the sizes of the first and second indices of the array.
    b.  a pointer to the storage 1. of the values.
3. On the stack, a pointer to 2.
4. Machine code that is written under the assumption that this array has two indices, that the size of an array element is 64 bits, that array elements can be copied directly to the CPU register for algebraic operations, etc.

When the array is passed to a function, only 3., the pointer is actually passed to the function, as a copy placed at the right spot on the stack.  Being on the stack, the pointer is immutable: the function cannot reallocate 2. (as then the pointer to 2. would change value, pointing to a new spot on th heap). But one can change the content of the array, and its sizes (`push!`). Upon return from the function, the caller still has the pointer to a place on same spot on the heap: any change made by the function to values inside the array is visible by the caller.

If function overwrites the array as a whole however, the interpretation is different

```julia
function foo(a)
    a = [0,0]
end
A = [1,2]
foo(A)
@show A

2-element Vector{Int64}:
 1
 2
```

Here `foo` was given a pointer to an array. By writing `a =`, `foo` discards its knowledge of the pointer, and uses the same tag-name to refer to a new spot of memory on the heap containing a pair of zeros.

## Strategies for performance

Allocating and deallocating memory on the heap takes time.  The amount of time for each allocation is tiny, but an ocean is just many drops.
The time required is not proportional to the amount of memory allocated, but to the amount of individual allocations made.  So let us consider as an example a function inside of a hot loop, that needs internal memory for its local computations, and to return the results it produces.  For the function to be fast it must not allocate/deallocate memory on the heap. Two different strategies can be used to this effect, which we could call "procedural" and "functional".

## Procedural strategy

In the procedural strategy, memory (say `Array`s) is preallocated "once and for all" on the heap, and passed to the function: arrays or slices of arrays are passed to the function.  The function uses this as work arrays and/or to return outputs: it is said to work "in place".  In that way there is one (or several) allocation[s] before the hot loop, and many calls to the function.  In Julia, the naming convention for functions that thus modify their input arguments is to end the function name with `!`.

```julia
function double!(a)
    for i∈eachindex(a)
        a[i] *= 2
    end
end

a = randn(10)
for i = 1:10   # "hot" loop
    double!(a)
end
```

This strategy is sometimes difficult to implement:

1. If the function requires some working space, passing arrays to it is a breach of separation of concern (a new algorithm might require less, or more memory).
2. The return type of some functions might be hard to predict: even for type-stable functions, evaluating the types of an algorithm's output given the types of its inputs is sometimes best left to the compiler.
3. Algebraic operations, for example, are built into larger expression. In an implementation of `*`, to be used in `a*b+c`, operating in place is not an option (see however syntaxes like `d .= a.*b.+c` that *are* operating in place).

## Functional strategy

In the functional strategy, the function creates new variables, both for intermediate results and for return value[s].  For performance, these variables are created on the stack, and must thus be immutable and of size known at compile time. 

```julia
function double(a)
    b = 2 .*a
    return b
end

using StaticArrays
a = SVector{10}(randn() for i=1:10)
for i = 1:10   # "hot" loop
    a = double(a)
end
```

Note that `double` is called with a `SVector`.  `SVector`s' size is a parameter to the type, and is thus known at compile time. `SVector`s' are immutable (`a[2] = 0` will throw an error), and thus `SVector`s live on the stack. If `double` was called with a `Vector` (which size is part of the value, and which is mutable), the operation `b = 2 .*a` would result in an allocation.

This strategy may also be difficult to implement:

1. If many array sizes may be used (leading to the compilation of many method instances).
2. If array sizes depend on input *values*, using types in which size is part of the type leads to type-unstable code.
3. If large arrays are used, compile times will become excessively high.

## Working with immutables

`SArray`s, `ntuple`s and `NamedTuple`s are useful immutable datastructures.

If a for example an `SArray` can not be modified after it is created, how can we fill it with values?  Julia provides two syntaxs to do this, comprehensions and `do`-loops.  A comprehension looks like this:

```julia
using StaticArrays

const N=3 
f(i) = i^2

a = SVector{N}(f(i) for i=1:N)
```

Importantly, `N` must be a compile-time constant.

A `do`-loop allows more complicated operations. We assume that the function `material` returns a `StaticArray` `σ` and a variable `χgp` of type difficult to predict. The `do`-loop creates the `ntuple` `accum`, which is then unpacked (using comprehensions) to sum the `σ` values into `r` and stack the `χgp` values into `χ`.

```julia
function residual(x,y)
    ngp   = 4
    accum = ntuple(ngp) do igp
        z = x[igp] + y[igp]
        σ, χgp = material(z)
        (σ=σ, χ=χgp) 
    end # do igp
    r = sum(        accum[igp].σ for igp=1:ngp)
    χ = NTuple{ngp}(accum[igp].χ for igp=1:ngp)
    return r,χ
end 
```

