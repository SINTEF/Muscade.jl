"""
    f  = FunctionFromVector(xs::AbstractRange,ys::AbstractVector)
    y  = f(x)

    Linear interpolation.  Fails if `x` is outside the range `xs`
"""
struct FunctionFromVector{X,Y} <:Function
    x::X  
    y::Y
end
FunctionFromVector(x::AbstractRange,y::AbstractVector{R}) where{R<:Real} = FunctionFromVector{typeof(x),typeof(y)}(x,y)
function (f::FunctionFromVector)(x)
    ix = (x-first(f.x))/step(f.x) +1
    @assert 1≤ix≤length(f.x)

    i  = min(trunc(𝕫,ix),length(f.x)-1)
    di = ix-i
    return f.y[i]*(1-di)+f.y[i+1]*di
end

struct Functor{name,Ta} <: Function
    captured::Ta
    function Functor{name}(;kwargs...) where{name}
        nt = NamedTuple(kwargs)
        return new{name,typeof(nt)}(nt)
    end
end
"""
    a = 3
    @functor (a,e=2) function f(x::Real)
        return a*x^e
    end

or

    a = 3
    @functor (a,e=2)  f(x::Real)=a*x^e

This is roughly equivalent to a closure defined as

    f(x::Real)=a*x^e

Functors are meant to facilitate the definition of "functions" in a Muscade input script, 
while avoiding several of the issues associated with defining a function (and in particular
a closure) in a script:    

- A closure captures a variable "by reference", while `@functor` captures it by value, which
might be more intuitive. 
- To ensure type stability, the variables captured by a closure
would have to be declared `const` - forbidding to update the input value without restarting Julia
- If the code of the function is not changed, the function is not parsed and compiled again, accelerating
the re-analysis.

It is not possible to associate multiple methods to a functor.

`Functor` is a subtype of the abstract type `Function`: functions that accept a `arg::Function` as an 
input will accept `arg to be a `Functor`.  Functions that require `arg:Functor` will not accept a classical
`Function`, allowing to ensure capture by value etc.

"""
macro functor(capturedargs,foo)
    @assert capturedargs.head == :tuple
    caparglist         = capturedargs.args
    ncaparg            = caparglist == Any[:($(Expr(:parameters)))] ? 0 : length(caparglist)
    capargnames        = Vector{Symbol}(undef,ncaparg)
    for iarg           = 1:ncaparg 
        arg            = caparglist[iarg]
        capargnames[iarg] = arg isa Symbol ? arg : arg.args[1]
    end

    foodic             = splitdef(foo)
    foodic[:body]      = MacroTools.postwalk(foodic[:body]) do ex
        ex isa Symbol && ex∈capargnames ? :(o.captured.$ex) : ex # prefix captured args with `o.captured.`
    end    
    functionname       = foodic[:name]
    functionsym        = QuoteNode(functionname)
    foodic[:name]      = :((o::Functor{$functionsym}))
    foo                = combinedef(foodic)
    ex                 = MacroTools.postwalk(rmlines,foo)
    ex                 = MacroTools.postwalk(unblock,ex)
    qex                = QuoteNode(ex) # unannotated code for the function 
    tag                = Symbol("tag_for_the_functor_macro_",functionname)
    return prettify(esc(quote
        $functionname = Functor{$functionsym}(;$(capturedargs)...)    # f = Functor{:f}(;(a,e=2)...)  
        if  ~@isdefined($tag) || $tag≠$qex
            $tag = $qex
            $foo
        end                                                          
    end))
end

macro nakedfunctor(capturedargs,foo)
    @assert capturedargs.head == :tuple
    caparglist         = capturedargs.args
    ncaparg            = caparglist == Any[:($(Expr(:parameters)))] ? 0 : length(caparglist)
    capargnames        = Vector{Symbol}(undef,ncaparg)
    for iarg           = 1:ncaparg 
        arg            = caparglist[iarg]
        capargnames[iarg] = arg isa Symbol ? arg : arg.args[1]
    end

    foodic             = splitdef(foo)
    foodic[:body]      = MacroTools.postwalk(foodic[:body]) do ex
        ex isa Symbol && ex∈capargnames ? :(o.captured.$ex) : ex # prefix captured args with `o.captured.`
    end    
    functionname       = foodic[:name]
    functionsym        = QuoteNode(functionname)
    foodic[:name]      = :((o::Functor{$functionsym}))
    foo                = combinedef(foodic)
    return prettify(esc(quote
        $functionname = Functor{$functionsym}(;$(capturedargs)...)    # f = Functor{:f}(;(a,e=2)...)  
        $foo
    end))
end

