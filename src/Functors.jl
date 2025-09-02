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


struct Functor{name,Ta} 
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

Defines `f::Functor{:f}`, which can be passed as an argument to some constructors.
This has advantages over a closure defined as

    f(x::Real)=a*x^e

- A closure captures a variable "by reference", while `@functor` captures it by value, which
might be more intuitive. 
- To ensure type stability, the variables captured by a closure
would have to be declared `const`.
- If the code of the function is not changed, the function is not parsed and compiled again.

It is not possible to associate multiple methods to a functor.

"""
macro functor(capturedargs,foo)
    @assert capturedargs.head == :tuple
    caparglist         = capturedargs.args
    ncaparg            = length(caparglist)
    capargnames        = Vector{Symbol}(undef,ncaparg)
    for iarg           = 1:ncaparg 
        arg            = caparglist[iarg]
        capargnames[iarg] = arg isa Symbol ? arg : arg.args[1]
    end

    @assert foo.head == :(=) || foo.head == :function
    functionheader     = foo.args[1]
    functionname       = functionheader.args[1]
    functionarglist    = functionheader.args[2]
    functionsym        = QuoteNode(functionname)
    functionbody       = foo.args[2]
    functionbody       = MacroTools.postwalk(functionbody) do ex
        ex isa Symbol && ex∈capargnames ? :(o.captured.$ex) : ex
    end
    ex                 = MacroTools.postwalk(rmlines,foo)
    ex                 = MacroTools.postwalk(unblock,ex)
    qex                = QuoteNode(ex) # unannotated code for the function 
    tag                = Symbol("tag_for_the_functor_macro_",functionname)
    return prettify(esc(quote
        $functionname = Functor{$functionsym}(;$(capturedargs)...)            # f = Functor{:f}(;(a,e=2)...)  
        if  ~@isdefined($tag) || $tag≠$qex
            $tag = $qex
            function (o::Functor{$functionsym})($functionarglist)     # function (o::Functor{:f})(x::Float64)
                $functionbody                                         #     o.captured.a*x^o.captured.e 
            end                                                       # end
        end                                                          
    end))
end

