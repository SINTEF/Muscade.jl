"""
    f  = Muscade.FunctionFromVector(xs::AbstractRange,ys::AbstractVector)
    y  = f(x)

    Linear interpolation.  Fails if `x` is outside the range `xs`
"""
struct FunctionFromVector{X,Y} <:Function
    x::X  
    y::Y
end
function (f::FunctionFromVector{<:AbstractRange})(x)
    @assert first(f.x)≤x≤last(f.x)
    ix = (x-first(f.x))/step(f.x) +1
    i  = min(trunc(𝕫,ix),length(f.x)-1)
    di = ix-i
    return f.y[i]*(1-di)+f.y[i+1]*di
end
function (f::FunctionFromVector{<:AbstractVector})(x)
    @assert first(f.x)≤x≤last(f.x)
    i = findlast(i->i≤x,f.x)
    if i==length(f.x) 
        i -=1
    end
    di = (x-f.x[i])/(f.x[i+1]-f.x[i])
    return f.y[i]*(1-di)+f.y[i+1]*di
end


"""
    function MyElement(nod::Vector{Node}; ... foo::Functor ...) 

`@functor` returns an object of type `Functor`. `Functor` is a subtype of the abstract type `Function`.

In an element constructor, requiring a function-like input to be a `Functor` forces the user of the 
element to use the `@functor` macro to define functions provided as input to an element.  This is designed
to safeguard the user from some arguably unintuitive behaviour of closures.    

See also: [`@functor`](@ref)
"""
struct Functor{name,Ta} <: Function
    captured::Ta
    function Functor{name}(;kwargs...) where{name}
        nt = NamedTuple(kwargs)
        return new{name,typeof(nt)}(nt)
    end
end
"""
    a = 3
    @functor with(a,e=2) function f(x::Real)
        return a*x^e
    end

or

    a = 3
    @functor with(a,e=2)  f(x::Real)=a*x^e
    e = 1
    @functor with(a,e)    f(x::Real)=a*x^e
    @functor with()       f(x::Real)=x^2

Creates a function-like object, of type `Functor`.

This is roughly equivalent to a closure defined as

    f(x::Real)=a*x^e

Functors are meant to facilitate the definition of "functions" in a Muscade input script, 
while avoiding several of the issues associated with defining a function (and in particular
a closure) in a script:    

- A closure captures a variable "by reference", while `@functor` captures it by value, which might be more intuitive. 
- To ensure type stability, the variables captured by a closure would have to be declared `const` - forbidding to update 
  the input value in a script without restarting Julia.
- If the code of the function is not changed, the function is not parsed and compiled again, accelerating the re-analysis.

It is not possible to associate multiple methods to a functor.

See also: [`Functor`](@ref)
"""
macro functor(capture,foo)
    # to debug, use 'Base.dump' on expressions
    # Build capargname, a vector of names of captured variables, to later replace a -> o.captured.a in the body of foo -> (o::Functor{:foo})
    if capture.head==:call # function call
        if length(capture.args)==1 # no captured arguments
            caparg     = Any[]
            capargname = Symbol[]
            ncaparg    = 0
        else
            if capture.args[2] isa Expr && capture.args[2].head isa Symbol && capture.args[2].head == :parameters # user not supposed to prefix captured args with ;, but I'm in a good mood
                error("Do not use ; in list of captured arguments")
            else
                caparg     = capture.args[2:end]
            end
            ncaparg    = length(caparg)
            capargname = Vector{Symbol}(undef,ncaparg)
            for icaparg = 1:ncaparg
                if caparg[icaparg] isa Symbol
                    capargname[icaparg] = caparg[icaparg]
                elseif caparg[icaparg] isa Expr
                    if caparg[icaparg].head == :kw
                        capargname[icaparg] = caparg[icaparg].args[1]
                    elseif caparg[icaparg].head == :parameters
                        muscadeerror("Invalid @functor definition 3")
                    end
                else
                    muscadeerror("Invalid @functor definition 4")
                end
            end
        end
    else
        muscadeerror("Invalid @functor definition 2")
    end

    # Build the code for the method associated to the functor 
    # TODO all variables must be either capturedargs or fooargs, no closure. Throw error otherwise
    foodict            = splitdef(foo)
    foodict[:body]     = MacroTools.postwalk(foodict[:body]) do ex
        ex isa Symbol && ex∈capargname ? :(o.captured.$ex) : ex # prefix captured args with `o.captured.` in method body
    end    
    fooname            = foodict[:name]                           # :foo
    functortype        = Expr(:curly,:Functor,QuoteNode(fooname)) # :(Functor{$fooname})                     # Functor{:foo}

    foodict[:name]     = Expr(:(::),:o,functortype)             # (o::Functor{:foo}), name of the method that implements foo(x)
    foo                = combinedef(foodict)                    # code of said method
    quotefoo           = MacroTools.postwalk(rmlines,foo)
    quotefoo           = MacroTools.postwalk(unblock,quotefoo)
    quotefoo           = QuoteNode(quotefoo)                    # quote of unannotated code for said method, to decide wether foo-code changed or not 

    # obscure variable name, to prevent reparsing of the foo definition
    tag                = Symbol("tag_for_the_functor_macro_",fooname)

    # build the code for the call to the functor constructor
    caparg           = Expr(:parameters,caparg...)     # place a ; in front of the argument list (any prefixed ; was cleaned earlier)
    constrcall       = Expr(:call,functortype,caparg)  # Functor{:foo}(;a,b=2)
    constructfunctor = Expr(:(=),fooname,constrcall)
    
    code = esc(quote 
        $constructfunctor
        if  ~@isdefined($tag) || $tag ≠ $quotefoo
           $tag = $quotefoo
           $foo
        end                                                          
    end)
    #@show prettify(code)
    return code
end





