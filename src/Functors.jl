abstract type Interpolator <: Function end

struct RangeInterpolator{Nc,X,Y} <:Interpolator
    x ::X
    y ::Y
    x0::Float64
    Δx⁻¹::Float64 # average time step
end
Nout(o::RangeInterpolator{Nc}) where{Nc} = Nc
function findinterval(o::RangeInterpolator{Nc,X,Y},x) where{Nc,X<:AbstractRange,Y}
    i⁻ = floor(Int64,(x-o.x0)*o.Δx⁻¹)+1
    i⁻ = max(1            ,i⁻)
    i⁻ = min(length(o.x)-1,i⁻)
    return i⁻
end
function findinterval(o::RangeInterpolator{Nc,X,Y},x) where{Nc,X<:AbstractVector,Y}
    i⁻ = floor(Int64,(x-o.x0)*o.Δx⁻¹)+1
    i⁻ = max(1            ,i⁻)
    i⁻ = min(length(o.x)-1,i⁻)
    while i⁻>1 && x<o.x[i⁻]  
        i⁻ -=1
    end
    while i⁻<length(o.x)-1 && x>o.x[i⁻+1]  
        i⁻ +=1
    end
    return i⁻
end

struct VectorInterpolator{Nc,X,Y} <:Interpolator
    x ::X
    y ::Y
end
Nout(o::VectorInterpolator{Nc}) where{Nc} = Nc
function findinterval(o::VectorInterpolator,x)
    x≤first(o.x) && return 1
    x≥last( o.x) && return length(o.x)-1
    i⁻ = 1
    i⁺ = length(o.x)
    Δi = i⁺-i⁻
    while Δi > 1
        i = div(Δi,2)+i⁻
        if x > o.x[i]
            i⁻ = i
        else
            i⁺ = i
        end
        Δi = i⁺-i⁻
    end
    return i⁻
end


struct linearinterp{Nc} end
function linearinterp{:scalar}(o,x,i⁻)
    y⁻, y⁺ = o.y[i⁻], o.y[i⁻+1]
    x⁻, x⁺ = o.x[i⁻], o.x[i⁻+1]
    return y⁻ + (y⁺-y⁻)*((x -x⁻)/(x⁺-x⁻))
end
function linearinterp{Nc}(o,x,i⁻) where{Nc}
    x⁻, x⁺ = o.x[i⁻], o.x[i⁻+1]
    r      = (x -x⁻)/(x⁺-x⁻)
    return SVector{Nc}(o.y[j,i⁻] + (o.y[j,i⁻+1]-o.y[j,i⁻])*r  for j∈1:Nc)
end
"""
    f = interpolator(X,Y;[quasirange=false])
    y = f(x)

Linear interpolation over a 1D domain.  Used to interpolate over time series
from logging data. Extrapolates linearly.

- `X` is either an `AbstractVector` or an `AbstractRange`, of length `Nx≥2` (the times of logging), containing strictly increasing values.
- `Y` is either a `Vector` of length `Nx`, or a `Matrix` of size `(Nc,Nx)` where `Nc` is the number of channels to interpolate.
- `quasirange` (default `false`).  If `X` is an `AbstractVector`, states that the logging time only have small deviations from values in an `AbstractRange`.

- `x` must be a scalar
- `y` will be either a scalar or a `SVector{Nc}`, depending on the type of `Y`.

"""
function interpolator(x,y;quasirange=false)
    Nc = isa(y,AbstractVector) ? :scalar : size(y,1)
    if quasirange || x isa AbstractRange 
        x0   = first(x)
        Δx⁻¹ = (length(x)-1)/(last(x)-x0)
        return RangeInterpolator{Nc,typeof(x),typeof(y)}(x,y,x0,Δx⁻¹)
    else
        return VectorInterpolator{Nc,typeof(x),typeof(y)}(x,y)
    end
end 
function (o::Interpolator)(x)
    i⁻ = findinterval(o,x)
    return linearinterp{Nout(o)}(o,x,i⁻)
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





