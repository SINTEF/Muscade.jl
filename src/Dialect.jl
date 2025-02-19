
# by Philippe Mainçon
#
# In VScode, use the "Fast Unicode math characters" plugging
# To interrupt Julia, CTRL-j,k
using Printf,StaticArrays
import Base.Threads.@spawn, Base.Threads.nthreads

## Basic types
# abstract types for dispatch
"""
    𝔹 (\\bbB)

an alias for `Bool`. For use in dispatching.
`𝔹1`... `𝔹4` are `AbstractArrays` of dimensions 1 to 4.
"""        
const 𝔹  = Bool      # \bbB
"""
    ℕ (\\bbN)

an alias for `UInt64`. For use in dispatching.
`ℕ1`... `ℕ4` are `AbstractArrays` of dimensions 1 to 4.
"""        
const ℕ  = UInt64
"""
    ℤ (\\bbZ)

an alias for abstract type `Integer`. For use in dispatching.
`ℤ1`... `ℤ4` are `AbstractArrays` of dimensions 1 to 4.
`ℤ11` is an `AbstractVector` of `AbstractVector`.
"""        
const ℤ  = Integer
"""
    ℝ (\\bbR)

an alias for abstract type `Real`. For use in dispatching.
`ℝ1`... `ℝ4` are `AbstractArrays` of dimensions 1 to 4.
`ℝ11` is an `AbstractVector` of `AbstractVector`.
"""        
const ℝ  = Real
# concrete types for allocation
"""
    𝕓 (\\bbb)

an alias for `Bool`. For use in `struct` definitions.
`𝕓1`... `𝕓4` are `Arrays` of dimensions 1 to 4.
"""
const 𝕓  = Bool
"""
    𝕟 (\\bbn)

an alias for `UInt64`. For use in `struct` definitions.
`𝕟1`... `𝕟4` are `Arrays` of dimensions 1 to 4.
"""
const 𝕟  = UInt64
"""
    𝕫 (\\bbz)

an alias for `Int64`. For use in `struct` definitions.
`𝕫1`... `𝕫4` are `Arrays` of dimensions 1 to 4.
`𝕫11` is a `Vector` of `Vector`.
"""
const 𝕫  = Int64
"""
    𝕣 (\\bbr)

an alias for `Float64`. For use in `struct` definitions.
`𝕣1`... `𝕣4` are `Arrays` of dimensions 1 to 4.
`𝕣11` is a `Vector` of `Vector`.
"""
const 𝕣  = Float64
"""
    ϵ (\\epsilon)

an alias for `Base.eps(𝕣)`. 
"""
const ϵ  = Base.eps(𝕣)
"""
    ∞ (\\infty)

an alias for `Base.inf`. 
"""
const ∞  = Base.Inf

# define arrays of these
for T in (:𝔹,:ℕ,:ℤ,:ℝ)
    #@eval export $T
    @eval const  $(Symbol(T,:x)) = AbstractArray{t} where {t<: $T}
    #@eval export $(Symbol(T,:x))  
    for N in (:1,:2,:3,:4)
        TN = Symbol(T,N)
        @eval const  $TN{t} = AbstractArray{t,$N} where {t<: $T}
        #@eval export $TN
    end
end
for T in (:𝕓,:𝕟,:𝕫,:𝕣)
    #@eval export $T
    Ts = Symbol(T,:s)
    for N in (:1,:2,:3,:4)
        TN = Symbol(T,N)
        @eval const  $TN = Array{$T,$N}
        #@eval export $TN
    end
end
const ℝ11 = AbstractVector{A} where {A<:ℝ1}
const ℤ11 = AbstractVector{A} where {A<:ℤ1}
const 𝕣11 = Vector{Vector{𝕣}}
const 𝕫11 = Vector{Vector{𝕫}}

## Miscellaneous
subtypeof(a::AbstractVector,b::AbstractVector) = a[a .<: Union{b...}]
# Given a variable, or its type, e.g. SMatrix{S,T}, get the name of the constructor, e.g. SMatrix
constructor(T::DataType)               = T.name.wrapper
constructor(x::T) where{T}             = T.name.wrapper
"""

    toggle(condition,a,b)

Typestable equivalent of `condition ? a : b`.  
Returns a value converted to `promote_type(typeof(a),typeof(b))`
"""
toggle(cond::Bool,a::Ta,b::Tb) where{Ta,Tb} = convert(promote_type(Ta,Tb), cond ? a : b)
# macro toggle(cond,a,b) # evaluate only a or only b
#     return :(convert(promote_type(typeof($a),typeof($b)), $cond ? $a : $b))
# end
getval(::Val{v}) where{v} = v

## Array handling
flat(a)                                = reshape(a,length(a))
# flatten a vector of vectors of identical size, but lead a matrix as-is
consolidate(a) = a
consolidate(a::AbstractVector{E}) where{E<:AbstractVector{T}} where{T} = reduce(hcat,a)

# same as 'unique', but returns also idx, a vector of vectors, index into v of the uniques
function uniques(v::AbstractVector{T}) where{T}
    u   = Vector{T        }()
    idx = Vector{Vector{𝕫}}()
    for (i,x) ∈ enumerate(v)
        if x ∉ u
            push!(u  ,x                      )
            push!(idx,findall([x==w for w∈v]))
        end
    end
    return u,idx
end



function showtime(t)
    return if t<1e-6
        @sprintf " %3d [ns]" round(Int,t*1e9)
    elseif t<1e-3
        @sprintf " %3d [μs]" round(Int,t*1e6)
    elseif t<9
        @sprintf " %3d [ms]" round(Int,t*1e3)
    elseif t<9*3600
        @sprintf "%4d [s] " round(Int,t)
    else
        @sprintf "%4d [h] " round(Int,t/3600)
    end
end

# if a function f is given the argument pointer= Ref{SomeType}()
# the function can then do e.g. vec=allocate(pointer,Vector...) and write to vec.
# and the caller retrievs the data with vec = pointer[] 
# advantage over "return vec" is if f throws, then vec still contains some data.

const Pointer = Base.RefValue
#function allocate(pointer::Pointer{T},target::T) where{T}
function allocate(pointer::Pointer,target) # TODO use line above
    pointer[]=target
    return target
end

copies(n,a::T) where{T} = NTuple{n,T}(deepcopy(a) for i∈1:n)

using MacroTools: postwalk,gensym_ids,rmlines,unblock 
"""
    @once f(x)= x^2 
    
do not parse the definition of function `f` again if not modified.
Using in a script, this prevents recompilations in `Muscade` or applications
based on it when they receive such functions as argument
"""    
macro once(ex)
    ex  = postwalk(rmlines,ex)
    ex  = postwalk(unblock,ex)
    qex = QuoteNode(ex)
    ex  = esc(ex)
    tag = gensym("tag")
    return quote
        if ~@isdefined($tag) || $tag≠$qex
            $tag = $qex
            $ex    
        end 
    end
end
"""
    default{:fieldname}(namedtuple,defval)

attempt to get a field `fieldname` from a `NamedTuple`. If `namedtuple` does not have 
such a field - or is not a `NamedTuple`, return `defval`.
"""
struct default{S} end
default{S}(t::T,d=nothing) where{S,T<:NamedTuple} = hasfield(T,S) ? getfield(t,S) : d
default{S}(t::T,d=nothing) where{S,T            } =                                 d


"""
An "identity vector"

id = IdVec
i  = 9834987
id[i] == i # true for any i

See also Julia's `identity` function.
"""
struct IdVec end
const idvec = IdVec()
@inline Base.getindex(::IdVec,i) = i
