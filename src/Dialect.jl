
# by Philippe Mainçon
#
# To get special characters in ATOM, type their "code", e.g. \partial, followed by ctrl-space (to get a suggestion) and TAB (or RETURN) to accept it.
# In VScode, use the "Fast Unicode math characters" plugging
# To interrupt Julia, CTRL-j,k
using Printf,StaticArrays
import Base.Threads.@spawn, Base.Threads.nthreads

## Basic types
# abstract types for dispatch
const 𝔹  = Bool      # \bbB
const ℕ  = UInt64
const ℤ  = Integer
const ℝ  = Real
# concrete types for allocation
const 𝕓  = Bool
const 𝕟  = UInt64
const 𝕫  = Int64
const 𝕣  = Float64
const ϵ  = Base.eps(𝕣)
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
# typestable equivalent of a ? b : c
toggle(cond::Bool,a::Ta,b::Tb) where{Ta,Tb} = convert(promote_type(Ta,Tb), cond ? a : b)
macro toggle(cond,a,b) # evaluate only a or only b
    return :(convert(promote_type(typeof($a),typeof($b)), $cond ? $a : $b))
end
getval(::Val{v}) where{v} = v

## Array handling
flat(a)                                = reshape(a,length(a))
# Take slice i from the d'th dimension of array a
@generated function Base.selectdim(a,::Val{d},i) where{d}
    precols = ()
    pstcols = ()
    for i = 1:d-1
        precols = (precols...,:)
    end
    for i = 1:ndims(a)-d
        pstcols = (pstcols...,:)
    end
    return quote
        return view(a,$(precols...),i,$(pstcols...))
    end
end
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

## Rear: indexing into the last index of an array
@generated function rearview(a,i)
    colons = ()
    for i = 1:ndims(a)-1
        colons = (colons...,:)
    end
    return quote
        return view(a,$(colons...),i)
    end
end
@generated function rearget(a,i)
    colons = ()
    for i = 1:ndims(a)-1
        colons = (colons...,:)
    end
    return quote
        return a[$(colons...),i]
    end
end
rearalloc(siz::NTuple{N, Any},el::E) where{E,N} = Array{E,N}(undef,siz)
@generated function rearalloc(siz::NTuple{Nsiz, Any},el::AbstractArray{E,Nel}) where{Nsiz,Nel,E}
    N = Nel+Nsiz
    return quote
        return Array{E,$N}(undef,(size(el)...,siz...))
    end
end
function rearset!(a::Array{E,Na},i::ℤ,b::Array{E,Nb}) where{E,Na,Nb}
    rearview(a,i)[:] = b
    return nothing
end
function rearset!(a::Vector{E},i::ℤ,b::E) where{E}
    rearview(a,i)[]  = b
    return nothing
end

function showtime(t)
    return if t<1e-6
        @sprintf " %3d [ns]" round(Int,t*1e9)
    elseif t<1e-3
        @sprintf " %3d [μs]" round(Int,t*1e6)
    elseif t<1
        @sprintf " %3d [ms]" round(Int,t*1e3)
    elseif t<3600
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

# @once f f(x)= x^2 # do not parse f again if not modified (prevent recompilation when passing function as arg from script)
using MacroTools: postwalk,gensym_ids,rmlines,unblock 
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


