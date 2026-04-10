
# by Philippe Mainçon
#
# In VScode, use the "Fast Unicode math characters" plugging
# To interrupt Julia, CTRL-j,k
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
"""
    ℂ (\\bbC)

an alias for abstract type `Complex{<:Real}`. For use in dispatching.
`ℂ1`... `ℂ4` are `AbstractArrays` of dimensions 1 to 4.
`ℂ11` is an `AbstractVector` of `AbstractVector`.
"""        
const ℂ  = Complex{<:Real}


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
    𝕔 (\\bbc)

an alias for `Complex{Float64}`. For use in `struct` definitions.
`𝕔1`... `𝕔4` are `Arrays` of dimensions 1 to 4.
`𝕔11` is a `Vector` of `Vector`.
"""
const 𝕔  = Complex{Float64}
"""
    ϵ (\\epsilon)

an alias for `Base.eps(𝕣)`. 
"""
const ϵ  = Base.eps(𝕣)
"""
    ∞ (\\infty)

an alias for `Base.inf`. 
"""
const ∞      = Base.Inf  # \infty
const 𝑖      = im        # \iti
const ℜ     =  real     # \Re 
const ℑ     =  imag     # \Im
const expπ𝑖 = cispi  
const exp𝑖  = cis
"""
    𝕫log2(i::𝕫)

Compute the integer `log2` of an integer, fast. Fails if `i` is not a power of two.
"""
function 𝕫log2(i::𝕫) 
    a = 63-leading_zeros(i)
    b = trailing_zeros(i) 
    a==b || error("Input must be a power of 2")
    return a
end

# define arrays of these
for T in (:𝔹,:ℕ,:ℤ,:ℝ,:ℂ)
    #@eval export $T
    @eval const  $(Symbol(T,:x)) = AbstractArray{t} where {t<: $T}
    #@eval export $(Symbol(T,:x))  
    for N in (:0,:1,:2,:3,:4)
        TN = Symbol(T,N)
        @eval const  $TN{t} = AbstractArray{t,$N} where {t<: $T}
        #@eval export $TN
    end
end
for T in (:𝕓,:𝕟,:𝕫,:𝕣,:𝕔)
    #@eval export $T
    for N in (:1,:2,:3,:4)
        TN = Symbol(T,N)
        @eval const  $TN = Array{$T,$N}
        #@eval export $TN
    end
end
for T in (:𝕓,:𝕟,:𝕫,:𝕣,:𝕔)
    T0 = Symbol(T,:0)
    @eval const  $T0 = Ref{$T}
end
const ℝ11      = AbstractVector{A} where {A<:ℝ1}
const ℤ11      = AbstractVector{A} where {A<:ℤ1}
const ℂ11      = AbstractVector{A} where {A<:ℂ1}
const 𝕣11      = Vector{Vector{𝕣}}
const 𝕫11      = Vector{Vector{𝕫}}
const 𝕔11      = Vector{Vector{𝕔}}
const Sparse𝕣2 = SparseMatrixCSC{𝕣,𝕫}
const Sparse𝕔2 = SparseMatrixCSC{𝕔,𝕫}
const LU𝕣      = SparseArrays.UMFPACK.UmfpackLU{𝕣,𝕫}

## Miscellaneous
subtypeof(a::AbstractVector,b::AbstractVector) = a[a .<: Union{b...}]
# Given a variable, or its type, e.g. SMatrix{S,T}, get the name of the constructor, e.g. SMatrix
constructor(T::DataType)               = T.name.wrapper
constructor(x::T) where{T}             = T.name.wrapper
"""

    Muscade.toggle(condition,a,b)

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

# if a function f is given the argument ref= Ref{SomeType}()
# the function can then do e.g. vec=allocate(ref,Vector...) and write to vec.
# and the caller retrievs the data with vec = ref[] 
# advantage over "return vec" is if f throws, then vec still contains some data.

function allocate(ref::Ref{R},target::T) where{R,T<:R}
    ref[]=target
    return target
end


copies(n,a::T) where{T} = NTuple{n,T}(deepcopy(a) for i∈1:n)

"""
    value      = default{:fieldname}(namedtuple,defaultvalue)
    namedtuple = default(inputnamedtuple,defaultnamedtuple)

The first syntax attempts to access field `fieldname` from `namedtuple`. If `namedtuple` does not have 
such a field - or is not a `NamedTuple`, return `defaultvalue`.

The second syntax creates `namedtuple` from `inputnamedtuple`, supplementing with fields and values from `defaultnamedtuple` where
there is no corresponding field in `inputnamedtuple`.  This a thin wrapper of Julia's `Base.merge`. 
"""
struct default{S} end
default{S}(t::T,d=nothing) where{S,T<:Base.Pairs} = default{S}((;t...),d)
default{S}(t::T,d=nothing) where{S,T<:NamedTuple} = hasfield(T,S) ? getfield(t,S) : d
default{S}(t::T,d=nothing) where{S,T            } =                                 d
default(in,def) = Base.merge(def,in)

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

"""
a*idmult == idmult*a = a
and this should compile as a no-op.
"""

struct IdMult end
const  idmult = IdMult()
@inline Base.:(*)(::IdMult,a) = a
@inline Base.:(*)(a,::IdMult) = a
@inline Base.:(/)(a,::IdMult) = a



"""
    Muscade.mod_onebased(i,n) = mod(i-1,n)+1

For `i::ℤ`, returns a value in `{1,...n}`.  This differs
from `mod` which return a value in `[0,n[`   
"""
mod_onebased(i::ℤ,n) = mod(i-1,n)+1

"""
    columnmatrix(v)

Reshape a vector into a matrix of size `(length(v),1)`    
"""
columnmatrix(v::Vector) = reshape(v,(length(v),1))
"""
    rowmatrix(v)

Reshape a vector into a matrix of size `(1,length(v))`    
"""
rowmatrix(   v::Vector) = reshape(v,(1,length(v)))

"""
    colnormalize(a)

Euclidian-normalize the columns of an SMatrix
"""    
function colnormalize(a::SMatrix{ndim,nvec,R}) where{ndim,nvec,R<:ℝ}
    n   = SVector{     nvec,R}(norm(a[:,ivec]) for ivec=1:nvec)
    out = SMatrix{ndim,nvec,R}(a[idim,ivec]/n[ivec] for idim=1:ndim,ivec=1:nvec)
end
