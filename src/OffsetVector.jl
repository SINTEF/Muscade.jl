"""
A simple offset vector that only implements "setindex!" and "getindex".
Based on Memory, so cannot be extended.
"""
struct OffsetVector{T}
    offset :: 𝕫
    a      :: Memory{T}
end
OffsetVector{T}(first::𝕫,last::𝕫) where{T} = OffsetVector(1-first,Memory{T}(undef,last-first+1))
OffsetVector{T}(         last::𝕫) where{T} = OffsetVector{T}(1,last)

function Base.setindex!(o::OffsetVector{T},x::T,i::ℤ) where{T}
    o.a[i +o.offset]  = x
end
function Base.setindex!(o::OffsetVector{T},x::T,i   ) where{T}
    o.a[i.+o.offset] .= x
end
Base.getindex(o::OffsetVector,i) = o.a[i.+o.offset]
Base.length(  o::OffsetVector  ) = length(o.a)

