

const noχ =nothing
const noFB=nothing

# MUST be used by elements to unpack X and U.  Today, the various derivatives are packed into tuples.  Would we use Adiff tomorrow, allowing
# correct computation of e.g. Coriolis terms in beam elements?
∂n(Y,n) = n+1≤lastindex(Y) ? Y[n+1] : zeros(eltype(Y[1]),size(Y[1])...)  # this implementation will be slow if zero is to be returned!
∂0(y)   = ∂n(y,0)
∂1(y)   = ∂n(y,1)
∂2(y)   = ∂n(y,2)

"""
`c = coord(node)`

Used by element constructors to obtain the coordinates of a vector of Nodes handed by
Muscade to the constructor.

See also: [`addnode!`](@ref), [`addelement!`](@ref), [`describe`](@ref), [`solve`](@ref)  
"""
coord(nod::AbstractVector{Node}) = [n.coord for n∈nod]


doflist(     ::Type{E}) where{E<:AbstractElement}  = muscadeerror(@sprintf("method 'Muscade.doflist' must be provided for elements of type '%s'\n",E))
### Not part of element API, not exported by Muscade




