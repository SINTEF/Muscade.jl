
"""
`noχ`

A constant, used by elements' `residual` or `lagrangian` as their 2nd output if they do not export
any material memory (for example, plastic strain).

Example:
`return L,noχ,noFB`

See also: [`noFB`](@ref)
"""
const noχ =nothing
"""
`noFB`

A constant, used by elements' `residual` or `lagrangian` as their 3rd output if they do provide 
any feedback to the solver (for example, on the reduction of the barrier parameter in interior point method).

Example:
`return L,noχ,noFB`

See also: [`noFB`](@ref)
"""
const noFB=nothing

#∂n(Y,n) = n+1≤lastindex(Y) ? Y[n+1] : zeros(eltype(Y[1]),size(Y[1])...)  # TODO this implementation will be slow if zero is to be returned!
"""
`position = ∂0(X)`

Used by elements' `residual` or `lagrangian` to extract the zero-th order time derivative
from the variables `X` and `U`.

See also: [`∂1`](@ref),[`∂2`](@ref)  
"""
∂0(y)   = y[1]

"""
`velocity = ∂1(X)`

Used by elements' `residual` or `lagrangian` to extract the first order time derivative
from the variables `X` and `U`. Where the solver does not provide this derivative (e.g.
a static solver), the output is a vector of zeros.

See also: [`∂0`](@ref),[`∂2`](@ref)  
"""
∂1(y)   = length(y) ≥2 ? y[2] : SVector(0. for x∈y[1])

"""
`position = ∂2(X)`

Used by elements' `residual` or `lagrangian` to extract the zero-th order time derivative
from the variables `X` and `U`. Where the solver does not provide this derivative (e.g.
a static solver), the output is a vector of zeros.

See also: [`∂0`](@ref),[`∂1`](@ref)  
"""
∂2(y)   = length(y) ≥3 ? y[2] : SVector(0. for x∈y[1])
∂n(n)   = (∂0,∂1,∂2)[n+1]

"""
`c = coord(node)`

Used by element constructors to obtain the coordinates of a vector of Nodes handed by
Muscade to the constructor.

See also: [`addnode!`](@ref), [`addelement!`](@ref), [`describe`](@ref), [`solve`](@ref)  
"""
coord(nod::AbstractVector{Node}) = [n.coord for n∈nod]

"""
`doflist(::Type{E<:AbstractElement})`

Elements must overload Muscade's `doflist` function.  
The method must take the element type as only input, and return
a `NamedTuple` with fieldnames `inod`,`class` and `field`.  The tuple-fields
are `NTuple`s of the same length.  For example
`Muscade.doflist( ::Type{<:Turbine}) = (inod =(1   ,1   ,2        ,2        ),
                                       class=(:X  ,:X  ,:A       ,:A       ),
                                       field=(:tx1,:tx2,:Δseadrag,:Δskydrag))`
In `δX` (aka `Λ`), `X`, `U` and `A` handed by Muscade to `residual` or `lagrandian`,
the dofs in the vectors will follow the order in the doflist. Element developers
are free to number their dofs by node, by field, or in any other way.


See also: [`∂0`](@ref),[`∂1`](@ref),[`∂2`](@ref)  
"""
doflist(     ::Type{E}) where{E<:AbstractElement}  = muscadeerror(@sprintf("method 'Muscade.doflist' must be provided for elements of type '%s'\n",E))




