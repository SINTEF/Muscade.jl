"""
    noFB

A constant, used by elements' `residual` or `lagrangian` as their 3rd output if they do provide 
any feedback to the solver (for example, on the reduction of the barrier parameter in interior point method).

Example:
`return L,noFB`

See also: [`noFB`](@ref)
"""
const noFB=nothing

zilch(a::SArray{S}) where{S} = zeros(SArray{S,𝕣})
zilch(a::AbstractArray)      = zeros(𝕣,size(a))
"""
    position = ∂0(X)

Used by elements' `residual` or `lagrangian` to extract the zero-th order time derivative
from the variables `X` and `U`.

See also: [`∂1`](@ref),[`∂2`](@ref),[`getsomedofs`](@ref)  
"""
∂0(y)   = y[1]

"""
    velocity = ∂1(X)

Used by elements' `residual` or `lagrangian` to extract the first order time derivative
from the variables `X` and `U`. Where the solver does not provide this derivative (e.g.
a static solver), the output is a vector of zeros.

See also: [`∂0`](@ref),[`∂2`](@ref),[`getsomedofs`](@ref)  
"""
∂1(y)   = length(y) ≥2 ? y[2] : zilch(y[1])

"""
    position = ∂2(X)

Used by elements' `residual` or `lagrangian` to extract the zero-th order time derivative
from the variables `X` and `U`. Where the solver does not provide this derivative (e.g.
a static solver), the output is a vector of zeros.

See also: [`∂0`](@ref),[`∂1`](@ref),[`getsomedofs`](@ref)  
"""
∂2(y)   = length(y) ≥3 ? y[3] : zilch(y[1])
∂n(n)   = (∂0,∂1,∂2)[n+1] # type unstable
∂n(y,ider) = length(y) ≥ider+1 ? y[ider+1] : zilch(y[1]) # slow

"""
    rotations = getsomedofs(X,[3,6])

Used by elements' `residual` or `lagrangian` to some degrees of freedom, and their
time derivatives, from the variables `X` and `U`. 

See also: [`∂0`](@ref),[`∂1`](@ref),[`∂2`](@ref)  
"""
getsomedofs(A::NTuple{Nder,SVector},ind) where{Nder} = ntuple(i->A[i][ind],Nder)

#const Dof{Ndof,Nder} = NTuple{Nder,SVector{Ndof}}

"""
    c = coord(node)

Used by element constructors to obtain the coordinates of a vector of Nodes handed by
Muscade to the constructor. `c` is accessed as 
```
c[inod][icoord]
```
where `inod` is the element-node number and `icoord` an index into a vector of coordinates.

Note that `c[inod]` points at the same memory as `nod[inod].coord`: do not mutate `c[inod]`!

See also: [`addnode!`](@ref), [`addelement!`](@ref), [`describe`](@ref), [`solve`](@ref)  
"""
coord(nod::AbstractVector{Node}) = [n.coord for n∈nod]


"""

    Muscade.doflist(::Type{E<:AbstractElement})

Elements must provide a method for `Muscade.doflist`.  

The method must take the element type as only input, and return
a `NamedTuple` with fieldnames `inod`,`class` and `field`.  The tuple-fields
are `NTuple`s of the same length.  For example
```
Muscade.doflist( ::Type{<:Turbine}) = (inod =(1   ,1   ,2        ,2        ),
                                       class=(:X  ,:X  ,:A       ,:A       ),
                                       field=(:tx1,:tx2,:Δseadrag,:Δskydrag))
```                                       
In `Λ`, `X`, `U` and `A` handed by Muscade to `residual` or `lagrangian`,
the dofs in the vectors will follow the order in the doflist. Element developers
are free to number their dofs by node, by field, or in any other way.

See also: [`Muscade.lagrangian`](@ref), [`Muscade.residual`](@ref), [`Muscade.nosecondorder`](@ref)  
"""
doflist(     ::Type{E}) where{E<:AbstractElement}  = muscadeerror(@sprintf("method 'Muscade.doflist' must be provided for elements of type '%s'\n",E))

"""
    nosecondorder(::Type{E<:AbstractElement})

Elements that define `residual` which would give excessive compilation and/or execution time if differentiated
to the second order can implement a method after the below pattern to limit differentiation to first order:   

    Muscade.nosecondorder(     ::Type{<:MyElementType}) = Val(true)
"""
nosecondorder(     ::Type{<:AbstractElement}) = Val(false)

"""
    @espy function Muscade.lagrangian(eleobj::MyElement,Λ,X,U,A,t,SP,dbg)
        ...
        return L,FB
    end

Elements must implement a method for `Muscade.lagrangian` or [`Muscade.residual`](@ref).

# Inputs
- `eleobj` an element object
- `Λ` a `SVector{nXdof,R} where{R<:Real}`, Lagrange multipliers (aka `δX` virtual displacements).
- `X` a `NTuple` of `SVector{nXdof,R} where{R<:Real}`, containing the Xdofs and, depending on the solver,
   their time derivatives. Use `x=∂0(X)`, `v=∂1(X)` and `a=∂2(X)` to safely obtain vectors of zeros
   where the solver leaves time derivatives undefined.
- `U` a `NTuple` of `SVector{nUdof,R} where{R<:Real}`, containing the Udofs and, depending on the solver,
   their time derivatives. Use `u=∂0(U)`, `̇u=∂1(U)` and `̈u=∂2(U)` to safely obtain vectors of zeros
   where the solver leaves time derivatives undefined.
- `A` a `SVector{nAdof,R} where{R<:Real}`.
- `t` a ``Real` containing the time.
- `SP` solver parameters (for example: the barrier parameter `γ` for 
  interior point methods).
- `dbg` a `NamedTuple` to be used _only_ for debugging purposes.

# Outputs
- `L` the lagrangian
- `FB` feedback from the element to the solver (for example: can `γ` be 
  reduced?). Return `noFB` of the element has no feedback to provide.

See also: [`Muscade.residual`](@ref), [`Muscade.doflist`](@ref), [`@espy`](@ref), [`∂0`](@ref), [`∂1`](@ref), [`∂2`](@ref), [`noFB`](@ref), 
"""
lagrangian()=nothing

"""
@espy function Muscade.residual(eleobj::MyElement,X,U,A,t,SP,dbg)
    ...
    return R,FB
end

Elements must implement a method for `Muscade.residual` or [`Muscade.lagrangian`](@ref).

# Inputs
- `eleobj` an element object
- `X` a `NTuple` of `SVector{nXdof,R} where{R<:Real}`, containing the Xdofs and, depending on the solver,
   their time derivatives. Use `x=∂0(X)`, `v=∂1(X)` and `a=∂2(X)` to safely obtain vectors of zeros
   where the solver leaves time derivatives undefined.
- `U` a `NTuple` of `SVector{nUdof,R} where{R<:Real}`, containing the Udofs and, depending on the solver,
   their time derivatives. Use `u=∂0(U)`, `̇u=∂1(U)` and `̈u=∂2(U)` to safely obtain vectors of zeros
   where the solver leaves time derivatives undefined.
- `A` a `SVector{nAdof,R} where{R<:Real}`.
- `t` a ``Real` containing the time.
- `SP` solver parameters (for example: the barrier parameter `γ` for 
  interior point methods).
- `dbg` a `NamedTuple` to be used _only_ for debugging purposes.

# Outputs
- `R` the residual
- `FB` feedback from the element to the solver (for example: can `γ` be 
  reduced?). Return `noFB` of the element has no feedback to provide.

See also: [`Muscade.lagrangian`](@ref), [`Muscade.nosecondorder`](@ref), [`Muscade.doflist`](@ref), [`@espy`](@ref), [`∂0`](@ref), [`∂1`](@ref), [`∂2`](@ref), [`noFB`](@ref)
"""
residual()=nothing
"""
    mut,opt = Muscade.allocate_drawdata(axis,eleobjs;kwargs...)

Elements that are to be displayed in graphical output must implement a method for `Muscade.allocate_drawdata`.

The method is to allocate `opt`, a `NamedTuple` of data that will not be mutated from frame to frame,
but are usefull in [`Muscade.update_drawdata`](@ref) or [`Muscade.draw!`](@ref).

The method is also to allocate `mut`, a `NamedTuple` of data that will be mutated from frame to frame.  When
implementing graphics with `GLMakie.jl`, the fields of `mut` must be exactly the updatable inputs provided to
`GLMakie.jl`'s drawing primitives: in [`Muscade.draw!`](@ref)
    lines!(axis,mut.x,mut.y)
is acceptable, but 
    lines!(axis,mut.x[:,s],mut.y[:,s])
    lines!(axis,mut.a.x,mut.a.y)
are not.

The content of `Arrays` in `opt` and `mut` can be `undef`-ined.

Inputs are:
- `axis` the "canvas" to draw on, typicaly a `GLMakie.jl` `Axis`.
- `eleobjs` an `AbstractVector` of element objects, of length `nel`.
- `kwargs` a `NamedTuple` containing the keyword arguments provided by the user. See [`default`](@ref).

See also: [`Muscade.update_drawdata`](@ref), [`Muscade.draw!`](@ref)
"""
allocate_drawdata(axis,::AbstractVector{E};kwargs...)                    where{E<:AbstractElement} = nothing,nothing # mut,opt
"""
    mut = Muscade.update_drawdata(  axis,::AbstractVector{E},oldmut,opt, Λ,X,U,A,t,SP,dbg)

Elements that are to be displayed in graphical output must implement a method for `Muscade.allocate_drawdata`.

For parametric element types
    Muscade.update_drawdata(axis,o::AbstractVector{Teleobj}, Λ,X,U,A,t,SP,dbg;kwargs...) where{Teleobj<:MyElement}
For non-parametric element types, one can simplify the above to:
    Muscade.update_drawdata(axis,o::AbstractVector{MyElement}, Λ,X,U,A,t,SP,dbg;kwargs...)

Inputs are:
- `axis` the "canvas" to draw on, typicaly a `GLMakie.jl` `Axis`.
- `eleobjs` an `AbstractVector` of element objects, of length `nel`.
- `oldmut` the output `mut` of `Muscade.allocate_drawdata` or of a previous call to `Muscade.update_drawdata`.
- `opt` the output `opt` of `Muscade.allocate_drawdata`
- `Λ` a matrix of size `(nXdof,nel)`
- `X` a `NTuple` (over the derivatives) of matrices of size `(nXdof,nel)`
- `U` a `NTuple` (over the derivatives) of matrices of size `(nUdof,nel)`
- `A` a matrix of size `(nAdof,nel)`
- `t` time
- `SP` solver parameters
- `dbg` debuging information
- `kwargs` a `NamedTuple` containing the keyword arguments provided by the user. See [`default`](@ref).

See also: [`Muscade.allocate_drawdata`](@ref), [`Muscade.draw!`](@ref)
"""
update_drawdata(axis,::AbstractVector{E},oldmut,opt, Λ,X,U,A,t,SP,dbg) where{E<:AbstractElement} = nothing         # mut
"""
    Muscade.draw!(axis,MyElement,mut,opt)

Elements that are to be displayed in graphical output must implement a method for `Muscade.draw!`.

Inputs are:
- `axis` the "canvas" to draw on, typicaly a `GLMakie.jl` `Axis`.
- `MyElement` used for dispatching to the right method.
- `mut` is a `NamedTuple`, as output by [`Muscade.update_drawdata`](@ref).  More specificaly, if implementing
   graphics with `GLMakie.jl`, `mut` has been
- `opt` is as returned by [`Muscade.allocate_drawdata`](@ref)

When implementing graphics with `GLMakie.jl`, the fields of `mut` must be exactly the updatable inputs provided to
`GLMakie.jl`'s drawing primitives: in [`Muscade.draw!`](@ref)
    lines!(axis,mut.x,mut.y)
is acceptable, but 
    lines!(axis,mut.x[:,s],mut.y[:,s])
or    
    lines!(axis,mut.a.x,mut.a.y)
are not. The reason is that when doing graphics with `GLMakie.jl`, `Muscade` will wrap each field of `mut`
into an `Observable` before calling the elements' methods `draw!`.  This allows `Muscade` to update the graphics 
by just calling [`Muscade.update_drawdata`](@ref) for each element.

See also: [`Muscade.allocate_drawdata`](@ref), [`Muscade.update_drawdata`](@ref)
"""
draw!(axis,::Type{E},obs,opt) where{E<:AbstractElement} = nothing         # nothing

