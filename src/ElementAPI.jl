"""
    noFB

A constant, used by elements' `residual` or `lagrangian` as their 3rd output if they do provide 
any feedback to the solver (for example, on the reduction of the barrier parameter in interior point method).

Example:
`return L,noFB`

See also: [`noFB`](@ref)
"""
const noFB=nothing

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
∂1(y)   = length(y) ≥2 ? y[2] : zeros(SVector{length(y[1])})

"""
    position = ∂2(X)

Used by elements' `residual` or `lagrangian` to extract the zero-th order time derivative
from the variables `X` and `U`. Where the solver does not provide this derivative (e.g.
a static solver), the output is a vector of zeros.

See also: [`∂0`](@ref),[`∂1`](@ref),[`getsomedofs`](@ref)  
"""
∂2(y)   = length(y) ≥3 ? y[3] : zeros(SVector{length(y[1])})
∂n(n)   = (∂0,∂1,∂2)[n+1] # type unstable
∂n(y,ider) = length(y) ≥ider+1 ? y[ider+1] : zeros(SVector{length(y[1])}) # slow

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

Elements must overload Muscade's `doflist` function.  
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

See also: [`lagrangian`](@ref), [`residual`](@ref)  
"""
doflist(     ::Type{E}) where{E<:AbstractElement}  = muscadeerror(@sprintf("method 'Muscade.doflist' must be provided for elements of type '%s'\n",E))

"""
    @espy function Muscade.lagrangian(eleobj::MyElement,Λ,X,U,A,t,SP,dbg)
        ...
        return L,FB
    end

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

See also: [`residual`](@ref), [`doflist`](@ref), [`@espy`](@ref), [`∂0`](@ref), [`∂1`](@ref), [`∂2`](@ref), [`noFB`](@ref)
"""
lagrangian()=nothing

"""
@espy function Muscade.residual(eleobj::MyElement,X,U,A,t,SP,dbg)
    ...
    return R,FB
end


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

See also: [`lagrangian`](@ref), [`doflist`](@ref), [`@espy`](@ref), [`∂0`](@ref), [`∂1`](@ref), [`∂2`](@ref), [`noFB`](@ref)
"""
residual()=nothing

"""
    P = constants(X,U,A,t)
    x = Muscade.motion{P}(X)

Transform a `NTuple` of `SVector`s, for example the vector `X` provided as an input to
`residual` or `Lagrangian` into a `SVector` of `∂ℝ`.  This can be used by an element to 
compute time derivatives, for example Euler, Coriolis and centrifugal accelerations, 
or strain rates.

Some principles of safe automatic differentiation must be adhered to:
- the function that uses `Muscade.motion` must also 'unpack' : no variable that is touched by 
  the output of `Muscade.motion` must be returned by the function without having been unpacked
  by `Muscade.position`, `Muscade.velocity` or `Muscade.acceleration`.
- The precendence `P` must be calculated using `constants` with all variables that are input to 
  the function and may be differentiated.
- If other levels of automatic differentiation are introduced within the function, unpack in reverse
  order of packing.    

See [`Muscade.position`](@ref), [`Muscade.velocity`](@ref), [`Muscade.acceleration`](@ref)
"""
struct motion{P,D}      end 
struct motion_{P,Q}     end 
motion{ P,D}(  a::NTuple{D,SV{N,R}}) where{D,P,N,R        } = SV{N}(motion_{P,P+D-1}(ntuple(j->a[j][i],D)) for i=1:N)
motion_{P,Q  }(a::NTuple{D,     R }) where{D,P  ,R<:Real,Q} = ∂ℝ{Q,1}(motion_{P,Q-1}(a),SV(motion_{P,Q-1}(a[2:D]))) 
motion_{P,P  }(a::NTuple{D,     R }) where{D,P  ,R<:Real  } = a[1]
struct position{P,Q}     end 
struct velocity{P,Q}     end 
struct acceleration{P,Q} end 

"""
    P = constants(X,U,A,t)
    N = length(X)
    x = Muscade.motion{P,ND}(X)
    E = f(x)    
    ε = Muscade.position{P,N}(E)

Extract the position from a variable that is a function of the output of `Muscade.motion`.

See [`Muscade.motion`](@ref), [`Muscade.velocity`](@ref), [`Muscade.acceleration`](@ref)
"""
function position{P,Q}(a::ℝ) where{P,Q}
    if     Q==1                            a 
    elseif Q==2               value{P+Q-1}(a)
    elseif Q==3  value{P+Q-2}(value{P+Q-1}(a))
    else muscadeerror((P=P,Q=Q,a=typeof(a)),"'position' requires Q∈{1,2,3}")
    end
end
  
"""
    P = constants(X,U,A,t)
    N = length(X)
    x = Muscade.motion{P,N}(X)
    E = f(x)    
    ̇ε = Muscade.velocity{P,N}(E)

Extract the velocity or rate from a variable that is a function of the output of `Muscade.motion`.

See [`Muscade.motion`](@ref), [`Muscade.position`](@ref), [`Muscade.acceleration`](@ref)
"""
function velocity{P,Q}(a::ℝ) where{P,Q}
    if     Q==1                          0.   
    elseif Q==2               ∂{P+Q-1,1}(a)[1]
    elseif Q==3  value{P+Q-2}(∂{P+Q-1,1}(a)[1]) 
    else muscadeerror((P=P,Q=Q,a=typeof(a)),"'velocity' requires Q∈{1,2,3}")
    end
end  
"""
    P = constants(X,U,A,t)
    N = length(X)
    x = Muscade.motion{P,N}(X)
    E = f(x)    
    ̈ε = Muscade.acceleration{P,N}(E)

Extract the velocity or rate from a variable that is a function of the output of `Muscade.motion`.

See [`Muscade.motion`](@ref), [`Muscade.position`](@ref), [`Muscade.velocity`](@ref)
"""
function acceleration{P,Q}(a::ℝ) where{P,Q}
    if     Q==1                        0.
    elseif Q==2                        0.
    elseif Q==3  ∂{P+Q-2,1}(∂{P+Q-1,1}(a)[1])[1]
    else muscadeerror((P=P,Q=Q,a=typeof(a)),"'acceleration' requires Q∈{1,2,3}")
    end
end  

#acceleration{P,Q}(a::ℝ) where{P,Q} = ∂{P+Q-2,1}(∂{P+Q-1,1}(a)[1])[1]
position{P,Q}(    a::AbstractVector) where{P,Q} = position{P,Q}.(a)
velocity{P,Q}(    a::AbstractVector) where{P,Q} = velocity{P,Q}.(a)
acceleration{P,Q}(a::AbstractVector) where{P,Q} = acceleration{P,Q}.(a)

