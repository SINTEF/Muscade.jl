@functor  with() gap(v,t,θ₀,r₀)= SVector(   v[4] - (v[1] + r₀*(cos(θ₀+v[3])-cos(θ₀))), 
                                            v[5] - (v[2] + r₀*(sin(θ₀+v[3])-sin(θ₀))))
"""
    ExcentricRigidConnection <: AbstractElement

An element to constrain a rigid body motion between a master node with fields (:t1,t2,:r3) and an emissary node with fields (:t1,:t2)

# Example
```
using Muscade
model = Model(:TestModel)
MasterNode  =   addnode!(model,𝕣[0.,    0.,     0.])
EmissaryNode  = addnode!(model,𝕣[1/√2,  1/√2,   0.])
e     = addelement!(model,ExcentricRigidConnection,[MasterNode EmissaryNode])
```    

See also:  [`Hold`](@ref), [`DofConstraint`](@ref)
"""
struct ExcentricRigidConnection <: AbstractElement end  
function ExcentricRigidConnection(nod::Vector{Node}) 
    c  = coord(nod)
    xₘ = SVector{3}(c[1]) # coordinates of the master node
    xₑ = SVector{3}(c[2]) # coordinates of the emissary node
    θ₀ = atan( xₑ[2]-xₘ[2], xₑ[1]-xₘ[1] )
    r₀ = norm([xₑ[2]-xₘ[2], xₑ[1]-xₘ[1]])
    return DofConstraint{:X,   2  ,5, 0, 0, 
                        (2,2),  (:λt1,:λt2), 
                        (1,1,1,2,2),(:t1,:t2,:r3,:t1,:t2),
                        (),   (),    
                        (),   (),    
                        typeof(gap),typeof((r₀,θ₀)),typeof(Muscade.equal)}(gap,(r₀,θ₀),Muscade.equal)
end