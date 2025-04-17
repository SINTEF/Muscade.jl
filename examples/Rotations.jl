# # 3D rotations
using LinearAlgebra
"""
    BeamElements.sinc1(x)

`BeamElements.sinc1(x) = sin(x)/x` - but  `BeamElements.sinc1(0.) = 1.`.  The function can be differentiated
to the fourth order.

This differs from Julia's `sinc(x) = sin(π*x)/(π*x)`.

See also [`scac`](@ref)
"""
sinc1(x) = sinc(x/π) 
function sinc1′(x)
    if abs(x)>1e-3
        s,c=sin(x),cos(x)
        c/x -s/x^2
    else
        x² = x*x
        x*(-1/3 +x²/30) 
    end
end
function sinc1″(x)
    if abs(x)>1e-1
        s,c=sin(x),cos(x)
        -s/x -2c/x^2 +2s/x^3
    else
        x² = x*x
        -1/3 +x²*(1/10 +x²*(-1/168 +x²*(1/6480))) 
    end
end
function sinc1‴(x)
    if abs(x)>0.4
        s,c=sin(x),cos(x)
        -c/x +3s/x^2 +6c/x^3 -6s/x^4
    else
        x² = x*x
        x*(1/5 +x²*(-1/42 +x²*(1/1080 +x²*(-1/55440 +x²*(1/4717440)))))
    end
end
function sinc1⁗(x) 
    x² = x*x
    1/5 +x²*(-1/14 +x²*(1/216 +x²*(-1/7920 +x²*(1/524160 +x²*(-1/54432000 +x²*(1/54432000 +x²*(-1/8143027200 +x²*(1/1656387532800))))))))
end
sinc1⁗′(x) = x*NaN
using Muscade
Muscade.@DiffRule1(sinc1,               sinc1′( a.x)                * a.dx )
Muscade.@DiffRule1(sinc1′,              sinc1″( a.x)                * a.dx )
Muscade.@DiffRule1(sinc1″,              sinc1‴( a.x)                * a.dx )
Muscade.@DiffRule1(sinc1‴,              sinc1⁗( a.x)                * a.dx )
Muscade.@DiffRule1(sinc1⁗,              sinc1⁗′(a.x)                * a.dx )


## sinc1(acos(x)), differentiable to fourth order over ]-1,1] 
"""
    BeamElements.scac(x)

`BeamElements.scac(x) = BeamElements.sinc1(acos(x)),`  The function can be differentiated
to the fourth order over ]-1,1] .

See also [`sinc1`](@ref)
"""
function scac(x)
    dx = x-1  
    if abs(dx)>1e-3 
        sinc1(acos(x))
    else  # deliberately a long Taylor series (5th order): this function will be adiffed at least to 2nd order, up to 4th order
        y = 1 + dx*(1/3 + dx*(-2/90 + dx*(0.0052911879917544626 + dx*(-0.0016229317117234072 + dx*(0.0005625))))) 
    end
end

const Mat33{R}   = SMatrix{3,3,R,9}
const Vec3{R}    = SVector{3,R}

"""
    M = BeamElements.spin(v::SVector{3})

Transform a rotation vector `v` into the cross product matrix `M`, such that
`M ∘₁ a = v × a`.

See also [`spin⁻¹`](@ref), [`Rodrigues`](@ref), [`Rodrigues⁻¹`](@ref).
"""
spin(  v::Vec3 ) = SMatrix{3,3}(0,v[3],-v[2],-v[3],0,v[1],v[2],-v[1],0)
"""
    v = BeamElements.spin⁻¹(M::SMatrix{3,3})

Transform a cross product matrix `M` into the rotation vector `v`, such that
`v × a = M ∘₁ a`.

See also [`spin`](@ref), [`Rodrigues`](@ref), [`Rodrigues⁻¹`](@ref).
"""
spin⁻¹(m::Mat33) = SVector{3}(m[3,2]-m[2,3],m[1,3]-m[3,1],m[2,1]-m[1,2])/2
"""
    t = BeamElements.trace(v::SMatrix{3,3})

Computes the trace of a matrix.
"""
trace( m::Mat33) = m[1,1]+m[2,2]+m[3,3] 
"""
    M = BeamElements.Rodrigues⁻¹(v::SVector{3})

Transform a rotation matrix `M` into the rotation vector `v`, such that
`|v| < π`. Undefined for rotations of angle `π`

See also [`spin`](@ref), [`spin⁻¹`](@ref), [`Rodrigues`](@ref), [`adjust`](@ref).
"""
Rodrigues⁻¹(m)   = spin⁻¹(m)/scac((trace(m)-1)/2)   # NB: is necessarily singular for π turn
"""
    M = BeamElements.Rodrigues(v::SVector{3})

Transform a rotation vector `v` into the rotation matrix `M`.

See also [`spin`](@ref), [`spin⁻¹`](@ref), [`Rodrigues⁻¹`](@ref), [`adjust`](@ref).
"""
function Rodrigues(v::Vec3) 
    S = spin(v)
    θ = norm(v)
    return LinearAlgebra.I + sinc1(θ)*S + sinc1(θ/2)^2/2*S*S  
end
# """
#     v1 = BeamElements.normalize(v::SVector{3})

# Compute a unit vector of same direction as `v`.  Fails
# if `|v|==0`.
# """
# normalize(v)     = v/norm(v)
## create a rotation vector that acts on u to make it colinear with v.  Fails if |u|=0, |v|=0 or θ=π
"""
    M = BeamElements.adjust(u::SVector{3},v::SVector{3})

Compute the matrix of the rotation with smallest angle that transforms `u` into a vector colinear with v.  
Fails if |u|=0, |v|=0 or if the angle of the rotation is π.

See also [`spin`](@ref), [`spin⁻¹`](@ref), [`Rodrigues`](@ref), [`Rodrigues⁻¹`](@ref).
"""
function adjust(u::Vec3{R},v::Vec3{R}) where{R}
    u,v = normalize.((u,v))
    c,w = dot(u,v), cross(u,v) 
    s   = norm(w)
    θ   = atan(s,c)
    return w/sinc1(θ)
end;