# # 3D rotations
using LinearAlgebra, StaticArrays
using Muscade: sinc1
"""
    Toolbox.scac(x)

`scac(x) = sinc1(acos(x)),`  The function can be differentiated
to the fourth order over ]-1,1] .

See also [`Toolbox.sinc1`](@ref)
"""
function scac(x)
    dx = x-1  
    if abs(dx)>1e-3 
        sinc1(acos(x))
    else  # deliberately a long Taylor series (5th order): this function will be adiffed at least to 2nd order, up to 4th order
        evalpoly(dx,(1,1/3,-2/90,0.0052911879917544626,-0.0016229317117234072,0.0005625))
    end
end

const Mat33{R}   = SMatrix{3,3,R,9}
const Vec3{R}    = SVector{3,R}

"""
    Toolbox.spin(v::SVector{3})

Transform a rotation vector `v` into the cross product matrix `M`, such that
`M ∘₁ a = v × a`.

See also [`Toolbox.spin⁻¹`](@ref), [`Toolbox.Rodrigues`](@ref), [`Toolbox.Rodrigues⁻¹`](@ref).
"""
spin(v::Vec3) = SMatrix{3,3}(0,v[3],-v[2],-v[3],0,v[1],v[2],-v[1],0)
"""
    toolbox.spin⁻¹(M::SMatrix{3,3})

Transform a cross product matrix `M` into the rotation vector `v`, such that
`v × a = M ∘₁ a`.

See also [`Toolbox.spin`](@ref), [`Toolbox.Rodrigues`](@ref), [`Toolbox.Rodrigues⁻¹`](@ref).
"""
spin⁻¹(m::Mat33) = SVector{3}(m[3,2]-m[2,3],m[1,3]-m[3,1],m[2,1]-m[1,2])/2
"""
    trace(v::SMatrix{3,3})

Computes the trace of a matrix.
"""
trace(m::Mat33) = m[1,1]+m[2,2]+m[3,3] 
"""
    Toolbox.Rodrigues⁻¹(v::SVector{3})

Transform a rotation matrix `M` into the rotation vector `v`, such that
`|v| < π`. Undefined for rotations of angle `π`

See also [`Toolbox.spin`](@ref), [`Toolbox.spin⁻¹`](@ref), [`Toolbox.Rodrigues`](@ref), [`Toolbox.adjust`](@ref).
"""
Rodrigues⁻¹(m)   = spin⁻¹(m)/scac((trace(m)-1.)/2.)   # NB: is necessarily singular for π turn
function LinearAlgebra.norm(v::Vec3)  # executes faster, COMPILES MUCH FASTER , and adiffs poorly at origin
    n = sqrt(v[1]*v[1]+v[2]*v[2]+v[3]*v[3])
    if n<1e-14
        n = zero(eltype(v))
    end
    return n
end

function spin²(S::Vec3) 
    ab  = S[2,3]*S[3,1] 
    ca  = S[1,2]*S[2,3] 
    bc  = S[3,1]*S[1,2]
    ma² = S[3,2]*S[2,3]
    mb² = S[3,1]*S[1,3]
    mc² = S[2,1]*S[1,2]
    return SMatrix{3,3}(mc²+mb²,ab,ca, ab,ma²+mc²,bc, ca,bc,ma²+mb²)
end     

"""
    Toolbox.Rodrigues(v::SVector{3})

Transform a rotation vector `v` into the rotation matrix `M`.

See also [`Toolbox.spin`](@ref), [`Toolbox.spin⁻¹`](@ref), [`Toolbox.Rodrigues⁻¹`](@ref), [`Toolbox.adjust`](@ref).
"""
@inline function Rodrigues(v::Vec3{R}) where{R} 
    # I + spin(v)*sinc1(norm3(v)) + spin²(v)*sinc1(norm3(v)/2)/√2 
    θ           = norm(v)
    s₁          = sinc1(θ  )
    s₂          = sinc1(θ/2)*(1/√2)
    a,b,c       = s₁*v[1],s₁*v[2],s₁*v[3] 
    vs₂         = v*s₂
    ab ,ac ,bc  = vs₂[1]*vs₂[2], vs₂[1]*vs₂[3], vs₂[2]*vs₂[3]
    ma²,mb²,mc² = vs₂[1]*vs₂[1], vs₂[2]*vs₂[2], vs₂[3]*vs₂[3]
    return Mat33{R}(1. -mc²-mb²,ab+c,ac-b,   ab-c,1. -ma²-mc²,bc+a,   ac+b,bc-a,1. -ma²-mb²)
end

"""
    Toolbox.adjust(u::SVector{3},v::SVector{3})

Compute the matrix of the rotation with smallest angle that transforms `u` into a vector colinear with v.  
Fails if |u|=0, |v|=0 or if the angle of the rotation is π.

See also [`Toolbox.spin`](@ref), [`Toolbox.spin⁻¹`](@ref), [`Toolbox.Rodrigues`](@ref), [`Toolbox.Rodrigues⁻¹`](@ref).
"""
function adjust(u::Vec3{R},v::Vec3{R}) where{R}
    u,v = normalize.((u,v))
    c,w = dot(u,v), cross(u,v) 
    s   = norm(w)
    θ   = atan(s,c)
    return w/sinc1(θ)
end
"""
    Toolbox.intrinsicrotationrates(rₑ::NTuple{ND,SMatrix{3,3}}) where{ND}

Transform a `NTuple` containing a rotation matrix and its extrinsic time derivatives,
into a `NTuple` containing a (zero) rotation vector and its intrinsic time derivatives.

See also [`Toolbox.spin`](@ref), [`Toolbox.spin⁻¹`](@ref), [`Toolbox.Rodrigues`](@ref), [`Toolbox.Rodrigues⁻¹`](@ref).
"""
function intrinsicrotationrates(rₑ::NTuple{ND,SMatrix{3,3}}) where{ND}
    vᵢ₀ =              (SVector{3,𝕣}(0,0,0),                                                                           )
    vᵢ₁ = ND<2 ? vᵢ₀ : (vᵢ₀...             , spin⁻¹(∂0(rₑ)' ∘₁ ∂1(rₑ))                                                 ) 
    vᵢ  = ND<3 ? vᵢ₁ : (vᵢ₁...                                        ,   spin⁻¹(∂1(rₑ)' ∘₁ ∂1(rₑ) + ∂0(rₑ)' ∘₁ ∂2(rₑ)))  
    return vᵢ
end



;