# 3D rotation

const Mat33{R}   = SMatrix{3,3,R,9}
const Vec3{R}    = SVector{3,R}

normalize(v)     = v/norm(v)
spin(  v::Vec3 ) = SMatrix{3,3}(0,v[3],-v[2],-v[3],0,v[1],v[2],-v[1],0)
spin⁻¹(m::Mat33) = SVector{3}(m[3,2]-m[2,3],m[1,3]-m[3,1],m[2,1]-m[1,2])/2
trace( m::Mat33) = sum(m[i,i] for i∈(1,2,3))
# sinc1(acos(x)), differentiable at x==1 
function scac(x)
    dx = x-1  # deliberately a long Taylor series: this function will be adiffed at least to 2nd order
    abs(dx)>1e-6 ? sinc1(acos(x)) :  1 + dx/3 - 2/90*dx^2 + 0.00529103*dx^3 - 0.0016223*dx^4 + 0.00056*dx^5 
end
Rodrigues⁻¹(m) = spin⁻¹(m)/scac((trace(m)-1)/2)
function Rodrigues(v::Vec3) 
    S = spin(v)
    θ = norm(v)
    return I + sinc1(θ)*S + sinc1(θ/2)^2/2*S*S  
end

# create a rotation vector that acts on u to make it colinear with v.  Fails if |u|=0, |v|=0 or θ=π
function adjust(u::Vec3{R},v::Vec3{R}) where{R}
    u,v = normalize.((u,v))
    c,w = dot(u,v), cross(u,v) 
    s   = norm(w)
    θ   = atan(s,c)
    return w/sinc1(θ)
end