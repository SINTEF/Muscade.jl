# 3D rotation

# The saga of sinc1 and its derivatives (NB: module Muscade defines these functions, in adiff)
# sinc1(x) = sin(x)/x, while Julia defines sinc(x) = sin(πx)/πx
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
Muscade.@DiffRule1(Elements.sinc1,               Elements.sinc1′( a.x)                * a.dx )
Muscade.@DiffRule1(Elements.sinc1′,              Elements.sinc1″( a.x)                * a.dx )
Muscade.@DiffRule1(Elements.sinc1″,              Elements.sinc1‴( a.x)                * a.dx )
Muscade.@DiffRule1(Elements.sinc1‴,              Elements.sinc1⁗( a.x)                * a.dx )
Muscade.@DiffRule1(Elements.sinc1⁗,              Elements.sinc1⁗′(a.x)                * a.dx )


# sinc1(acos(x)), differentiable to fourth order over ]-1,1] 
function scac(x)
    dx = x-1  
    if abs(dx)>1e-2 
        sinc1(acos(x))
    else  # deliberately a long Taylor series (5th order): this function will be adiffed at least to 2nd order, up to 4th order
        y = 1 + dx*(1/3 + dx*(-2/90 + dx*(0.0052911879917544626 + dx*(-0.0016229317117234072 + dx*(0.0005625))))) 
    end
end
         

const Mat33{R}   = SMatrix{3,3,R,9}
const Vec3{R}    = SVector{3,R}

spin(  v::Vec3 ) = SMatrix{3,3}(0,v[3],-v[2],-v[3],0,v[1],v[2],-v[1],0)
spin⁻¹(m::Mat33) = SVector{3}(m[3,2]-m[2,3],m[1,3]-m[3,1],m[2,1]-m[1,2])/2
trace( m::Mat33) = m[1,1]+m[2,2]+m[3,3] 
Rodrigues⁻¹(m)   = spin⁻¹(m)/scac((trace(m)-1)/2)   # NB: is necessarily singular for π turn
function Rodrigues(v::Vec3) 
    S = spin(v)
    θ = norm(v)
    return I + sinc1(θ)*S + sinc1(θ/2)^2/2*S*S  
end

normalize(v)     = v/norm(v)
# create a rotation vector that acts on u to make it colinear with v.  Fails if |u|=0, |v|=0 or θ=π
function adjust(u::Vec3{R},v::Vec3{R}) where{R}
    u,v = normalize.((u,v))
    c,w = dot(u,v), cross(u,v) 
    s   = norm(w)
    θ   = atan(s,c)
    return w/sinc1(θ)
end