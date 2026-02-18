# Define the soil contact forces
struct SoilContact <: AbstractElement
    z₀ :: 𝕣
    Kh :: 𝕣
    Kv :: 𝕣
    Ch :: 𝕣
    Cv :: 𝕣
end
SoilContact(nod::Vector{Node};z₀=0.::𝕣,Kh=0.::𝕣,Kv=0.::𝕣,Ch=0.::𝕣,Cv=0.::𝕣) = SoilContact(z₀,Kh,Kv,Ch,Cv)
@espy function Muscade.residual(o::SoilContact, X,U,A, t,SP,dbg) 
    x,x′ = ∂0(X)[1], ∂1(X)[1]
    y,y′ = ∂0(X)[2], ∂1(X)[2]
    z,z′ = ∂0(X)[3], ∂1(X)[3]
    if z < o.z₀ 
        R         = SVector(o.Kh*x +o.Ch*x′,o.Kh*y +o.Ch*y′,o.Kv*(z-o.z₀)+o.Cv*z′)
    else 
        R         = SVector(0,0,0)
    end
    return R,noFB
end
Muscade.doflist( ::Type{SoilContact})  = (inod =(1 ,1, 1), class=(:X,:X,:X), field=(:t1,:t2,:t3))
