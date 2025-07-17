include("Rotations.jl")
using StaticArrays, LinearAlgebra
using Muscade

struct FollowerPointLoad{Tf,Tm} <: AbstractElement
    f :: Tf
    m :: Tm
end
Muscade.nosecondorder(::Type{<:FollowerPointLoad}) = Val(true)
Muscade.doflist(::Type{<:FollowerPointLoad}) = (inod = (1,1,1,1,1,1),class= ntuple(i->:X,ndof),field= (:t1,:t2,:t3,:r1,:r2,:r3) )
FollowerPointLoad(nod::Vector{Node};f,m) = FollowerPointLoad(f,m)
@espy function Muscade.residual(o::FollowerPointLoad,   X,U,A,t,SP,dbg) 
    i = SVector(4,5,6)
    X_                  = motion{P}(X)
    V_                  = X_[i]


    P,ND                = constants(X),length(X)
    X₀                  = ∂0(X)
    TX₀                 = revariate{1}(X₀)
    Tgp,Tε,Tvₛₘ,Trₛₘ,Tvₗ₂ = kinematics(o,TX₀)
    X_                  = motion{P}(X)
    ☼ε ,ε∂X₀            = composewithJacobian{P,ND}(Tε  ,X_)
    vₛₘ∂X₀               = composeJacobian{    P,  }(Tvₛₘ,X₀)
    rₛₘ                  = composevalue{       P,ND}(Trₛₘ,X_)
    ♢κ                  = composevalue{       P,ND}(Tvₗ₂,X_).*(2/o.L)  # ♢: evaluate only on request
    vᵢ₀                 = (SVector(0,0,0),)
    vᵢ₁                 = ND≥1 ? (vᵢ₀...,   spin⁻¹(∂0(rₛₘ)' ∘₁ ∂1(rₛₘ))) : vᵢ₀ 
    vᵢ                  = ND≥2 ? (vᵢ₁...,   spin⁻¹(∂1(rₛₘ)' ∘₁ ∂1(rₛₘ) + ∂0(rₛₘ)' ∘₁ ∂2(rₛₘ))) : vᵢ₁
    gp                  = ntuple(ngp) do igp
        Tx,Tκ           = Tgp[igp].x, Tgp[igp].κ
        ☼x,x∂X₀         = composewithJacobian{P,ND}(Tx  ,X_)
        ☼κ,κ∂X₀         = composewithJacobian{P,ND}(Tκ  ,X_)
        fᵢ,mᵢ,fₑ,mₑ     = ☼resultants(o.mat,ε,κ,x,rₛₘ,vᵢ)          
        R               = (fᵢ ∘₀ ε∂X₀ + mᵢ ∘₁ κ∂X₀ + fₑ ∘₁ x∂X₀ + mₑ ∘₁ vₛₘ∂X₀) * o.dL[igp]     
        @named(R)
    end
    R                   = sum(gpᵢ.R for gpᵢ∈gp) 
    return R,noFB  
end
function kinematics(o::FollowerPointLoad,X₀)  
    cₘ,rₘ,tgₘ,tgₑ,ζnod,ζgp,L  = o.cₘ,o.rₘ,o.tgₘ,o.tgₑ,o.ζnod,o.ζgp,o.L   # As-meshed element coordinates and describing tangential vector
    ## transformation to corotated system
    uᵧ₁,vᵧ₁,uᵧ₂,vᵧ₂        = vec3(X₀,1:3), vec3(X₀,4:6), vec3(X₀,7:9), vec3(X₀,10:12)
    vₗ₂,rₛₘ,vₛₘ              = fast(SVector(vᵧ₁...,vᵧ₂...)) do v
        vᵧ₁,vᵧ₂            = vec3(v,1:3), vec3(v,4:6)
        rₛ₁                = fast(Rodrigues,vᵧ₁)
        rₛ₂                = fast(Rodrigues,vᵧ₂)
        vₗ₂                = 0.5*Rodrigues⁻¹(rₛ₂ ∘₁ rₛ₁')
        rₛₘ                = fast(Rodrigues,vₗ₂) ∘₁ rₛ₁ ∘₁ o.rₘ  
        vₛₘ                = Rodrigues⁻¹(rₛₘ)              
        return vₗ₂,rₛₘ,vₛₘ
    end   
    cₛ                     = 0.5*(uᵧ₁+uᵧ₂)
    uₗ₂                    = rₛₘ'∘₁(uᵧ₂+tgₘ*ζnod[2]-cₛ)-tgₑ*ζnod[2]    #Local displacement of node 2
    ## interpolation
    ε                     = √((uₗ₂[1]+L/2)^2+uₗ₂[2]^2+uₗ₂[3]^2)*2/L - 1.       
    gp                    = ntuple(ngp) do igp
        yₐ,yᵤ,yᵥ,κₐ,κᵤ,κᵥ = o.yₐ[igp],o.yᵤ[igp],o.yᵥ[igp],o.κₐ[igp],o.κᵤ[igp],o.κᵥ[igp]
        κ                 = SVector(         κₐ*vₗ₂[1], κᵤ*uₗ₂[2]+κᵥ*vₗ₂[3], κᵤ*uₗ₂[3]-κᵥ*vₗ₂[2])  
        y                 = SVector(yₐ*uₗ₂[1]         , yᵤ*uₗ₂[2]+yᵥ*vₗ₂[3], yᵤ*uₗ₂[3]-yᵥ*vₗ₂[2])                              
        x                 = rₛₘ∘₁(tgₑ*ζgp[igp]+y)+cₛ+cₘ 
        (κ=κ,x=x)
    end
    return gp,ε,vₛₘ,rₛₘ,vₗ₂
end

