include("Rotations.jl")

# # Euler beam element

using StaticArrays, LinearAlgebra, Muscade

# Data structure containing the cross section material properties
struct BeamCrossSection
    EA :: 𝕣
    EI :: 𝕣
    GJ :: 𝕣
    ## ρ  :: 𝕣 
    ## μ  :: 𝕣 
    ## Add moment of inertia about x for dynamic torque
    ## Cd :: SVector{3,𝕣}
    ## Ca :: SVector{3,𝕣}
    ## A  :: SVector{3,𝕣}
end
BeamCrossSection(;EA=EA,EI=EI,GJ=GJ) = BeamCrossSection(EA,EI,GJ);

# Resultant function that computes the internal loads from the strains and curvatures, and external loads on the element. 
@espy function resultants(o::BeamCrossSection,ε,κ,xᵧ,rₛₘ,vᵢ) 

    r₀  = ∂0(rₛₘ)  # orientation of the element's local refsys
    vᵢ₁ = ∂1(vᵢ)  # intrinsic rotation rate         of the element's local refsys
    vᵢ₂ = ∂2(vᵢ)  # intrinsic rotation acceleration of the element's local refsys
    ☼mₑ = SVector(0.,0.,0.) # external couples at Gauss point. mₑ is in local coordinates 


    xᵧ₀,xᵧ₁,xᵧ₂ = ∂0(xᵧ),∂1(xᵧ),∂2(xᵧ)
    xₗ₁          = xᵧ₁ ∘₁ r₀
    xₗ₂          = xᵧ₂ ∘₁ r₀
    ## Compute drag force (hard-coded parameters so far)
    # ρ = 1025.0
    # A  = SVector(0.0,1.0,1.0)
    # Cd = SVector(0.0,1.0,1.0) # SVector(0.0,0.0,0.0)
    # fd = .5 * ρ * A .* Cd .* xₗ₁ #.* abs.(xₗ₁) #mind the sign: forces exerted by element on its environment
    # ## Compute inertia force (hard-coded parameter so far)
    μ   = (1.0,1.0,1.0)
    fi = μ .* xₗ₂ 
    ## Compute added mass force (hard-coded parameter so far)
    # Ca = SVector(0.0,0.0,0.0)
    # fa = ρ * Ca .* xₗ₂
    
    # ☼fₑ = fd+fa+
    ☼fₑ = fi #SVector(0.,0.,0.) # external forces at Gauss point (no external moment/torque/... so far). fₑ is in local coordinates 
    ☼fᵢ = o.EA*∂0(ε)

    ## WARNING: curvatures are defined as rate of rotation along the element, not second derivatives of deflection.  
    ## Hence κ[3]>0 implies +2 direction is inside curve, 
    ##       κ[2]>0 implies -3 direction is inside curve.
    ☼mᵢ  = SVector(o.GJ*∂0(κ)[1],o.EI*∂0(κ)[2],o.EI*∂0(κ)[3])# replace by κ₀ 
    return fᵢ,mᵢ,fₑ,mₑ
end;

## Static Euler beam element, with two nodes, two Gauss points and 12 degrees of freedom. 
const ngp        = 4
const ndim       = 3
const ndof       = 12
const nnod       = 2;

# Shape functions for a beam element with support ζ∈[-1/2,1/2]. Though the shape function matrices are sparse, do not "unroll" them.  That would be faster but considerably clutter the code                          
yₐ(ζ) =            2ζ       # differential axial displacement or roll field
yᵤ(ζ) =  -4ζ^3    +3ζ       # deflection due to differential nodal transverse translation
yᵥ(ζ) =        ζ^2   - 1/4  # deflection due to differenttial rotation (bending, not torsion)
κₐ(ζ) =                2    # torsion  . κₐ = yₐ′ . Divide by L .    
κᵤ(ζ) =  -24ζ               # curvature. κᵤ = yᵤ′′. Divide by L².
κᵥ(ζ) =                2;   # curvature. κᵥ = yᵥ′′. Divide by L .

# Data structure describing an EulerBeam3D element as meshed
struct EulerBeam3D{Mat} <: AbstractElement
    cₘ       :: SVector{3,𝕣}    # Position of the middle of the element, as meshed
    rₘ       :: Mat33{𝕣}        # Orientation of the element, as meshed, represented by a rotation matrix (from global to local)
    ζgp      :: SVector{ngp,𝕣}  # Location of the Gauss points for the normalized element with length 1
    ζnod     :: SVector{nnod,𝕣} # Location of the nodes for the normalized element with length 1
    tgₘ      :: SVector{ndim,𝕣} # Vector connecting the nodes of the element in the global coordinate system
    tgₑ      :: SVector{ndim,𝕣} # Vector connecting the nodes of the element in the local coordinate system
    yₐ       :: SVector{ngp,𝕣}  # Value at gp of shape function for differential axial displacement or roll field
    yᵤ       :: SVector{ngp,𝕣}  # Value at gp of shape function for deflection due to differential nodal transverse translation
    yᵥ       :: SVector{ngp,𝕣}  # Value at gp of shape function for deflection due to differenttial rotation (bending, not torsion)
    κₐ       :: SVector{ngp,𝕣}  # Value at gp of shape function for torsion  . κₐ = yₐ′ . Divided by L .    
    κᵤ       :: SVector{ngp,𝕣}  # Value at gp of shape function for curvature. κᵤ = yᵤ′′. Divided by L².
    κᵥ       :: SVector{ngp,𝕣}  # Value at gp of shape function for curvature. κᵥ = yᵥ′′. Divided by L .
    L        :: 𝕣
    dL       :: SVector{ngp,𝕣}  # length associated to each Gauss point
    mat      :: Mat # used to store material properties (BeamCrossSection, for example)
end;

# For performance, `residual` will only accept differentiation to first order
Muscade.nosecondorder(::Type{<:EulerBeam3D}) = Val(true)

# Define nodes, classes, and field names of dofs
Muscade.doflist(     ::Type{<:EulerBeam3D}) = 
        (inod = (1,1,1,1,1,1, 2,2,2,2,2,2), 
         class= ntuple(i->:X,ndof), 
         field= (:t1,:t2,:t3,:r1,:r2,:r3, :t1,:t2,:t3,:r1,:r2,:r3) )

# Constructor for the EulerBeam3D element. Arguments: node list, material, and direction of the first bending axis in the global coordinate system.  
function EulerBeam3D(nod::Vector{Node};mat,orient2::SVector{ndim,𝕣}=SVector(0.,1.,0.))
    c       = coord(nod)
    ## Position of the middle of the element in the global coordinate system (as-meshed)
    cₘ      = SVector{ndim}((c[1]+c[2])/2)
    ## Length and tangential vector to the element in the global coordinate system  
    tgₘ     = SVector{ndim}( c[2]-c[1]   )
    L       = norm(tgₘ)
    t       = tgₘ/L
    ## Create t, n, b which are the longitudinal and two transverse unit vectors to the element (as-meshed). 
    ## NB: orient2, provided by the user, will define the first bending axis. 
    orient2/= norm(orient2)
    n       = orient2 - t*dot(orient2,t) 
    nn      = norm(n) 
    nn>1e-3 || muscadeerror("Provide a 'orient' input that is not nearly parallel to the element")
    n      /= nn
    b       = cross(t,n)
    rₘ      = SMatrix{ndim,ndim}(t...,n...,b...)
    ## Tangential vector and node coordinates in the local coordinate system
    tgₑ     = SVector{ndim}(L,0,0)
    ## Weight associated to each Gauss point
    # dL      = SVector{ngp }(L/2   , L/2 )
    dL    = SVector{ngp}(L/2*(18-sqrt(30))/36,L/2*(18+sqrt(30))/36  ,L/2*(18+sqrt(30))/36,L/2*(18-sqrt(30))/36  ) 
    ## Location ζgp of the Gauss points for a unit-length beam element, with nodes at ζnod=±1/2. 
    # ζgp     = SVector{ngp }(-1/2√3,1/2√3) 
    ζgp     = SVector{ngp }(-1/2*sqrt(3/7+2/7*sqrt(6/5)),-1/2*sqrt(3/7-2/7*sqrt(6/5)), +1/2*sqrt(3/7-2/7*sqrt(6/5)),+1/2*sqrt(3/7+2/7*sqrt(6/5))) 
    ζnod    = SVector{nnod }(-1/2  ,1/2  )
    shapes  = (yₐ.(ζgp), yᵤ.(ζgp), yᵥ.(ζgp)*L, κₐ.(ζgp)/L, κᵤ.(ζgp)/L^2, κᵥ.(ζgp)/L)
    return EulerBeam3D(cₘ,rₘ,ζgp,ζnod,tgₘ,tgₑ,shapes...,L,dL,mat)
end;

# Define now the residual function for the EulerBeam3D element.
vec3(v,ind) = SVector{3}(v[i] for i∈ind)


@espy function Muscade.residual(o::EulerBeam3D,   X,U,A,t,SP,dbg) 
    P,ND        = constants(X),length(X)
    X_          = motion{P}(X)
    gp_,ε_,vₛₘ_,rₛₘ_,vₗ₂_ = kinematics(o,X_,justinvoke)
    ☼ε          = motion⁻¹{P,ND}(ε_  )  
    rₛₘ          = motion⁻¹{P,ND}(rₛₘ_)  
    ♢κ          = motion⁻¹{P,ND}(vₗ₂_).*(2/o.L) 
    vᵢ₀         =              (SVector(0,0,0),                                                                           )
    vᵢ₁         = ND<1 ? vᵢ₀ : (vᵢ₀...        , spin⁻¹(∂0(rₛₘ)' ∘₁ ∂1(rₛₘ))                                                 ) 
    vᵢ          = ND<2 ? vᵢ₁ : (vᵢ₁...                                   ,   spin⁻¹(∂1(rₛₘ)' ∘₁ ∂1(rₛₘ) + ∂0(rₛₘ)' ∘₁ ∂2(rₛₘ)))  

    X₀          = ∂0(X)
    TX₀         = revariate{1}(X₀)
    Tgp,Tε,Tvₛₘ,Trₛₘ,Tvₗ₂ = kinematics(o,TX₀,fast)
    ε∂X₀        = composeJacobian{P}(Tε ,X₀)
    vₛₘ∂X₀       = composeJacobian{P}(Tvₛₘ,X₀)

    vₗ₂ = motion⁻¹{P,ND}(vₗ₂_).*(2/o.L)
    gp          = ntuple(ngp) do igp
        x_,κ_   = gp_[igp].x, gp_[igp].κ   
        ☼x      = motion⁻¹{P,ND}(x_  )
        ☼κ      = motion⁻¹{P,ND}(κ_  )

        Tx,Tκ   = Tgp[igp].x, Tgp[igp].κ
        x∂X₀    = composeJacobian{P}(Tx,X₀)
        κ∂X₀    = composeJacobian{P}(Tκ,X₀)
        fᵢ,mᵢ,fₑ,mₑ = ☼resultants(o.mat,ε,κ,x,rₛₘ,vᵢ)          # call the "resultant" function to compute loads (local coordinates) from strains/curvatures/etc. using material properties. Note that output is dual of input. 
        R       = (fᵢ ∘₀ ε∂X₀ + mᵢ ∘₁ κ∂X₀ + fₑ ∘₁ x∂X₀ + mₑ ∘₁ vₛₘ∂X₀) * o.dL[igp]     # Contribution to the local nodal load of this Gauss point  [ndof] = scalar*[ndof] + [ndim]⋅[ndim,ndof] + [ndim]⋅[ndim,ndof]
        @named(R)
    end
    R           = sum(gpᵢ.R for gpᵢ∈gp) 
    return R,noFB  
end;

# Transformation to corotated system and interpolation
function kinematics(o::EulerBeam3D,X₀,fast)  
    cₘ,rₘ,tgₘ,tgₑ,ζnod,ζgp,L  = o.cₘ,o.rₘ,o.tgₘ,o.tgₑ,o.ζnod,o.ζgp,o.L   # As-meshed element coordinates and describing tangential vector

    ## transformation to corotated system
    uᵧ₁,vᵧ₁,uᵧ₂,vᵧ₂  = vec3(X₀,1:3), vec3(X₀,4:6), vec3(X₀,7:9), vec3(X₀,10:12)
    vₗ₂,rₛₘ,vₛₘ        = fast(SVector(vᵧ₁...,vᵧ₂...)) do v
        vᵧ₁,vᵧ₂      = vec3(v,1:3), vec3(v,4:6)
        rₛ₁          = fast(Rodrigues,vᵧ₁)
        rₛ₂          = fast(Rodrigues,vᵧ₂)
        vₗ₂          = 0.5*Rodrigues⁻¹(rₛ₂ ∘₁ rₛ₁')
        rₛₘ          = fast(Rodrigues,vₗ₂) ∘₁ rₛ₁ ∘₁ o.rₘ  
        vₛₘ          = Rodrigues⁻¹(rₛₘ)              
        return vₗ₂,rₛₘ,vₛₘ
    end   
    cₛ               = 0.5*(uᵧ₁+uᵧ₂)
    uₗ₂              = rₛₘ'∘₁(uᵧ₂+tgₘ*ζnod[2]-cₛ)-tgₑ*ζnod[2]    #Local displacement of node 2
    ## interpolation
    ε               = √((uₗ₂[1]+L/2)^2+uₗ₂[2]^2+uₗ₂[3]^2)*2/L - 1.      
    gp              = ntuple(ngp) do igp
        yₐ,yᵤ,yᵥ,κₐ,κᵤ,κᵥ = o.yₐ[igp],o.yᵤ[igp],o.yᵥ[igp],o.κₐ[igp],o.κᵤ[igp],o.κᵥ[igp]
        κ           = SVector(         κₐ*vₗ₂[1], κᵤ*uₗ₂[2]+κᵥ*vₗ₂[3], κᵤ*uₗ₂[3]-κᵥ*vₗ₂[2])  
        y           = SVector(yₐ*uₗ₂[1]         , yᵤ*uₗ₂[2]+yᵥ*vₗ₂[3], yᵤ*uₗ₂[3]-yᵥ*vₗ₂[2])  
        x           = rₛₘ∘₁(tgₑ*ζgp[igp]+y)+cₛ+cₘ 
        (κ=κ,x=x)
    end
    return gp,ε,vₛₘ,rₛₘ,vₗ₂
end

