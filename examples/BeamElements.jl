include("Rotations.jl")

# # Euler beam element

using StaticArrays, LinearAlgebra
using Muscade

# Data structure containing the cross section material properties
struct BeamCrossSection
    EA :: 𝕣
    EI :: 𝕣
    GJ :: 𝕣
end
    # ρ  :: 𝕣 
    # μ  :: 𝕣 
    # Rotation about x.... moment of inertia?
    # Cd :: SVector{3,𝕣}
    # Ca :: SVector{3,𝕣}
    # A  :: SVector{3,𝕣}

BeamCrossSection(;EA=EA,EI=EI,GJ=GJ) = BeamCrossSection(EA,EI,GJ);

# Resultant function that computes the internal loads from the strains and curvatures, and external loads on the element. 
@espy function resultants(o::BeamCrossSection,ε,κ,x,vₛₘ) 
    # WARNING: curvatures are defined as rate of rotation along the element, not second derivatives of deflection.  
    # Hence κ[3]>0 implies +2 direction is inside curve, 
    #       κ[2]>0 implies -3 direction is inside curve.
    ☼f₁ = o.EA*∂0(ε)
    ☼m  = SVector(o.GJ*∂0(κ)[1],o.EI*∂0(κ)[2],o.EI*∂0(κ)[3])# replace by κ₀ 
    ☼fₑ = SVector(0.,0.,0.) # external forces  at Gauss point. fₑ is in local coordinates # add inertia and drag
    ☼mₑ = SVector(0.,0.,0.) # external couples at Gauss point. mₑ is in local coordinates 
    return f₁,m,fₑ,mₑ
end;



## Static Euler beam element, with two nodes, two Gauss points and 12 degrees of freedom. 

const ngp        = 2
const ndim       = 3
const ndof       = 12
const nnod       = 2;

# Shape functions for a beam element with support ζ∈[-1/2,1/2]. Though the shape function matrices are sparse, do not "unroll" them.  That would be faster but considerably clutter the code                          
yₐ(ζ) =            2ζ       # differential axial displacement or roll field
yᵤ(ζ) =  -4ζ^3    +3ζ       # deflection due to differential nodal transverse translation
yᵥ(ζ) =        ζ^2   - 1/4  # deflection due to differenttial rotation (bending, not torsion)
κₐ(ζ) =                2    # torsion  . κₐ = yₐ′ . Divide by L .    
κᵤ(ζ) =  -24ζ               # curvature. κᵤ = yᵤ′′. Divide by L².
κᵥ(ζ) =                2    # curvature. κᵥ = yᵥ′′. Divide by L .

# Data structure describing an EulerBeam3D element as meshed
struct EulerBeam3D{Mat} <: AbstractElement
    cₘ       :: SVector{3,𝕣}    # Position of the middle of the element
    rₘ       :: Mat33{𝕣}        # Orientation of the element (see code)
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
    mat      :: Mat # Used to store material properties (BeamCrossSection, for example)
end

# Define nodes, classes, and field names of dofs
Muscade.doflist(::Type{<:EulerBeam3D}) = (inod = (1,1,1,1,1,1, 2,2,2,2,2,2), class= ntuple(i->:X,ndof), field= (:t1,:t2,:t3,:r1,:r2,:r3, :t1,:t2,:t3,:r1,:r2,:r3) )

# Define now the constructor for the EulerBeam3D element. Arguments: node coordinates and direction of the first bending axis in the global coordinate system.  
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
    ## Length associated to each Gauss point
    dL      = SVector{ngp }(L/2   , L/2 )
    ## Location of the Gauss points for a unit-length beam element, with nodes at ±1/2. 
    ζgp     = SVector{ngp }(-1/2√3,1/2√3) # ζ∈[-1/2,1/2]
    ζnod    = SVector{nnod}(-1/2  ,1/2  ) # ζ∈[-1/2,1/2]
    shapes  = (yₐ.(ζgp), yᵤ.(ζgp), yᵥ.(ζgp), κₐ.(ζgp)/L, κᵤ.(ζgp)/L^2, κᵥ.(ζgp)/L)
    return EulerBeam3D(cₘ,rₘ,ζgp,ζnod,tgₘ,tgₑ,shapes...,L,dL,mat)
end

# Define now the residual function for the EulerBeam3D element.

@espy function Muscade.residual(o::EulerBeam3D,   X,U,A,t,SP,dbg) 
    X₀          = ∂0(X)
    TX₀         = revariate{1}(X₀)
    Tgp,Tε,Tvₛₘ  = kinematics(o,TX₀)

    P,ND        = precedence(X₀),length(X)
    X_          = motion{P}(X)

    ☼ε          = motion⁻¹{P,ND}(compose(value{P+1}( Tε  ),X_))
    ☼vₛₘ         = motion⁻¹{P,ND}(compose(value{P+1}( Tε  ),X_))
    ε∂X₀        =                compose(∂{P+1,ndof}(Tε  ),X₀ )
    vₛₘ∂X₀        =              compose(∂{P+1,ndof}(Tvₛₘ  ),X₀ )
    gp          = ntuple(ngp) do igp
        Tx,Tκ   = Tgp[igp].x, Tgp[igp].κ
        ☼x      = motion⁻¹{P,ND}(compose(value{P+1}( Tx  ),X_))
        ☼κ      = motion⁻¹{P,ND}(compose(value{P+1}( Tκ  ),X_))
        x∂X₀    =                compose(∂{P+1,ndof}(Tx  ),X₀ )
        κ∂X₀    =                compose(∂{P+1,ndof}(Tκ  ),X₀ )
        f₁,mᵢ,fₑ,mₑ = ☼resultants(o.mat,ε,κ,x,vₛₘ)          # call the "resultant" function to compute loads (local coordinates) from strains/curvatures/etc. using material properties. Note that output is dual of input. 
        R       = (f₁ ∘₀ ε∂X₀ + mᵢ ∘₁ κ∂X₀ + fₑ ∘₁ x∂X₀ + mₑ ∘₁ vₛₘ∂X₀) * o.dL[igp]     # Contribution to the local nodal load of this Gauss point  [ndof] = scalar*[ndof] + [ndim]⋅[ndim,ndof] + [ndim]⋅[ndim,ndof]
        @named(R)
    end
    R               = sum(gpᵢ.R for gpᵢ∈gp) 
    return R,noFB  
end
function kinematics(o::EulerBeam3D,X₀)  
    cₘ,rₘ,tgₘ,tgₑ,ζnod,ζgp,L  = o.cₘ,o.rₘ,o.tgₘ,o.tgₑ,o.ζnod,o.ζgp,o.L   # As-meshed element coordinates and describing tangential vector

    # transformation to corotated system
    uᵧ₁,vᵧ₁,uᵧ₂,vᵧ₂  = SVector{3}(X₀[i] for i∈1:3), SVector{3}(X₀[i] for i∈4:6),SVector{3}(X₀[i] for i∈7:9),SVector{3}(X₀[i] for i∈10:12)
    vₛ               = (vᵧ₁+vᵧ₂)/2
    rₛₘ              = Rodrigues(vₛ) ∘₁ o.rₘ
    vₛₘ              = Rodrigues⁻¹(rₛₘ)
    cₛ               = (uᵧ₁+uᵧ₂)/2
    uₗ₂              = rₛₘ'∘₁(uᵧ₂+tgₘ*ζnod[2]-cₛ)-tgₑ*ζnod[2]    #Local displacement of node 2
    vₗ₂              = Rodrigues⁻¹(rₛₘ'∘₁Rodrigues(vᵧ₂)∘₁rₘ)     #Local rotation of node 2
    
    # interpolation
    ε               = √((uₗ₂[1]+L/2)^2+uₗ₂[2]^2+uₗ₂[3]^2)*2/L - 1.       
    gp              = ntuple(ngp) do igp
        yₐ,yᵤ,yᵥ,κₐ,κᵤ,κᵥ = o.yₐ[igp],o.yᵤ[igp],o.yᵥ[igp],o.κₐ[igp],o.κᵤ[igp],o.κᵥ[igp]
        κ           = SVector(         κₐ*vₗ₂[1], κᵤ*uₗ₂[2]+κᵥ*vₗ₂[3], κᵤ*uₗ₂[3]-κᵥ*vₗ₂[2])  
        y           = SVector(yₐ*uₗ₂[1]         , yᵤ*uₗ₂[2]+yᵥ*vₗ₂[3], yᵤ*uₗ₂[3]-yᵥ*vₗ₂[2])                              
        x           = rₛₘ∘₁(tgₑ*ζgp[igp]+y)+cₛ+cₘ 
        (κ=κ,x=x)
    end
    return gp,ε,vₛₘ
end

