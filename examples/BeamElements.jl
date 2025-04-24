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
@espy function resultants(o::BeamCrossSection,x,ε,κ,rot) 
    ☼f₁ = o.EA*∂0(ε)
    ☼m  = SVector(o.GJ*∂0(κ)[1],o.EI*∂0(κ)[2],o.EI*∂0(κ)[3]) 
    μ   =  SVector(1.0,1.0,1.0)
    xᵧ₂ = ∂2(x)
    # @show VALUE(xᵧ₂)
    xₗ₂ = xᵧ₂ ∘₁ rot
    # @show VALUE(xₗ₂)
    fi = -μ .* xₗ₂ 
    # @show VALUE(fi)
    ☼fₑ = -fi #SVector(0.,0.,0.) # external forces at Gauss point (no external moment/torque/... so far). fₑ is in local coordinates # add inertia and drag
    return f₁,m,fₑ
end;

# Static Euler beam element, with two nodes, two Gauss points and 12 degrees of freedom. 
const ngp        = 2
const ndim       = 3
const ndof       = 12
const nnod       = 2;

# Shape functions for a beam element with support ζ∈[-1/2,1/2]. Though the shape function matrices are sparse, do not "unroll" them.  That would be faster but considerably clutter the code                          
# Nₐᵢ is the displacement field induced by node i's axial translations (also the rotation induced by node i's rotation about the element's direction/torsion)
Nₐ₁(ζ) =                    -ζ +1/2                          
Nₐ₂(ζ) =                     ζ +1/2;                 
# Nᵤᵢ is the deflection field induced by node i's transverse displacements
Nᵤ₁(ζ) =  2ζ^3          -3/2*ζ +1/2              
Nᵤ₂(ζ) = -2ζ^3          +3/2*ζ +1/2;          
# Nᵥᵢ is the deflection field induced by node i's rotations (bending, not torsion)
Nᵥ₁(ζ) =   ζ^3 -1/2*ζ^2 -1/4*ζ +1/8          
Nᵥ₂(ζ) =   ζ^3 +1/2*ζ^2 -1/4*ζ -1/8;          

# First derivatives ∂N/∂ζ used to compute strain and torsion. For an element of length L, use ∂N/∂x=∂N/∂ζ/L.
Bₐ₁(ζ) = -1        
Bₐ₂(ζ) =  1;
# Second derivatives ∂²N/∂ζ² used to compute curvature. For an element of length L, use ∂²N/∂x²=∂²N/∂ζ²/L².
Bᵤ₁(ζ) =   12ζ
Bᵥ₁(ζ) =    6ζ-1
Bᵤ₂(ζ) =  -12ζ  
Bᵥ₂(ζ) =    6ζ+1;

# Data structure describing an EulerBeam3D element as meshed
struct EulerBeam3D{Mat} <: AbstractElement
    cₘ       :: SVector{3,𝕣}    # Position of the middle of the element
    rₘ       :: Mat33{𝕣}        # Orientation of the element (see code)
    ζgp      :: SVector{ngp,𝕣}  # Location of the Gauss points for the normalized element with length 1
    ζnod     :: SVector{nnod,𝕣} # Location of the nodes for the normalized element with length 1
    tgₘ      :: SVector{ndim,𝕣} # Vector connecting the nodes of the element in the global coordinate system
    tgₑ      :: SVector{ndim,𝕣} # Vector connecting the nodes of the element in the local coordinate system
    Nε       :: SVector{ngp,SVector{     ndof,𝕣}}           # strain at the Gauss points
    Nκ       :: SVector{ngp,SMatrix{ndim,ndof,𝕣,ndim*ndof}} # curvatures at the Gauss points
    Ny      :: SVector{ngp,SMatrix{ndim,ndof,𝕣,ndim*ndof}} # coordinates of the Gauss points
    dL       :: SVector{ngp,𝕣}  # length associated to each Gauss point
    mat      :: Mat # Used to store material properties (BeamCrossSection, for example)
end

# Define nodes, classes, and field names for Muscade
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
    ζnod    = SVector{ngp }(-1/2  ,1/2  ) # ζ∈[-1/2,1/2]
    L²      = L^2
    ## Using the first derivative of the shape function to get the strain at Gauss points 
    Nε      = SVector{ngp}(@SVector [Bₐ₁(ζᵢ)/L,0,         0,         0,         0,          0,         Bₐ₂(ζᵢ)/L,0,         0,          0,          0,          0         ] for ζᵢ∈ζgp)  # Nε[igp][idof]
    ## Using the first and second derivatives of the shape function to get the torsion and curvature at Gauss points
    Nκ      = SVector{ngp}(@SMatrix [0         0          0          Bₐ₁(ζᵢ)/L  0           0          0         0          0           Bₐ₂(ζᵢ)/L   0           0         ;
                                     0         Bᵤ₁(ζᵢ)/L² 0          0          0           Bᵥ₁(ζᵢ)/L 0         Bᵤ₂(ζᵢ)/L² 0           0           0           Bᵥ₂(ζᵢ)/L;
                                     0         0          Bᵤ₁(ζᵢ)/L² 0          -Bᵥ₁(ζᵢ)/L 0          0         0          Bᵤ₂(ζᵢ)/L²  0           -Bᵥ₂(ζᵢ)/L  0         ] for ζᵢ∈ζgp) # Nκ[igp][idim,idof]
    ## Using the shape functions to get the coordinates of the Gauss points
    Ny      = SVector{ngp}(@SMatrix [Nₐ₁(ζᵢ)   0          0          0          0           0          Nₐ₂(ζᵢ)   0          0           0           0           0         ;
                                     0         Nᵤ₁(ζᵢ)    0          0          0           Nᵥ₁(ζᵢ)    0         Nᵤ₂(ζᵢ)    0           0           0           Nᵥ₂(ζᵢ)   ;
                                     0         0          Nᵤ₁(ζᵢ)    0          -Nᵥ₁(ζᵢ)    0          0         0          Nᵤ₂(ζᵢ)     0           -Nᵥ₂(ζᵢ)    0         ] for ζᵢ∈ζgp) # Ny[igp][idim,idof]
    return EulerBeam3D(cₘ,rₘ,ζgp,ζnod,tgₘ,tgₑ,Nε,Nκ,Ny,dL,mat)
end

const saco = StaticArrays.sacollect
const v3   = SVector{3};

# Define now the residual function for the EulerBeam3D element.

@espy function Muscade.residual(o::EulerBeam3D,   X,U,A,t,SP,dbg) 
    cₘ,tgₑ          = o.cₘ,o.tgₑ   # As-meshed element coordinates and describing tangential vector
    Nε,Nκ,Ny        = o.Nε,o.Nκ,o.Ny           # From shape functions
    ζgp,dL          = o.ζgp,o.dL        # Gauss points coordinates, node coordinates and length associated to each Gauss point

    X₀              = ∂0(X)
    # OD              = min(2,precedence(X₀)+1) 
    # TY₀,_,_         = Taylor{OD}(X->X₀2Y₀(o,X),X₀)
    # Y₀∂X₀           = ∂(TY₀)(X₀)

    P,ND            = constants(X,U,A,t),length(X) 
    X_              = motion{P}(X)
    Y₀_,rₛₘ_,cₛ_     = X₀2Y₀(o::EulerBeam3D,X_)
    _,  rₛₘ ,_      = X₀2Y₀(o::EulerBeam3D,X₀)
    unpack          = motion⁻¹{P,ND  }
    gp              = ntuple(ngp) do igp
        ☼ε          = motion⁻¹{P,ND  }(                   Nε[igp]∘₁Y₀_       ) # TODO type unstable because P,ND are no longer compile-time consts in closure?
        ☼κ          = motion⁻¹{P,ND  }(                   Nκ[igp]∘₁Y₀_       ) 
        ☼x          = motion⁻¹{P,ND  }(rₛₘ_∘₁(tgₑ*ζgp[igp]+Ny[igp]∘₁Y₀_)+cₛ_+cₘ) 
        ☼rₛₘ         = motion⁻¹{P,ND,0}(rₛₘ_                                   )
        f₁,m,fₑ     = ☼resultants(o.mat,x,ε,κ,rₛₘ)          # call the "resultant" function to compute loads (local coordinates) from strains/curvatures/etc. using material properties. Note that output is dual of input. 
        Rₗ           = (f₁ ∘₀ Nε[igp] + m∘₁Nκ[igp] + fₑ∘₁Ny[igp]) * dL[igp]     # Contribution to the local nodal load of this Gauss point  [ndof] = scalar*[ndof] + [ndim]⋅[ndim,ndof] + [ndim]⋅[ndim,ndof]
        # Rₗ           = (fₑ∘₁Ny[igp]) * dL[igp]     # Contribution to the local nodal load of this Gauss point  [ndof] = scalar*[ndof] + [ndim]⋅[ndim,ndof] + [ndim]⋅[ndim,ndof]
        # @show VALUE(fₑ∘₁Ny[igp])
        @named(Rₗ)
    end
    # R               = sum(gpᵢ.Rₗ for gpᵢ∈gp) ∘₁ Y₀∂X₀ 
    Rₗ               = sum(gpᵢ.Rₗ for gpᵢ∈gp) 
    r₁,r₂,r₃,r₄  = SVector{3}(Rₗ[i] for i∈1:3), SVector{3}(Rₗ[i] for i∈4:6),SVector{3}(Rₗ[i] for i∈7:9),SVector{3}(Rₗ[i] for i∈10:12)
    R₁,R₂,R₃,R₄  = r₁∘₁rₛₘ,r₂∘₁rₛₘ,r₃∘₁rₛₘ,r₄∘₁rₛₘ
    R            = SVector(R₁...,R₂...,R₃...,R₄...)           
    # @show VALUE(gp[1].Rₗ + gp[2].Rₗ)
    # @show VALUE(Y₀∂X₀)
    @show VALUE(R)
    return R,noFB  
end
function X₀2Y₀(o::EulerBeam3D,X₀)  
    rₘ,tgₘ,tgₑ,ζnod  = o.rₘ,o.tgₘ,o.tgₑ,o.ζnod   # As-meshed element coordinates and describing tangential vector
    uᵧ₁,vᵧ₁,uᵧ₂,vᵧ₂  = SVector{3}(X₀[i] for i∈1:3), SVector{3}(X₀[i] for i∈4:6),SVector{3}(X₀[i] for i∈7:9),SVector{3}(X₀[i] for i∈10:12)
    rₛ               = Rodrigues((vᵧ₁+vᵧ₂)/2)
    rₛ               = Rodrigues(adjust(rₛ∘₁tgₘ,tgₘ+uᵧ₂-uᵧ₁))∘₁rₛ   
    rₛₘ              = rₛ∘₁o.rₘ
    cₛ               = (uᵧ₁+uᵧ₂)/2
    uₗ₁              = rₛₘ'∘₁(uᵧ₁+tgₘ*ζnod[1]-cₛ)-tgₑ*ζnod[1]    #Local displacement of node 1
    uₗ₂              = rₛₘ'∘₁(uᵧ₂+tgₘ*ζnod[2]-cₛ)-tgₑ*ζnod[2]    #Local displacement of node 2
    vₗ₁              = Rodrigues⁻¹(rₛₘ'∘₁Rodrigues(vᵧ₁)∘₁rₘ)      #Local rotation of node 1
    vₗ₂              = Rodrigues⁻¹(rₛₘ'∘₁Rodrigues(vᵧ₂)∘₁rₘ)      #Local rotation of node 2
    Y₀              = SVector(uₗ₁...,vₗ₁...,uₗ₂...,vₗ₂...)           #  Y₀ , Y₀∂X₀ = ∂Y / ∂X
    return Y₀,rₛₘ,cₛ
end

