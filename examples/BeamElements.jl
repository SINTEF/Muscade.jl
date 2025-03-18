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
@espy function resultants(o::BeamCrossSection,ε,κ,xᵧ,rot,::Val{P},::Val{ND}) where{P,ND}
    ε₀          = Muscade.position{P,ND}(ε) 
    κ₀          = Muscade.position{P,ND}(κ) 
    xᵧ₀,xᵧ₁,xᵧ₂ = Muscade.posVelAcc{P,ND}(xᵧ)
    xₗ₁ = xᵧ₁ ∘ rot
    xₗ₂ = xᵧ₂ ∘ rot
    ρ = 1025.0
    Cd = SVector(0.0,1.0,1.0)
    A  = SVector(0.0,1.0,1.0)
    fd = .5 * ρ .* A .* Cd .* xₗ₁ .* abs.(xₗ₁) #mind the sign: forces exerted by element on its environment
    μ   = 1.0
    fi = μ * xₗ₂ 
    Ca = SVector(0.0,1.0,1.0)
    fa = ρ * Ca .* xₗ₂
    ☼f₁ = o.EA*ε₀ # replace by ε₀
    ☼m  = SVector(o.GJ*κ₀[1],o.EI*κ₀[2],o.EI*κ₀[3])# replace by κ₀ 
    ☼fₑ = fd+fi+fa # SVector(0.,0.,0.) # external forces at Gauss point (no external moment/torque/... so far). fₑ is in local coordinates # add inertia and drag
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
    Nδx      :: SVector{ngp,SMatrix{ndim,ndof,𝕣,ndim*ndof}} # coordinates of the Gauss points
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
    Nδx      = SVector{ngp}(@SMatrix [Nₐ₁(ζᵢ)   0          0          0          0           0          Nₐ₂(ζᵢ)   0          0           0           0           0         ;
                                     0         Nᵤ₁(ζᵢ)    0          0          0           Nᵥ₁(ζᵢ)    0         Nᵤ₂(ζᵢ)    0           0           0           Nᵥ₂(ζᵢ)   ;
                                     0         0          Nᵤ₁(ζᵢ)    0          -Nᵥ₁(ζᵢ)    0          0         0          Nᵤ₂(ζᵢ)     0           -Nᵥ₂(ζᵢ)    0         ] for ζᵢ∈ζgp) # Nδx[igp][idim,idof]
    return EulerBeam3D(cₘ,rₘ,ζgp,ζnod,tgₘ,tgₑ,Nε,Nκ,Nδx,dL,mat)
end

const saco = StaticArrays.sacollect
const v3   = SVector{3};

# Define now the residual function for the EulerBeam3D element.
@espy function Muscade.residual(o::EulerBeam3D,   X,U,A,t,SP,dbg) 
    ## Fetch the element properties 
    cₘ,rₘ,tgₘ,tgₑ     = o.cₘ,o.rₘ,o.tgₘ,o.tgₑ   # As-meshed element coordinates and describing tangential vector
    Nε,Nκ,Nδx         = o.Nε,o.Nκ,o.Nδx           # From shape functions
    ζgp,ζnod,dL      = o.ζgp,o.ζnod,o.dL        # Gauss points coordinates, node coordinates and length associated to each Gauss point
    ## In the following, the goal is to compute the Jacobian T transforming quantities from/to local/global coordinate systems using automatic differentiation
    P                = constants(X)
    ND               = length(X)
    X_               = Muscade.motion{P}(X)
    ## δX_l and T contain time derivatives, cₛ,rₛₘ do not
    δXₗ,T,cₛ,rₛₘ = coordinateTransform(X_,o)
    ## Compute local load contributions at each Gauss point
    gp              = ntuple(ngp) do igp
        ☼ε,☼κ,☼δxₗ   = Nε[igp]∘δXₗ, Nκ[igp]∘δXₗ, Nδx[igp]∘δXₗ   # axial strain, curvatures, displacement - all local (including their time derivatives)
        ☼x          = rₛₘ∘(tgₑ*ζgp[igp]+δxₗ)+cₛ+cₘ             # [ndim], global coordinates of Gauss points
        f₁,m,fₑ     = ☼resultants(o.mat,ε,κ,x,rₛₘ,Val(P),Val(ND))          # call the "resultant" function to compute loads (local coordinates) from strains/curvatures/etc. using material properties. Note that output is dual of input. 
        Rₗ           = (f₁ ∘₀ Nε[igp] + m∘Nκ[igp] + fₑ∘Nδx[igp])*dL[igp]     # Contribution to the local nodal load of this Gauss point  [ndof] = scalar*[ndof] + [ndim]⋅[ndim,ndof] + [ndim]⋅[ndim,ndof]
        @named(Rₗ)
    end
    ## Summation of local load contributions from each Gauss point, and transformation to the global coordinate system. 
    R  = sum(gpᵢ.Rₗ for gpᵢ∈gp) ∘ T 
    return R,noFB  
end


function coordinateTransform(X,o::EulerBeam3D)
    ## P is an integer that enables variate to keep track with respect to what X,U,A,t have been differentated before. Note that P is defined at compilation time. No run time. 
    P                = constants(X) 
    ## We are going do differentiate wrt X (to get the Jacobian T for example). 
    ## Describe here the content of ΔX contains (zeros and ones)
    ΔX               = variate{P,ndof}(X)
    ## Note that X is a tuple containing (positions, velocities, accelerations) and ∂0(X) returns only positions
    ## X is not an adiff with respect to time. Use motions to go from tuple to adiff. Do not forget constants, ses motions doc. 
    ## Fetch the nodal displacements uᵧ₁ uᵧ₂ and rotations vᵧ₁, vᵧ₂ from X, expressed in the global coordinate system
    uᵧ₁,vᵧ₁,uᵧ₂,vᵧ₂  = SVector{3}(ΔX[i] for i∈1:3), SVector{3}(ΔX[i] for i∈4:6),SVector{3}(ΔX[i] for i∈7:9),SVector{3}(ΔX[i] for i∈10:12)
    ## Conversion to the local coordinate system
    cₛ               = (uᵧ₁+uᵧ₂)/2
    rₛ               = Rodrigues((vᵧ₁+vᵧ₂)/2)
    rₛ               = Rodrigues(adjust(rₛ∘o.tgₘ,o.tgₘ+uᵧ₂-uᵧ₁))∘rₛ   
    rₛₘ              = rₛ∘o.rₘ
    uₗ₁              = rₛₘ'∘(uᵧ₁+o.tgₘ*o.ζnod[1]-cₛ)-o.tgₑ*o.ζnod[1]    #Local displacement of node 1
    uₗ₂              = rₛₘ'∘(uᵧ₂+o.tgₘ*o.ζnod[2]-cₛ)-o.tgₑ*o.ζnod[2]    #Local displacement of node 2
    vₗ₁              = Rodrigues⁻¹(rₛₘ'∘Rodrigues(vᵧ₁)∘o.rₘ)      #Local rotation of node 1
    vₗ₂              = Rodrigues⁻¹(rₛₘ'∘Rodrigues(vᵧ₂)∘o.rₘ)      #Local rotation of node 2
    ## δXₗ contains all local displacements ("value") and partial derivatives ("δ") with respect to ΔX
    δXₗ,T            = value_∂{P,ndof}(SVector(uₗ₁...,vₗ₁...,uₗ₂...,vₗ₂...))
    cₛ               = value{P}(cₛ)
    rₛₘ              = value{P}(rₛₘ)
    return δXₗ,T,cₛ,rₛₘ 
end