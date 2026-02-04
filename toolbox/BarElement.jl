# # Bar element

using StaticArrays, LinearAlgebra, Muscade

# Data structure containing the cross section material properties
struct AxisymmetricBarCrossSection
    EA  :: 𝕣 # Axial stiffness [N]
    μ   :: 𝕣 # Mass per unit length [kg/m]
    w   :: 𝕣 # Weight per unit length [N/m]
    Caₜ :: 𝕣 # Tangential added mass per unit length [kg/m]
    Clₜ :: 𝕣 # Tangential linear damping coefficient per unit length [N/m/(m/s)]
    Cqₜ :: 𝕣 # Tangential quadratic damping coefficient per unit length [N/m/(m/s)^2], for example from drag
    Caₙ :: 𝕣 # Normal added mass per unit length [kg/m] for motions along second axis
    Clₙ :: 𝕣 # Normal linear damping coefficient per unit length [N/m/(m/s)] for motions along second axis
    Cqₙ :: 𝕣 # Normal quadratic damping coefficient per unit length [N/m/(m/s)^2], for motions along second axis
    # TODO: add gravity field to bar properties (time dependent), and use it to compute the weight. This to enable static analyses. 
end
AxisymmetricBarCrossSection(;EA,μ=μ,w=0.,Caₜ=0.,Clₜ=0.,Cqₜ=0.,Caₙ=0.,Clₙ=0.,Cqₙ=0.) = AxisymmetricBarCrossSection(EA,μ,w,Caₜ,Clₜ,Cqₜ,Caₙ,Clₙ,Cqₙ);

const ngp        = 4 # Number of Gauss points
const ndim       = 3 # Number of dimensions
const nXdof      = 6 # Number of X-class degrees of freedom
const nXnod      = 2 # Number of X-class nodes
const nUdof      = 3 # Number of U-class degrees of freedom   

# Shape functions
ψ₁(ζ) = -ζ + 1/2          
ψ₂(ζ) =  ζ + 1/2         

# Data structure describing an Bar3D element as meshed
"""
    Bar3D

A bar element
"""
struct Bar3D{Mat,Uforce} <: AbstractElement
    cₘ       :: SVector{ndim,𝕣}  # Position of the middle of the element, as meshed
    tgₘ      :: SVector{ndim,𝕣}  # Vector connecting the nodes of the element in the global coordinate system (global)
    tgₑ      :: SVector{ndim,𝕣}  # Vector connecting the nodes of the element in the local coordinate system  (local)
    L₀       :: 𝕣                # As-meshed length of the element
    Lₛ        :: 𝕣                # Stress-free length of the element (by default equal to L₀)
    mat      :: Mat              # Used to store material properties (AxisymmetricBarCrossSection, for example)
    wgp      :: SVector{ngp,𝕣}   # Weight associated to each Gauss point
    ζgp      :: SVector{ngp,𝕣}   # Location of the Gauss points for the normalized element defined on [-1/2,1/2]
    ζnod     :: SVector{nXnod,𝕣} # Location of the nodes for the normalized element defined on [-1/2,1/2]
    ψ₁       :: SVector{ngp,𝕣}   # Value at gp of shape function
    ψ₂       :: SVector{ngp,𝕣}   # Value at gp of shape function
end;

# For performance, `residual` will only accept differentiation to first order
Muscade.no_second_order(::Type{<:Bar3D}) = Val(true)

# Define nodes, classes, and field names of dofs for the element, in absence/presence of U-dofs, respectively
Muscade.doflist(     ::Type{Bar3D{Mat,false}}) where{Mat} = 
        (inod = (1,1,1,         2,2,2), 
         class= (:X,:X,:X,      :X,:X,:X), 
         field= (:t1,:t2,:t3,   :t1,:t2,:t3) )
Muscade.doflist(     ::Type{Bar3D{Mat,true}}) where{Mat} = 
        (inod = (1,1,1,         2,2,2,          3,3,3), 
         class= (:X,:X,:X,      :X,:X,:X,       :U,:U,:U),  
         field= (:t1,:t2,:t3,   :t1,:t2,:t3,    :t1,:t2,:t3) )

# Constructor of the Bar3D element. 
Bar3D(nod;kwargs...) = Bar3D{false}(nod;kwargs...) # by default, Bar3D does not have Udof.
function Bar3D{Udof}(nod::Vector{Node};mat,Lₛ=0.) where {Udof}
    c       = coord(nod)
    # Position of the middle of the element in the global coordinate system (as-meshed)
    cₘ      = SVector{3}((c[1]+c[2])/2)
    # Tangential vector to the element in the local and global coordinate system, and its length (as-meshed)
    tgₘ     = SVector{ndim}(c[2]-c[1])
    L₀ =  norm(tgₘ)
    if Lₛ == 0.
        Lₛ =  L₀
    end
    tgₑ     = SVector{ndim}(L₀,0,0)
    # Location ζgp of the Gauss points associated weigths, and values of the shape functions, for a unit-length bar element, with nodes at ζnod=±1/2. 
    wgp    = SVector{ngp}(      L₀/2*(18-sqrt(30))/36,          L₀/2*(18+sqrt(30))/36  ,        L₀/2*(18+sqrt(30))/36,          L₀/2*(18-sqrt(30))/36       ) 
    ζgp     = SVector{ngp  }(   -1/2*sqrt(3/7+2/7*sqrt(6/5)),   -1/2*sqrt(3/7-2/7*sqrt(6/5)),   +1/2*sqrt(3/7-2/7*sqrt(6/5)),   +1/2*sqrt(3/7+2/7*sqrt(6/5))) 
    ζnod    = SVector{nXnod}(   -1/2  ,1/2  )
    shapes  = (ψ₁.(ζgp), ψ₂.(ζgp))
    return Bar3D{typeof(mat),Udof}(cₘ,tgₘ,tgₑ,L₀,Lₛ,mat,wgp,ζgp,ζnod,shapes...)
end;

# Internal and external loads at a given Gauss point with coordinates x, and strain ε. 
@espy function resultants(o::AxisymmetricBarCrossSection,ε,x,u) 
    # Unit vector tangential to the element
    δ      = ∂0(u)
    # Velocity and acceleration
    v,a      = ∂1(x),∂2(x)    
    # Inertia force
    fi      = o.μ * a
    # Weight
    fw =  SVector(0,0,o.w)
    # Added mass
    aₜ = a ∘₁ δ         # Tangential acceleration (scalar)
    aₙ = a - aₜ * δ     # Normal acceleration (vector)
    fa  = SVector{3}(o.Caₜ * aₜ,o.Caₙ * aₙ[2],o.Caₙ* aₙ[3])
    # Linear and quadratic damping
    vₜ = v ∘₁ δ         # Tangential acceleration (scalar)
    vₙ = v - vₜ * δ     # Normal acceleration (vector)
    fqₜ = o.Cqₜ * vₜ^2;          if vₜ < 0; fqₜ = -fqₜ end
    fqₙ2 = o.Cqₙ * vₙ[2]^2;     if vₙ[2] < 0; fqₙ2 = -fqₙ2 end
    fqₙ3 = o.Cqₙ * vₙ[3]^2;     if vₙ[3] < 0; fqₙ3 = -fqₙ3 end
    fd  = SVector{3}(o.Clₜ * vₜ + fqₜ, o.Clₙ * vₙ[2] + fqₙ2 ,o.Clₙ* vₙ[3]+ fqₙ3 )
    # Sum of external forces
    ☼fe      =   fi+fw+fa+fd
    # Internal forces
    ☼fᵢ      = o.EA*∂0(ε)
    return fᵢ,fe
end;

# The function below is already defined for beam elements
# vec3(v,ind) = SVector{3}(v[i] for i∈ind);

# Define now the residual function for the Bar3D element.
@espy function Muscade.residual(o::Bar3D{Mat,Udof},   X,U,A,t,SP,dbg) where{Mat,Udof}
    # Obtain motions (i.e. including velocity and accelerations) from X
    P,ND    = constants(X),length(X)
    x_      = motion{P}(X)
    # Motions of the nodes, center of the element
    uᵧ₁,uᵧ₂   = vec3(x_,1:3), vec3(x_,4:6) 
    c        = o.cₘ + 0.5*(uᵧ₁+uᵧ₂) 
    # Element direction and length
    tg      = o.tgₘ + uᵧ₂ - uᵧ₁
    L       = √((o.L₀+uᵧ₂[1]-uᵧ₁[1])^2+(uᵧ₂[2]-uᵧ₁[2])^2+(uᵧ₂[3]-uᵧ₁[3])^2)
    δ_       = tg/L
    # Strains
    ε_       = L/o.Lₛ - 1
    # Compute how strains vary with nodal displacements (will be used in the Princple of Virtual Work, PVW)
    ε,δ = motion⁻¹{P,ND}(ε_,δ_); δ₀ = ∂0(δ)
    ε∂X₀ = 1/o.L₀*SVector{6}(-δ₀[1],-δ₀[2],-δ₀[3],δ₀[1],δ₀[2],δ₀[3]) 
    # Compute Gauss point kinematics
    gp = ntuple(ngp) do igp; 
        x = c + tg * o.ζgp[igp]; 
        @named(x); 
    end
    # Compute loads at Gauss points
    gpContrib = ntuple(ngp) do igp
        ζ = o.ζgp[igp]                          # Coordinate of the Gauss point along [-1/2,1/2]
        # Compute how motions of Gauss point vary with nodal displacements (used in PVW below)
        x∂X₀ = SMatrix{3,6}(ψ₁(ζ),0,0, 0,ψ₁(ζ),0, 0,0,ψ₁(ζ), ψ₂(ζ),0,0, 0,ψ₂(ζ),0, 0,0,ψ₂(ζ))   
        x = motion⁻¹{P,ND}(gp[igp].x)          # Physical location of the Gauss point 
        fᵢ,fₑ     = ☼resultants(o.mat,ε,x,δ)   # Compute loads from strains/motions, etc.
        fₑ        = Udof ? fₑ-∂0(U) : fₑ       # If there are unknown loads, they're added here (U is per unit length)
        #  Application of PVW, local contribution of the integral over the element
        R_        = ( fᵢ ∘₀ ε∂X₀ + fₑ ∘₁ x∂X₀ ) * o.wgp[igp]   
        @named(R_);
    end
    R                   = sum(gpᵢ.R_ for gpᵢ∈gpContrib)
return R,noFB  
end;

# The following functions explain how the bar element should be drawn
using GLMakie
"""

Drawing a `Bar3D`.

    draw!(axis,state)
 
Optional arguments (and their default values) are
- `line_color = :black` color of the line
- `Udof` (`true` iff element has Udofs) wether to draw U-forces.
- `Uscale = 1.` How many meter is a Newton per meter?
"""
function Muscade.allocate_drawing(axis,o::AbstractVector{Bar3D{Tmat,Udof}};kwargs...) where{Tmat,Udof}
    args                 = default{:Bar3D     }(kwargs,(;)     )  
    section              = default{:section         }(args,zeros(2,0))  
    nsec                 = size(section,2)                            
    opt = (default(args,(line_color=:black,Uscale=1.,Udof=Udof))...,nel= length(o))
    nel_udof          = opt.Udof          ? opt.nel   : 0
    mut=(
            node         = 𝕣2(undef,3,3*opt.nel)                        ,
            shape_x      = 𝕣2(undef,3,3*opt.nel)           ,   
            ucrest       = 𝕣2(undef,3,5*nel_udof)                       , # idim, 6point-lift,iel
        )   
    return mut,opt
end

function Muscade.update_drawing(axis,o::AbstractVector{Bar3D{Tmat,Udof}},oldmut,opt, Λ,X,U,A,t,SP,dbg) where{Tmat,Udof} 
    mut               = oldmut 
    X₀                = ∂0(X) # Nodal displacements
    U₀                = ∂0(U) # External forces
    it1,it2           = SVector{3}(1:3),SVector{3}(4:6)
    node              = reshape(mut.node,       (3,3,opt.nel))
    shape_x           = reshape(mut.shape_x,    (3,3,opt.nel))
    if opt.Udof       
        ucrest        = reshape(mut.ucrest,      (3,5         ,opt.nel)) 
    end
    for (iel,oᵢ) = enumerate(o)
        node[:,1,iel] = oᵢ.cₘ - oᵢ.tgₘ/2 + X₀[it1,iel]
        node[:,2,iel] = oᵢ.cₘ + oᵢ.tgₘ/2 + X₀[it2,iel]
        node[:,3,iel].= NaN  
        shape_x[:,1,iel] = oᵢ.cₘ - oᵢ.tgₘ/2 + X₀[it1,iel]
        shape_x[:,2,iel] = oᵢ.cₘ + oᵢ.tgₘ/2 + X₀[it2,iel]
        shape_x[:,3,iel].= NaN  
        if opt.Udof
        ucrest[:,1,iel] = node[:,1,iel]
        ucrest[:,2,iel] = node[:,1,iel] +  view(U₀,:,iel) * opt.Uscale
        ucrest[:,3,iel] = node[:,2,iel] +  view(U₀,:,iel) * opt.Uscale
        ucrest[:,4,iel] = node[:,2,iel]
        ucrest[:,5,iel].= NaN
        end
    end
    return mut
end

function Muscade.display_drawing!(axis,::Type{Bar3D{Tmat,Udof}},obs,opt) where{Tmat,Udof}
    scatter!(           axis, obs.node       ,color = opt.line_color , marker=:circle,markersize=3)  
    lines!(             axis, obs.shape_x    ,color = opt.line_color ,linewidth=.5                )
    opt.Udof  && lines!(axis, obs.ucrest     ,color = :red           ,linewidth=.5                )    
end



