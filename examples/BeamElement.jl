include("Rotations.jl")

# # Euler beam element

using StaticArrays, LinearAlgebra, Muscade

# Data structure containing the cross section material properties
struct BeamCrossSection
    EA  :: 𝕣  # axial stiffness 
    EI₂ :: 𝕣 # bending stiffness about second axis
    EI₃ :: 𝕣 # bending stiffness about third axis
    GJ  :: 𝕣 # torsional stiffness (about longitudinal axis)
    μ   :: 𝕣 # mass per unit length
    ι₁  :: 𝕣 # (mass) moment of inertia for rotation about the element longitudinal axis per unit length
end
BeamCrossSection(;EA=EA,EI₂=EI₂,EI₃=EI₃,GJ=GJ,μ=μ,ι₁=ι₁) = BeamCrossSection(EA,EI₂,EI₃,GJ,μ,ι₁);

# Resultant function that computes the internal loads from the strains and curvatures, and external loads on the element. 
@espy function resultants(o::BeamCrossSection,ε,κ,xᵧ,rₛₘ,vᵢ) 
    r₀  = ∂0(rₛₘ)  # orientation of the element's local refsys
    vᵢ₁ = ∂1(vᵢ)  # intrinsic rotation rate         of the element's local refsys
    vᵢ₂ = ∂2(vᵢ)  # intrinsic rotation acceleration of the element's local refsys
    xᵧ₀,xᵧ₁,xᵧ₂ = ∂0(xᵧ),∂1(xᵧ),∂2(xᵧ)
    xₗ₁          = xᵧ₁ ∘₁ r₀
    xₗ₂          = xᵧ₂ ∘₁ r₀
    ## Compute drag force (example) and added-mass force (example)
    ## fa = ρ * Ca .* xₗ₂
    ## fd = .5 * ρ * A .* Cd .* xₗ₁ #.* abs.(xₗ₁)
    ## Compute translational inertia force 
    fi = o.μ * xᵧ₂ 
    ☼fₑ = fi # external forces at Gauss point.
    ## Compute roll inertia moment 
    m₁ₗ = o.ι₁*vᵢ₂[1] #local 
    mᵧ = ∂0(rₛₘ)[:,1] * m₁ₗ #global
    ☼mₑ = mᵧ  # external couples at Gauss point. 
    ## Compute internal loads
    ☼fᵢ = o.EA*∂0(ε)
    ## WARNING: curvatures are defined as rate of rotation along the element, not second derivatives of deflection.  
    ## Hence κ[3]>0 implies +2 direction is inside curve, 
    ##       κ[2]>0 implies -3 direction is inside curve.
    ☼mᵢ  = SVector(o.GJ*∂0(κ)[1],o.EI₃*∂0(κ)[2],o.EI₂*∂0(κ)[3])
    return fᵢ,mᵢ,fₑ,mₑ
end;

## Static Euler beam element, with two nodes, two Gauss points and 12 degrees of freedom. 
const ngp        = 4
const ndim       = 3
const nXdof      = 12
const nUdof      = 3
const nXnod      = 2;

# Shape functions for a beam element with support ζ∈[-1/2,1/2]. Though the shape function matrices are sparse, do not "unroll" them.  That would be applyer but considerably clutter the code                          
yₐ(ζ) =            2ζ       # differential axial displacement or roll field
yᵤ(ζ) =  -4ζ^3    +3ζ       # deflection due to differential nodal transverse translation
yᵥ(ζ) =        ζ^2   - 1/4  # deflection due to differenttial rotation (bending, not torsion)
κₐ(ζ) =                2    # torsion  . κₐ = yₐ′ . Divide by L .    
κᵤ(ζ) =  -24ζ               # curvature. κᵤ = yᵤ′′. Divide by L².
κᵥ(ζ) =                2;   # curvature. κᵥ = yᵥ′′. Divide by L .

# Data structure describing an EulerBeam3D element as meshed
struct EulerBeam3D{Mat,Uforce} <: AbstractElement
    cₘ       :: SVector{3,𝕣}     # Position of the middle of the element, as meshed
    rₘ       :: Mat33{𝕣}         # Orientation of the element, as meshed, represented by a rotation matrix (from global to local)
    ζgp      :: SVector{ngp,𝕣}   # Location of the Gauss points for the normalized element with length 1
    ζnod     :: SVector{nXnod,𝕣} # Location of the nodes for the normalized element with length 1
    tgₘ      :: SVector{ndim,𝕣}  # Vector connecting the nodes of the element in the global coordinate system
    tgₑ      :: SVector{ndim,𝕣}  # Vector connecting the nodes of the element in the local coordinate system
    yₐ       :: SVector{ngp,𝕣}   # Value at gp of shape function for differential axial displacement or roll field
    yᵤ       :: SVector{ngp,𝕣}   # Value at gp of shape function for deflection due to differential nodal transverse translation
    yᵥ       :: SVector{ngp,𝕣}   # Value at gp of shape function for deflection due to differenttial rotation (bending, not torsion)
    κₐ       :: SVector{ngp,𝕣}   # Value at gp of shape function for torsion  . κₐ = yₐ′ . Divided by L .    
    κᵤ       :: SVector{ngp,𝕣}   # Value at gp of shape function for curvature. κᵤ = yᵤ′′. Divided by L².
    κᵥ       :: SVector{ngp,𝕣}   # Value at gp of shape function for curvature. κᵥ = yᵥ′′. Divided by L .
    L        :: 𝕣                # as meshed length of the element
    dL       :: SVector{ngp,𝕣}   # length associated to each Gauss point
    mat      :: Mat              # used to store material properties (BeamCrossSection, for example)
end;

# For performance, `residual` will only accept differentiation to first order
Muscade.nosecondorder(::Type{<:EulerBeam3D}) = Val(true)

# Define nodes, classes, and field names of dofs
Muscade.doflist(     ::Type{EulerBeam3D{Mat,false}}) where{Mat} = 
        (inod = (1,1,1,1,1,1, 2,2,2,2,2,2), 
         class= ntuple(i->:X,nXdof), 
         field= (:t1,:t2,:t3,:r1,:r2,:r3, :t1,:t2,:t3,:r1,:r2,:r3) )
Muscade.doflist(     ::Type{EulerBeam3D{Mat,true}}) where{Mat} = 
        (inod = (1,1,1,1,1,1, 2,2,2,2,2,2, 3,3,3), 
         class= (ntuple(i->:X,nXdof)...,ntuple(i->:U,nUdof)...), 
         field= (:t1,:t2,:t3,:r1,:r2,:r3, :t1,:t2,:t3,:r1,:r2,:r3,  :t1,:t2,:t3) )
# ElementType for the EulerBeam3D element. Arguments: node list, material, and direction of the first bending axis in the global coordinate system.  
EulerBeam3D(nod;kwargs...) = EulerBeam3D{false}(nod;kwargs...) # by default, EulerBeam3D does not have Udof.
function EulerBeam3D{Udof}(nod::Vector{Node};mat,orient2::SVector{ndim,𝕣}=SVector(0.,1.,0.)) where {Udof}
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
    dL    = SVector{ngp}(L/2*(18-sqrt(30))/36,L/2*(18+sqrt(30))/36  ,L/2*(18+sqrt(30))/36,L/2*(18-sqrt(30))/36  ) 
    ## Location ζgp of the Gauss points for a unit-length beam element, with nodes at ζnod=±1/2. 
    ζgp     = SVector{ngp  }(-1/2*sqrt(3/7+2/7*sqrt(6/5)),-1/2*sqrt(3/7-2/7*sqrt(6/5)), +1/2*sqrt(3/7-2/7*sqrt(6/5)),+1/2*sqrt(3/7+2/7*sqrt(6/5))) 
    ζnod    = SVector{nXnod}(-1/2  ,1/2  )
    shapes  = (yₐ.(ζgp), yᵤ.(ζgp), yᵥ.(ζgp)*L, κₐ.(ζgp)/L, κᵤ.(ζgp)/L^2, κᵥ.(ζgp)/L)
    return EulerBeam3D{typeof(mat),Udof}(cₘ,rₘ,ζgp,ζnod,tgₘ,tgₑ,shapes...,L,dL,mat)
end;

# Define now the residual function for the EulerBeam3D element.
@espy function Muscade.residual(o::EulerBeam3D{Mat,Udof},   X,U,A,t,SP,dbg) where{Mat,Udof}
    P,ND                = constants(X),length(X)
    ## Compute all quantities at Gauss point, their time derivatives, including intrinsic roll rate and acceleration
    gp_,ε_,vₛₘ_,rₛₘ_,vₗ₂_,_,_ = kinematics{:direct}(o,motion{P}(X))
#    gp_,ε_,vₛₘ_,rₛₘ_,vₗ₂_,_,_ = kinematics(o,motion{P}(X),(f,x)->f(x))
    gpval,☼ε , rₛₘ       = motion⁻¹{P,ND}(gp_,ε_,rₛₘ_  ) 
    vᵢ                  = intrinsicrotationrates(rₛₘ)
    ## compute all Jacobians of the above quantities with respect to X₀
    X₀                  = ∂0(X)
    TX₀                 = revariate{1}(X₀)  # check type
    Tgp,Tε,Tvₛₘ,_,_,_,_  = kinematics{:compose}(o,TX₀) # the crux
    gp∂X₀,ε∂X₀,vₛₘ∂X₀    = composeJacobian{P}((Tgp,Tε,Tvₛₘ),X₀)
    ## Quadrature loop: compute resultants
    gp                  = ntuple(ngp) do igp
        ☼x,☼κ           = gpval[igp].x, gpval[igp].κ   
        x∂X₀,κ∂X₀       = gp∂X₀[igp].x, gp∂X₀[igp].κ
        fᵢ,mᵢ,fₑ,mₑ     = ☼resultants(o.mat,ε,κ,x,rₛₘ,vᵢ)          # call the "resultant" function to compute loads (local coordinates) from strains/curvatures/etc. using material properties. Note that output is dual of input. 
        fₑ              = Udof ? fₑ-∂0(U) : fₑ                    # U is per unit length
        R_               = (fᵢ ∘₀ ε∂X₀ + mᵢ ∘₁ κ∂X₀ + fₑ ∘₁ x∂X₀ + mₑ ∘₁ vₛₘ∂X₀) * o.dL[igp]     # Contribution to the local nodal load of this Gauss point  [nXdof] = scalar*[nXdof] + [ndim]⋅[ndim,nXdof] + [ndim]⋅[ndim,nXdof]
        @named(R_)
    end
    R                   = sum(gpᵢ.R_ for gpᵢ∈gp) 
    ♢κ                  = motion⁻¹{P,ND}(SVector(vₗ₂_[1],vₗ₂_[3],-vₗ₂_[2])).*(2/o.L) 
    ♢rₛₘ                 = motion⁻¹{P,ND}(rₛₘ_)
    return R,noFB  
end;
struct kinematics{Mode} end
function kinematics{Mode}(o::EulerBeam3D,X₀)  where{Mode}
    cₘ,rₘ,tgₘ,tgₑ,ζnod,ζgp,L  = o.cₘ,o.rₘ,o.tgₘ,o.tgₑ,o.ζnod,o.ζgp,o.L   # As-meshed element coordinates and describing tangential vector
    vₛₘ,rₛₘ,uₗ₂,vₗ₂,cₛₘ  = corotated{Mode}(o,X₀)
    ε                = √((uₗ₂[1]+L/2)^2+uₗ₂[2]^2+uₗ₂[3]^2)*2/L - 1.      
    gp               = ntuple(ngp) do igp  # gp[igp].κ, gp[igp].x
        yₐ,yᵤ,yᵥ,κₐ,κᵤ,κᵥ = o.yₐ[igp],o.yᵤ[igp],o.yᵥ[igp],o.κₐ[igp],o.κᵤ[igp],o.κᵥ[igp]
        κ            = SVector(         κₐ*vₗ₂[1], κᵤ*uₗ₂[2]+κᵥ*vₗ₂[3], κᵤ*uₗ₂[3]-κᵥ*vₗ₂[2])  
        y            = SVector(yₐ*uₗ₂[1]         , yᵤ*uₗ₂[2]+yᵥ*vₗ₂[3], yᵤ*uₗ₂[3]-yᵥ*vₗ₂[2])  
        x            = rₛₘ∘₁(tgₑ*ζgp[igp]+y)+cₛₘ 
        (κ=κ,x=x)  
    end
    return gp,ε,vₛₘ,rₛₘ,vₗ₂,uₗ₂,cₛₘ
end

vec3(v,ind) = SVector{3}(v[i] for i∈ind);
struct corotated{Mode} end 
function corotated{Mode}(o::EulerBeam3D,X₀)  where{Mode}
    cₘ,rₘ,tgₘ,tgₑ,ζnod,ζgp,L  = o.cₘ,o.rₘ,o.tgₘ,o.tgₑ,o.ζnod,o.ζgp,o.L   # As-meshed element coordinates and describing tangential vector
    uᵧ₁,uᵧ₂,vᵧ             = vec3(X₀,1:3), vec3(X₀,7:9), SVector(X₀[4],X₀[5],X₀[6],X₀[10],X₀[11],X₀[12])
    Δvᵧ,rₛₘ,vₛₘ             = apply{Mode}(vᵧ) do v
        vᵧ₁,vᵧ₂            = vec3(v,1:3), vec3(v,4:6)
        rₛ₁                = apply{Mode}(Rodrigues,vᵧ₁)
        rₛ₂                = apply{Mode}(Rodrigues,vᵧ₂)
        Δvᵧ_         = 0.5*Rodrigues⁻¹(rₛ₂ ∘₁ rₛ₁')
        rₛₘ_          = apply{Mode}(Rodrigues,Δvᵧ_) ∘₁ rₛ₁ ∘₁ o.rₘ  
        vₛₘ_          = Rodrigues⁻¹(rₛₘ_)              
        return Δvᵧ_,rₛₘ_,vₛₘ_
    end   
    cₛ               = 0.5*(uᵧ₁+uᵧ₂)
    uₗ₂              = rₛₘ' ∘₁ (uᵧ₂+tgₘ*ζnod[2]-cₛ)-tgₑ*ζnod[2]    #Local displacement of node 2
    vₗ₂              = rₛₘ' ∘₁ Δvᵧ
    return vₛₘ,rₛₘ,uₗ₂,vₗ₂,cₛ+cₘ
end;

#Assumes that BeamElement.jl has been included previously, and that "using GLMakie" has been invoked
using GLMakie
"""

Drawing a `EulerBeam3D`.

    draw!(axis,state)

    draw!(axis,state;EulerBeam3D=(;style=:shape))

    α      = 2π*(0:19)/20
    circle = 0.1*[cos.(α) sin.(α)]'
    draw!(axis,state;EulerBeam3D=(;style=:solid,section = circle))

`style=:shape` shows the deformed neutral axis of the element. It has optional arguments `frame=true` 
(draws the element's corotated frame of reference)
and `nseg=10` (number of points to show the deflected shape of each element). 

`style=:solid` shows the deformed shape of the element. It requires the input `section=...` to be given
a matrix of size `(2,nsec)` describing `nsec` points around the cross section of the element (no need to close 
the circumference by repeating the first point at the end).  It has optional arguments `nseg=10` as above, `marking=true`
to draw a longitudinal marking and `solid_color=:yellow`.
 
Other optional arguments (and their default values) are
- `Udof` (`true` iff element has Udofs) wether to draw U-forces.
- `draw_frame = false` wether to draw the local reference frame of each element
- `draw_marking = true` wether to draw "longitudinal marking" along the element.  Will only draw if style=:solid.
- `nseg = 1` number of segments to display the shape of a deformed element
- `solid_color = :yellow` color of the surface if `style=:solid`
- `line_color = :black` color of the line if `style=:sshape`
- `Uscale = 1.` How many meter is a Newton per meter?
"""
function Muscade.allocate_drawing(axis,o::AbstractVector{EulerBeam3D{Tmat,Udof}};kwargs...) where{Tmat,Udof}
    args                 = default{:EulerBeam3D     }(kwargs,(;)     )  
    section              = default{:section         }(args,zeros(2,0))  
    nsec                 = size(section,2)                            
    opt = (default(args,(style=:shape,draw_frame=false,draw_marking=true,nseg=1,
                  solid_color=:yellow,line_color=:black,Uscale=1.,Udof=Udof))...,
            nel          = length(o)                                  ,
            nsec         = nsec                                       ,                    
            section      = section                                    ,
            markrad      = nsec==0 ? 0. : 1.01*maximum(section[1,:])      
        )
    opt.style==:solid && nsec<2 && muscadeerror("An section description must be provided for 'solid' plot")
    nel_shape         = opt.style==:shape ? opt.nel   : 0
    nel_shape_frame   = opt.draw_frame    ? nel_shape : 0
    nel_solid         = opt.style==:solid ? opt.nel   : 0 
    nel_solid_marking = opt.draw_marking  ? nel_solid : 0
    nel_udof          = opt.Udof          ? opt.nel   : 0

    mut=(
            node         = 𝕣2(undef,3,3*opt.nel)                        ,
            shape_x      = 𝕣2(undef,3,(opt.nseg+2)*nel_shape)           ,   
            shape_frame  = 𝕣2(undef,3,3*3*nel_shape_frame)              , # idim, point-point-lift, ivec, iel
            solid_vertex = 𝕣2(undef,3,opt.nsec*(opt.nseg+1)*nel_solid)  , 
            solid_face   = 𝕫2(undef,2*opt.nsec* opt.nseg   *nel_solid,3),
            solid_mark   = 𝕣2(undef,3,(opt.nseg+2)*nel_solid_marking)   ,     
            ucrest       = 𝕣2(undef,3,5*nel_udof)                       , # idim, 6point-lift,iel
        )   
    return mut,opt
end

function Muscade.update_drawing(axis,o::AbstractVector{EulerBeam3D{Tmat,Udof}},oldmut,opt, Λ,X,U,A,t,SP,dbg) where{Tmat,Udof} 
    mut               = oldmut 
    X₀                = ∂0(X)
    U₀                = ∂0(U)
    it1,ir1,it2,ir2   = SVector{3}(1:3),SVector{3}(4:6),SVector{3}(7:9),SVector{3}(10:12)
    nsec              = size(opt.section,2) 
    node = reshape(mut.node,(3,3,opt.nel))
    for (iel,oᵢ) = enumerate(o)
        node[:,1,iel] = oᵢ.cₘ - oᵢ.tgₘ/2 + X₀[it1,iel]
        node[:,2,iel] = oᵢ.cₘ + oᵢ.tgₘ/2 + X₀[it2,iel]
        node[:,3,iel].= NaN  
    end

    if opt.style==:shape
        ζ = range(-1/2,1/2,opt.nseg+1)
        if opt.draw_frame shape_frame  = reshape(mut.shape_frame ,(3,3,3       ,opt.nel)) end
        if opt.Udof       ucrest       = reshape(mut.ucrest,      (3,5         ,opt.nel)) end
        shape_x                        = reshape(mut.shape_x     ,(3,opt.nseg+2,opt.nel))
        for (iel,oᵢ) = enumerate(o)
            cₘ,rₘ,tgₘ,tgₑ,ζnod,ζgp,L  = oᵢ.cₘ,oᵢ.rₘ,oᵢ.tgₘ,oᵢ.tgₑ,oᵢ.ζnod,oᵢ.ζgp,oᵢ.L   
            X₀ₑ = view(X₀,:,iel)
            vₛₘ,rₛₘ,uₗ₂,vₗ₂,cₛₘ = corotated{:direct}(oᵢ,X₀ₑ) 
            if opt.draw_frame
                for ivec = 1:3
                    shape_frame[:,1,ivec,iel] = cₛₘ
                    shape_frame[:,2,ivec,iel] = cₛₘ + oᵢ.L/3*rₛₘ[:,ivec]
                    shape_frame[:,3,ivec,iel].= NaN
                end
            end
            if opt.Udof
                ucrest[:,1,iel] = node[:,1,iel]
                ucrest[:,2,iel] = node[:,1,iel] + rₛₘ ∘₁ view(U₀,:,iel) * opt.Uscale
                ucrest[:,3,iel] = node[:,2,iel] + rₛₘ ∘₁ view(U₀,:,iel) * opt.Uscale
                ucrest[:,4,iel] = node[:,2,iel]
                ucrest[:,5,iel].= NaN
            end
            for (i,ζᵢ) ∈ enumerate(ζ)
                y          = SVector(yₐ(ζᵢ)*uₗ₂[1] , yᵤ(ζᵢ)*uₗ₂[2]+L*yᵥ(ζᵢ)*vₗ₂[3], yᵤ(ζᵢ)*uₗ₂[3]-L*yᵥ(ζᵢ)*vₗ₂[2])  
                shape_x[:,i         ,iel] = rₛₘ∘₁(tgₑ*ζᵢ+y)+cₛₘ 
                shape_x[:,opt.nseg+2,iel].= NaN
            end        
        end
    elseif opt.style==:solid
        ζ = range(-1/2,1/2,opt.nseg+1)
        idx(iel,iseg,isec) = mod_onebased(isec,opt.nsec)+opt.nsec*(iseg-1+(opt.nseg+1)*(iel-1)) # 1st index into rvertex
        if opt.Udof         ucrest         = reshape(mut.ucrest       ,(3,5          ,opt.nel)) end
        if opt.draw_marking solid_mark     = reshape(mut.solid_mark  ,(3,opt.nseg+2 ,opt.nel)) end
        solid_face                         = reshape(mut.solid_face  ,(2,opt.nsec, opt.nseg   ,opt.nel,3))
        solid_vertex                       = reshape(mut.solid_vertex,(3,opt.nsec, opt.nseg+1 ,opt.nel))
        for (iel,oᵢ) = enumerate(o)
            cₘ,rₘ,tgₘ,tgₑ,ζnod,ζgp,L  = oᵢ.cₘ,oᵢ.rₘ,oᵢ.tgₘ,oᵢ.tgₑ,oᵢ.ζnod,oᵢ.ζgp,oᵢ.L   
            X₀ₑ = view(X₀,:,iel)
            vₛₘ,rₛₘ,uₗ₂,vₗ₂,cₛₘ = corotated{:direct}(oᵢ,X₀ₑ) 
            vᵧ₁,vᵧ₂          = vec3(X₀ₑ,4:6), vec3(X₀ₑ,10:12)
            rₛ₁              = Rodrigues(vᵧ₁)
            rₛ₂              = Rodrigues(vᵧ₂)
            if opt.Udof
                ucrest[:,1,iel] = node[:,1,iel]
                ucrest[:,2,iel] = node[:,1,iel] + rₛₘ ∘₁ view(U₀,:,iel) * opt.Uscale
                ucrest[:,3,iel] = node[:,2,iel] + rₛₘ ∘₁ view(U₀,:,iel) * opt.Uscale
                ucrest[:,4,iel] = node[:,2,iel]
                ucrest[:,5,iel].= NaN
            end
            Δv = Rodrigues⁻¹(rₛ₂ ∘₁ rₛ₁')/opt.nseg
            for (iseg,ζᵢ) ∈ enumerate(ζ) # actualy iterating over nseg+1 segment boundaries
                y  = SVector(yₐ(ζᵢ)*uₗ₂[1] , yᵤ(ζᵢ)*uₗ₂[2]+L*yᵥ(ζᵢ)*vₗ₂[3], yᵤ(ζᵢ)*uₗ₂[3]-L*yᵥ(ζᵢ)*vₗ₂[2])  # interpolate
                xn = rₛₘ∘₁(tgₑ*ζᵢ+y)+cₛₘ # point on neutral axis
                r  = Rodrigues((iseg-1)*Δv) ∘₁ rₛ₁ ∘₁ rₘ  
                if opt.draw_marking 
                    solid_mark[:,    iseg  ,iel] = xn .+ r[:,2]*opt.markrad 
                    solid_mark[:,opt.nseg+2,iel].= NaN 
                end
                for isec = 1:opt.nsec
                    solid_vertex[:,isec,iseg,iel] = xn .+ r[:,2]*opt.section[1,isec] + r[:,3]*opt.section[2,isec] 
                    if iseg≤opt.nseg
                        i1,i2,i3,i4 = idx(iel,iseg,isec),idx(iel,iseg  ,isec+1),idx(iel,iseg+1,isec  ),idx(iel,iseg+1,isec+1)
                        solid_face[1,isec,iseg,iel,:] = SVector(i1,i2,i4)    
                        solid_face[2,isec,iseg,iel,:] = SVector(i1,i4,i3)   
                    end
                end
            end  
        end
    end
    return mut
end

function Muscade.display_drawing!(axis,::Type{EulerBeam3D{Tmat,Udof}},obs,opt) where{Tmat,Udof}
    scatter!(                                          axis, obs.node                         ,color = opt.line_color , marker=:circle,markersize=3)  
    opt.style==:shape  &&                     lines!(  axis, obs.shape_x                      ,color = opt.line_color ,linewidth=.5                )
    opt.style==:shape  && opt.draw_frame   && lines!(  axis, obs.shape_frame                  ,color = :grey          ,linewidth=.5                )    
    opt.style==:solid  &&                     mesh!(   axis, obs.solid_vertex, obs.solid_face ,color = opt.solid_color                             )  
    opt.style==:solid  && opt.draw_marking && lines!(  axis, obs.solid_mark                   ,color = opt.line_color                              )    
    opt.Udof           &&                     lines!(  axis, obs.ucrest                       ,color = :red           ,linewidth=.5                )    
end



