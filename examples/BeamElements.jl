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
    L        :: 𝕣               # as meshed length of the element
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
    dL    = SVector{ngp}(L/2*(18-sqrt(30))/36,L/2*(18+sqrt(30))/36  ,L/2*(18+sqrt(30))/36,L/2*(18-sqrt(30))/36  ) 
    ## Location ζgp of the Gauss points for a unit-length beam element, with nodes at ζnod=±1/2. 
    ζgp     = SVector{ngp }(-1/2*sqrt(3/7+2/7*sqrt(6/5)),-1/2*sqrt(3/7-2/7*sqrt(6/5)), +1/2*sqrt(3/7-2/7*sqrt(6/5)),+1/2*sqrt(3/7+2/7*sqrt(6/5))) 
    ζnod    = SVector{nnod }(-1/2  ,1/2  )
    shapes  = (yₐ.(ζgp), yᵤ.(ζgp), yᵥ.(ζgp)*L, κₐ.(ζgp)/L, κᵤ.(ζgp)/L^2, κᵥ.(ζgp)/L)
    return EulerBeam3D(cₘ,rₘ,ζgp,ζnod,tgₘ,tgₑ,shapes...,L,dL,mat)
end;

# Define now the residual function for the EulerBeam3D element.
vec3(v,ind) = SVector{3}(v[i] for i∈ind)

# Il semble que la perfection soit atteinte non quand il n’y a plus rien à ajouter, mais quand il n’y a plus rien à retrancher. Antoine de Saint-Exupéry.
@espy function Muscade.residual(o::EulerBeam3D,   X,U,A,t,SP,dbg) 
    P,ND                = constants(X),length(X)
    ## Compute all quantities at Gauss point, their time derivatives, including intrinsic roll rate and acceleration
    gp_,ε_,vₛₘ_,rₛₘ_,vₗ₂_,_ = kinematics(o,motion{P}(X))
    gpval,☼ε , rₛₘ       = motion⁻¹{P,ND}(gp_,ε_,rₛₘ_  ) 
    vᵢ                  = intrinsicrotationrates(rₛₘ)
    ## compute all Jacobians of the above quantities with respect to X₀
    X₀                  = ∂0(X)
    TX₀                 = revariate{1}(X₀)
    Tgp,Tε,Tvₛₘ,_,_,_    = kinematics(o,TX₀,fast)
    gp∂X₀,ε∂X₀,vₛₘ∂X₀    = composeJacobian{P}((Tgp,Tε,Tvₛₘ),X₀)
    ## Quadrature loop: compute resultants, and 
    gp                  = ntuple(ngp) do igp
        ☼x,☼κ           = gpval[igp].x, gpval[igp].κ   
        x∂X₀,κ∂X₀       = gp∂X₀[igp].x, gp∂X₀[igp].κ
        fᵢ,mᵢ,fₑ,mₑ     = ☼resultants(o.mat,ε,κ,x,rₛₘ,vᵢ)          # call the "resultant" function to compute loads (local coordinates) from strains/curvatures/etc. using material properties. Note that output is dual of input. 
        R               = (fᵢ ∘₀ ε∂X₀ + mᵢ ∘₁ κ∂X₀ + fₑ ∘₁ x∂X₀ + mₑ ∘₁ vₛₘ∂X₀) * o.dL[igp]     # Contribution to the local nodal load of this Gauss point  [ndof] = scalar*[ndof] + [ndim]⋅[ndim,ndof] + [ndim]⋅[ndim,ndof]
        @named(R)
    end
    R                   = sum(gpᵢ.R for gpᵢ∈gp) 
    ♢κ                  = motion⁻¹{P,ND}(vₗ₂_).*(2/o.L) 
    return R,noFB  
end;
function kinematics(o::EulerBeam3D,X₀,fast=justinvoke)  
    cₘ,rₘ,tgₘ,tgₑ,ζnod,ζgp,L  = o.cₘ,o.rₘ,o.tgₘ,o.tgₑ,o.ζnod,o.ζgp,o.L   # As-meshed element coordinates and describing tangential vector
    vₛₘ,rₛₘ,uₗ₂,vₗ₂,cₛₘ  = corotated(o,X₀,fast)
    ε                = √((uₗ₂[1]+L/2)^2+uₗ₂[2]^2+uₗ₂[3]^2)*2/L - 1.      
    gp               = ntuple(ngp) do igp  # gp[igp].κ, gp[igp].x
        yₐ,yᵤ,yᵥ,κₐ,κᵤ,κᵥ = o.yₐ[igp],o.yᵤ[igp],o.yᵥ[igp],o.κₐ[igp],o.κᵤ[igp],o.κᵥ[igp]
        κ            = SVector(         κₐ*vₗ₂[1], κᵤ*uₗ₂[2]+κᵥ*vₗ₂[3], κᵤ*uₗ₂[3]-κᵥ*vₗ₂[2])  
        y            = SVector(yₐ*uₗ₂[1]         , yᵤ*uₗ₂[2]+yᵥ*vₗ₂[3], yᵤ*uₗ₂[3]-yᵥ*vₗ₂[2])  
        x            = rₛₘ∘₁(tgₑ*ζgp[igp]+y)+cₛₘ 
        (κ=κ,x=x)  
    end
    return gp,ε,vₛₘ,rₛₘ,vₗ₂,uₗ₂
end
function corotated(o::EulerBeam3D,X₀,fast=justinvoke)  
    cₘ,rₘ,tgₘ,tgₑ,ζnod,ζgp,L  = o.cₘ,o.rₘ,o.tgₘ,o.tgₑ,o.ζnod,o.ζgp,o.L   # As-meshed element coordinates and describing tangential vector
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
    cₛ               = 0.5*(uᵧ₁+uᵧ₂)
    uₗ₂              = rₛₘ'∘₁(uᵧ₂+tgₘ*ζnod[2]-cₛ)-tgₑ*ζnod[2]    #Local displacement of node 2
    return vₛₘ,rₛₘ,uₗ₂,vₗ₂,cₛ+cₘ
end
"""

Drawing a `EulerBeam3D`.

    draw(axe,state)

    draw(axe,state;EulerBeam3D=(;style=:simple))

    draw(axe,state;EulerBeam3D=(;style=:shape))

    α = 2π*(0:19)/20
    circle = 0.1*[cos.(α) sin.(α)]'
    draw(axe,state;EulerBeam3D=(;style=:solid,section = circle))

`style=:simple` (default) shows a straight line between visible nodes.

`style=:shape` shows the deformed neutral axis of the element. It has optional arguments `frame=true` 
(draws the element's corotated frame of reference)
and `nseg=10` (number of points to show the deflected shape of each element). 

`style=:solid` shows the deformed shape of the element. It requires the input `section=...` to be given
a matrix of size `(2,nsec)` describing `nsec` points around the cross section of the element (no need to close 
the circumference by repeating the first point at the end).  It has optional arguments `nseg=10` as above, `marking=true`
to draw a longitudinal marking and `solid_color=:yellow`.
 
All above options share the optional argument `line_color=:black`.

"""
function Muscade.draw(axe,o::Vector{T}, Λ,X,U,A,t,SP,dbg;kwargs...) where{T<:EulerBeam3D}
    nel           = length(o)
    args          = default{:EulerBeam3D}(kwargs,(;))
    style         = default{:style   }(args,:simple)
    draw_frame    = default{:frame   }(args,true  )
    draw_marking  = default{:marking }(args,true  )
    nseg          = default{:nseg    }(args,10     )
    section       = default{:section }(args,zeros(2,0))
    solid_color   = default{:color   }(args,:yellow)
    line_color    = default{:color   }(args,:black)
    nsec          = size(section,2)
    X₀            = ∂0(X)
    it1,ir1,it2,ir2 = SVector{3}(1:3),SVector{3}(4:6),SVector{3}(7:9),SVector{3}(10:12)
    if     style==:simple
        line = Array{𝕣,3}(undef,3,3,nel)
        for (iel,oᵢ) = enumerate(o)
            line[:,1,iel] = oᵢ.cₘ - oᵢ.tgₘ/2 + X₀[it1,iel]
            line[:,2,iel] = oᵢ.cₘ + oᵢ.tgₘ/2 + X₀[it2,iel]
            line[:,3,iel].= NaN
        end
        rline = reshape(line,(3,3nel))
        lines!(  axe,rline,color = line_color                )    
        scatter!(axe,rline,color = line_color, marker=:circle)    
    elseif style==:shape
        ζ = range(-1/2,1/2,nseg+1)
        x = Array{𝕣,3}(undef,3,nseg+2,nel)
        if draw_frame  
            frame = 𝕣2(undef,3,9nel)
        end    
        for (iel,oᵢ) = enumerate(o)
            cₘ,rₘ,tgₘ,tgₑ,ζnod,ζgp,L  = oᵢ.cₘ,oᵢ.rₘ,oᵢ.tgₘ,oᵢ.tgₑ,oᵢ.ζnod,oᵢ.ζgp,oᵢ.L   
            X₀ₑ = view(X₀,:,iel)
            vₛₘ,rₛₘ,uₗ₂,vₗ₂,cₛₘ = corotated(oᵢ,X₀ₑ) 
            if draw_frame
                ℓ = oᵢ.L/3
                for i = 1:3
                    frame[:,9*(iel-1)+3*i-2] = cₛₘ
                    frame[:,9*(iel-1)+3*i-1] = cₛₘ + ℓ*rₛₘ[:,i]
                    frame[:,9*(iel-1)+3*i-0].= NaN
                end
            end
            for (i,ζᵢ) ∈ enumerate(ζ)
                y                       = SVector(yₐ(ζᵢ)*uₗ₂[1] , yᵤ(ζᵢ)*uₗ₂[2]+L*yᵥ(ζᵢ)*vₗ₂[3], yᵤ(ζᵢ)*uₗ₂[3]-L*yᵥ(ζᵢ)*vₗ₂[2])  
                x[:,i,iel] = rₛₘ∘₁(tgₑ*ζᵢ+y)+cₛₘ 
            end        
            x[:,nseg+2,iel] .= NaN
        end
        lines!(  axe,reshape(x,(3,(nseg+2)*nel)),color = line_color)
        xnod = x[:,[1,nseg+1],:]
        scatter!(axe,reshape(xnod,(3,2*nel)),color = line_color, marker=:circle)    
        if draw_frame  # move to "draw shape"
            lines!(axe,frame,color = :grey,linewidth=1)    
        end
    elseif style==:solid
        nsec≥2 || muscadeerror("An section description must be provided for 'solid' plot")
        ζ = range(-1/2,1/2,nseg+1)
        vertex             = Array{𝕣,4}(undef,3,  nsec, nseg+1 ,nel  ) 
        face               = Array{𝕫,5}(undef,  2,nsec, nseg   ,nel,3) 
        rvertex            = reshape(vertex,(3,   nsec*(nseg+1)*nel  ))
        rface              = reshape(face,  (   2*nsec* nseg   *nel,3))
        idx(iel,iseg,isec) = imod(isec,nsec)+nsec*(iseg-1+(nseg+1)*(iel-1)) # 1st index into rvertex
        if draw_marking
            mark   = Array{𝕣,3}(undef,3, nseg+2 ,nel  )     
            rmark  = reshape(mark,(3,   (nseg+2)*nel  ))
            markrad = 1.01*maximum(section[1,:])
        end
        for (iel,oᵢ) = enumerate(o)
            cₘ,rₘ,tgₘ,tgₑ,ζnod,ζgp,L  = oᵢ.cₘ,oᵢ.rₘ,oᵢ.tgₘ,oᵢ.tgₑ,oᵢ.ζnod,oᵢ.ζgp,oᵢ.L   
            X₀ₑ = view(X₀,:,iel)
            vₛₘ,rₛₘ,uₗ₂,vₗ₂,cₛₘ = corotated(oᵢ,X₀ₑ) 
            vᵧ₁,vᵧ₂          = vec3(X₀ₑ,4:6), vec3(X₀ₑ,10:12)
            rₛ₁              = Rodrigues(vᵧ₁)
            rₛ₂              = Rodrigues(vᵧ₂)
            for (iseg,ζᵢ) ∈ enumerate(ζ)
                y  = SVector(yₐ(ζᵢ)*uₗ₂[1] , yᵤ(ζᵢ)*uₗ₂[2]+L*yᵥ(ζᵢ)*vₗ₂[3], yᵤ(ζᵢ)*uₗ₂[3]-L*yᵥ(ζᵢ)*vₗ₂[2])  
                xn = rₛₘ∘₁(tgₑ*ζᵢ+y)+cₛₘ # point on neutral axis
                v  = (iseg-1)/nseg*Rodrigues⁻¹(rₛ₂ ∘₁ rₛ₁')
                r  = Rodrigues(v) ∘₁ rₛ₁ ∘₁ rₘ  
                if draw_marking 
                    mark[:,iseg,iel] = xn .+ r[:,2]*markrad 
                end
                for isec = 1:nsec
                    vertex[:,isec,iseg,iel] = xn .+ r[:,2]*section[1,isec] + r[:,3]*section[2,isec] 
                    if iseg≤nseg
                        i1,i2,i3,i4 = idx(iel,iseg,isec),idx(iel,iseg  ,isec+1),idx(iel,iseg+1,isec  ),idx(iel,iseg+1,isec+1)
                        face[1,isec,iseg,iel,:] = SVector(i1,i2,i4)    
                        face[2,isec,iseg,iel,:] = SVector(i1,i4,i3)   
                    end
                end
            end  
        end
        mesh!(   axe,rvertex, rface    , color = solid_color) 
        if draw_marking
            mark[:,nseg+2,:] .= NaN 
            lines!(  axe,rmark,color = line_color)    
        end
    end
end
