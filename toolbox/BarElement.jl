# # Bar element

using StaticArrays, LinearAlgebra, Muscade

# Data structure containing the cross section material properties
struct AxisymmetricBarCrossSection
    EA  :: 𝕣 # Axial stiffness [N]
    μ   :: 𝕣 # Mass per unit length [kg/m]
    w   :: 𝕣 # Weight per unit length [N/m]
    Ca₁ :: 𝕣 # Tangential added mass per unit length [kg/m]
    Cl₁ :: 𝕣 # Tangential linear damping coefficient per unit length [N/m/(m/s)]
    Cq₁ :: 𝕣 # Tangential quadratic damping coefficient per unit length [N/m/(m/s)^2], for example from drag
    Ca₂ :: 𝕣 # Tranvserse added mass per unit length [kg/m] for motions along second axis
    Cl₂ :: 𝕣 # Transverse linear damping coefficient per unit length [N/m/(m/s)] for motions along second axis
    Cq₂ :: 𝕣 # Transverse quadratic damping coefficient per unit length [N/m/(m/s)^2], for motions along second axis
    # TODO: add gravity field to bar properties (time dependent), and use it to compute the weight. This to enable static analyses. 
end
AxisymmetricBarCrossSection(;EA,μ=μ,w=0.,Ca₁=0.,Cl₁=0.,Cq₁=0.,Ca₂=0.,Cl₂=0.,Cq₂=0.) = AxisymmetricBarCrossSection(EA,μ,w,Ca₁,Cl₁,Cq₁,Ca₂,Cl₂,Cq₂);




const ngp        = 4
const ndim       = 3
const nXdof      = 6
const nUdof      = 3
const nXnod      = 2;
ψ₁(ζ) = -ζ + 1/2          
ψ₂(ζ) =  ζ + 1/2         

# Data structure describing an Bar3D element as meshed
"""
    Bar3D

An Euler bar element
"""
struct Bar3D{Mat,Uforce} <: AbstractElement
    cₘ       :: SVector{3,𝕣}     # Position of the middle of the element, as meshed
    tgₘ      :: SVector{3,𝕣}     # Vector connecting the nodes of the element in the global coordinate system
    L₀       :: 𝕣                # As-meshed length of the element
    mat      :: Mat              # Used to store material properties (AxisymmetricBarCrossSection, for example)
    wgp      :: SVector{ngp,𝕣}   # weight associated to each Gauss point
    ζgp      :: SVector{ngp,𝕣}   # Location of the Gauss points for the normalized element with length 1
    ζnod     :: SVector{nXnod,𝕣} # Location of the nodes for the normalized element with length 1
    ψ₁       :: SVector{ngp,𝕣}   # Value at gp of shape function for differential axial displacement
    ψ₂       :: SVector{ngp,𝕣}   # Value at gp of shape function for differential axial displacement
end;

# For performance, `residual` will only accept differentiation to first order
Muscade.no_second_order(::Type{<:Bar3D}) = Val(true)
 

# Define nodes, classes, and field names of dofs
Muscade.doflist(     ::Type{Bar3D{Mat,false}}) where{Mat} = 
        (inod = (1,1,1, 2,2,2), 
         class= (:X,:X,:X,:X,:X,:X), 
         field= (:t1,:t2,:t3, :t1,:t2,:t3) )
Muscade.doflist(     ::Type{Bar3D{Mat,true}}) where{Mat} = 
        (inod = (1,1,1, 2,2,2, 3,3,3), 
         class= (:X,:X,:X,:X,:X,:X,:U,:U,:U),  
         field= (:t1,:t2,:t3, :t1,:t2,:t3, :t1,:t2,:t3) )

# ElementType for the Bar3D element. Arguments: node list, material, and direction of the first bending axis in the global coordinate system.  
Bar3D(nod;kwargs...) = Bar3D{false}(nod;kwargs...) # by default, Bar3D does not have Udof.
function Bar3D{Udof}(nod::Vector{Node};mat) where {Udof}
    c       = coord(nod)
    ## Position of the middle of the element in the global coordinate system (as-meshed)
    cₘ      = SVector{3}((c[1]+c[2])/2)
    ## Tangential vector to the element in the global coordinate system, and its length (as-meshed)
    tgₘ     = SVector{3}( c[2]-c[1]   )
    L₀      = norm(tgₘ)
    ## Location ζgp of the Gauss points for a unit-length beam element, with nodes at ζnod=±1/2, and weigths. 
    wgp    = SVector{ngp}(      L₀/2*(18-sqrt(30))/36,          L₀/2*(18+sqrt(30))/36  ,        L₀/2*(18+sqrt(30))/36,          L₀/2*(18-sqrt(30))/36       ) 
    ζgp     = SVector{ngp  }(   -1/2*sqrt(3/7+2/7*sqrt(6/5)),   -1/2*sqrt(3/7-2/7*sqrt(6/5)),   +1/2*sqrt(3/7-2/7*sqrt(6/5)),   +1/2*sqrt(3/7+2/7*sqrt(6/5))) 
    ζnod    = SVector{nXnod}(   -1/2  ,1/2  )
    shapes  = (ψ₁.(ζgp), ψ₂.(ζgp))
    return Bar3D{typeof(mat),Udof}(cₘ,tgₘ,L₀,mat,wgp,ζgp,ζnod,shapes...)
end;

@espy function resultants(o::AxisymmetricBarCrossSection,ε,x,u) 
    # Inertia 
    a      = ∂2(x)
    fi      = o.μ * a
    
    δ      = ∂0(u)
    # Added mass
    a₁ = a ∘₁ δ
    a₂ = a - a₁ * δ
    fa  = SVector{3}(o.Ca₁ * a₁,o.Ca₂ * a₂[2],o.Ca₂* a₂[3])
    # Sum of external forces
    ☼fe      =   fi+fa
    # Internal forces
    ☼fᵢ      = o.EA*∂0(ε)
    return fᵢ,fe
end;

vec3(v,ind) = SVector{3}(v[i] for i∈ind);

# Define now the residual function for the Bar3D element.
@espy function Muscade.residual(o::Bar3D{Mat,Udof},   X,U,A,t,SP,dbg) where{Mat,Udof}
    P,ND    = constants(X),length(X)
    x_      = motion{P}(X) 
    uᵧ₁,uᵧ₂   = vec3(x_,1:3), vec3(x_,4:6) 
    c        = o.cₘ + 0.5*(uᵧ₁+uᵧ₂) 
    tg      = o.tgₘ + uᵧ₂ - uᵧ₁
    L       = √((o.L₀+uᵧ₂[1]-uᵧ₁[1])^2+(uᵧ₂[2]-uᵧ₁[2])^2+(uᵧ₂[3]-uᵧ₁[3])^2)
    ε_       = L/o.L₀ - 1
    u_       = tg/L
    ε,u = motion⁻¹{P,ND}(ε_,u_)
    uval = ∂0(u)
    ε∂X₀ = 1/o.L₀*SVector{6}(-uval[1],-uval[2],-uval[3],uval[1],uval[2],uval[3]) # how strains vary with nodal displacements (will be used in the Princple of Virtual Work, PVW)

    # Compute Gauss point kinematics
    gp = ntuple(ngp) do igp; 
        x = c + tg * o.ζgp[igp]; 
        @named(x); 
    end
    # Compute loads at Gauss points
    gpContrib = ntuple(ngp) do igp
        ζ = o.ζgp[igp] # Coordinate of the Gauss point along [-1/2,1/2]
        x∂X₀ = SMatrix{3,6}(ψ₁(ζ),0,0, 0,ψ₁(ζ),0, 0,0,ψ₁(ζ), ψ₂(ζ),0,0, 0,ψ₂(ζ),0, 0,0,ψ₂(ζ))   # how motions of Gauss point vary with nodal displacements
        
        x = motion⁻¹{P,ND}(gp[igp].x)     # Physical location of the Gauss point 
        fᵢ,fₑ     = ☼resultants(o.mat,ε,x,u)  # compute loads from strains/motions, etc.
        
        fₑ        = Udof ? fₑ-∂0(U) : fₑ                       # U is per unit length
        R_        = ( fᵢ ∘₀ ε∂X₀ + fₑ ∘₁ x∂X₀ ) * o.wgp[igp]   #  Application of PVW, local contribution of the integral over the element
        @named(R_);
    end
    R                   = sum(gpᵢ.R_ for gpᵢ∈gpContrib) 
    return R,noFB  
end;


# The following functions explain how the bar element should be drawn
# using GLMakie
# """

# Drawing a `Bar3D`.

#     draw!(axis,state)

#     draw!(axis,state;Bar3D=(;style=:shape))

#     α      = 2π*(0:19)/20
#     circle = 0.1*[cos.(α) sin.(α)]'
#     draw!(axis,state;Bar3D=(;style=:solid,section = circle))

# `style=:shape` shows the deformed neutral axis of the element. It has optional arguments `frame=true` 
# (draws the element's corotated frame of reference)
# and `nseg=10` (number of points to show the deflected shape of each element). 

# `style=:solid` shows the deformed shape of the element. It requires the input `section=...` to be given
# a matrix of size `(2,nsec)` describing `nsec` points around the cross section of the element (no need to close 
# the circumference by repeating the first point at the end).  It has optional arguments `nseg=10` as above, `marking=true`
# to draw a longitudinal marking and `solid_color=:yellow`.
 
# Other optional arguments (and their default values) are
# - `Udof` (`true` iff element has Udofs) wether to draw U-forces.
# - `draw_frame = false` wether to draw the local reference frame of each element
# - `draw_marking = true` wether to draw "longitudinal marking" along the element.  Will only draw if style=:solid.
# - `nseg = 1` number of segments to display the shape of a deformed element
# - `solid_color = :yellow` color of the surface if `style=:solid`
# - `line_color = :black` color of the line if `style=:sshape`
# - `Uscale = 1.` How many meter is a Newton per meter?
# """
# function Muscade.allocate_drawing(axis,o::AbstractVector{Bar3D{Tmat,Udof}};kwargs...) where{Tmat,Udof}
#     args                 = default{:Bar3D     }(kwargs,(;)     )  
#     section              = default{:section         }(args,zeros(2,0))  
#     nsec                 = size(section,2)                            
#     opt = (default(args,(style=:shape,draw_frame=false,draw_marking=true,nseg=1,
#                   solid_color=:yellow,line_color=:black,Uscale=1.,Udof=Udof))...,
#             nel          = length(o)                                  ,
#             nsec         = nsec                                       ,                    
#             section      = section                                    ,
#             markrad      = nsec==0 ? 0. : 1.01*maximum(section[1,:])      
#         )
#     opt.style==:solid && nsec<2 && muscadeerror("An section description must be provided for 'solid' plot")
#     nel_shape         = opt.style==:shape ? opt.nel   : 0
#     nel_shape_frame   = opt.draw_frame    ? nel_shape : 0
#     nel_solid         = opt.style==:solid ? opt.nel   : 0 
#     nel_solid_marking = opt.draw_marking  ? nel_solid : 0
#     nel_udof          = opt.Udof          ? opt.nel   : 0

#     mut=(
#             node         = 𝕣2(undef,3,3*opt.nel)                        ,
#             shape_x      = 𝕣2(undef,3,(opt.nseg+2)*nel_shape)           ,   
#             shape_frame  = 𝕣2(undef,3,3*3*nel_shape_frame)              , # idim, point-point-lift, ivec, iel
#             solid_vertex = 𝕣2(undef,3,opt.nsec*(opt.nseg+1)*nel_solid)  , 
#             solid_face   = 𝕫2(undef,2*opt.nsec* opt.nseg   *nel_solid,3),
#             solid_mark   = 𝕣2(undef,3,(opt.nseg+2)*nel_solid_marking)   ,     
#             ucrest       = 𝕣2(undef,3,5*nel_udof)                       , # idim, 6point-lift,iel
#         )   
#     return mut,opt
# end

# function Muscade.update_drawing(axis,o::AbstractVector{Bar3D{Tmat,Udof}},oldmut,opt, Λ,X,U,A,t,SP,dbg) where{Tmat,Udof} 
#     mut               = oldmut 
#     X₀                = ∂0(X)
#     U₀                = ∂0(U)
#     it1,ir1,it2,ir2   = SVector{3}(1:3),SVector{3}(4:6),SVector{3}(7:9),SVector{3}(10:12)
#     nsec              = size(opt.section,2) 
#     node = reshape(mut.node,(3,3,opt.nel))
#     for (iel,oᵢ) = enumerate(o)
#         node[:,1,iel] = oᵢ.cₘ - oᵢ.tgₘ/2 + X₀[it1,iel]
#         node[:,2,iel] = oᵢ.cₘ + oᵢ.tgₘ/2 + X₀[it2,iel]
#         node[:,3,iel].= NaN  
#     end

#     if opt.style==:shape
#         ζ = range(-1/2,1/2,opt.nseg+1)
#         if opt.draw_frame shape_frame  = reshape(mut.shape_frame ,(3,3,3       ,opt.nel)) end
#         if opt.Udof       ucrest       = reshape(mut.ucrest,      (3,5         ,opt.nel)) end
#         shape_x                        = reshape(mut.shape_x     ,(3,opt.nseg+2,opt.nel))
#         for (iel,oᵢ) = enumerate(o)
#             cₘ,rₘ,tgₘ,tgₑ,ζnod,ζgp,L  = oᵢ.cₘ,oᵢ.rₘ,oᵢ.tgₘ,oᵢ.tgₑ,oᵢ.ζnod,oᵢ.ζgp,oᵢ.L   
#             X₀ₑ = view(X₀,:,iel)
#             vₛₘ,rₛₘ,uₗ₂,vₗ₂,cₛₘ = corotated{:direct}(oᵢ,X₀ₑ) 
#             if opt.draw_frame
#                 for ivec = 1:3
#                     shape_frame[:,1,ivec,iel] = cₛₘ
#                     shape_frame[:,2,ivec,iel] = cₛₘ + oᵢ.L/3*rₛₘ[:,ivec]
#                     shape_frame[:,3,ivec,iel].= NaN
#                 end
#             end
#             if opt.Udof
#                 ucrest[:,1,iel] = node[:,1,iel]
#                 ucrest[:,2,iel] = node[:,1,iel] + rₛₘ ∘₁ view(U₀,:,iel) * opt.Uscale
#                 ucrest[:,3,iel] = node[:,2,iel] + rₛₘ ∘₁ view(U₀,:,iel) * opt.Uscale
#                 ucrest[:,4,iel] = node[:,2,iel]
#                 ucrest[:,5,iel].= NaN
#             end
#             for (i,ζᵢ) ∈ enumerate(ζ)
#                 y          = SVector(yₐ(ζᵢ)*uₗ₂[1] , yᵤ(ζᵢ)*uₗ₂[2]+L*yᵥ(ζᵢ)*vₗ₂[3], yᵤ(ζᵢ)*uₗ₂[3]-L*yᵥ(ζᵢ)*vₗ₂[2])  
#                 shape_x[:,i         ,iel] = rₛₘ∘₁(tgₑ*ζᵢ+y)+cₛₘ 
#                 shape_x[:,opt.nseg+2,iel].= NaN
#             end        
#         end
#     elseif opt.style==:solid
#         ζ = range(-1/2,1/2,opt.nseg+1)
#         idx(iel,iseg,isec) = mod_onebased(isec,opt.nsec)+opt.nsec*(iseg-1+(opt.nseg+1)*(iel-1)) # 1st index into rvertex
#         if opt.Udof         ucrest         = reshape(mut.ucrest       ,(3,5          ,opt.nel)) end
#         if opt.draw_marking solid_mark     = reshape(mut.solid_mark  ,(3,opt.nseg+2 ,opt.nel)) end
#         solid_face                         = reshape(mut.solid_face  ,(2,opt.nsec, opt.nseg   ,opt.nel,3))
#         solid_vertex                       = reshape(mut.solid_vertex,(3,opt.nsec, opt.nseg+1 ,opt.nel))
#         for (iel,oᵢ) = enumerate(o)
#             cₘ,rₘ,tgₘ,tgₑ,ζnod,ζgp,L  = oᵢ.cₘ,oᵢ.rₘ,oᵢ.tgₘ,oᵢ.tgₑ,oᵢ.ζnod,oᵢ.ζgp,oᵢ.L   
#             X₀ₑ = view(X₀,:,iel)
#             vₛₘ,rₛₘ,uₗ₂,vₗ₂,cₛₘ = corotated{:direct}(oᵢ,X₀ₑ) 
#             vᵧ₁,vᵧ₂          = vec3(X₀ₑ,4:6), vec3(X₀ₑ,10:12)
#             rₛ₁              = Rodrigues(vᵧ₁)
#             rₛ₂              = Rodrigues(vᵧ₂)
#             if opt.Udof
#                 ucrest[:,1,iel] = node[:,1,iel]
#                 ucrest[:,2,iel] = node[:,1,iel] + rₛₘ ∘₁ view(U₀,:,iel) * opt.Uscale
#                 ucrest[:,3,iel] = node[:,2,iel] + rₛₘ ∘₁ view(U₀,:,iel) * opt.Uscale
#                 ucrest[:,4,iel] = node[:,2,iel]
#                 ucrest[:,5,iel].= NaN
#             end
#             Δv = Rodrigues⁻¹(rₛ₂ ∘₁ rₛ₁')/opt.nseg
#             for (iseg,ζᵢ) ∈ enumerate(ζ) # actualy iterating over nseg+1 segment boundaries
#                 y  = SVector(yₐ(ζᵢ)*uₗ₂[1] , yᵤ(ζᵢ)*uₗ₂[2]+L*yᵥ(ζᵢ)*vₗ₂[3], yᵤ(ζᵢ)*uₗ₂[3]-L*yᵥ(ζᵢ)*vₗ₂[2])  # interpolate
#                 xn = rₛₘ∘₁(tgₑ*ζᵢ+y)+cₛₘ # point on neutral axis
#                 r  = Rodrigues((iseg-1)*Δv) ∘₁ rₛ₁ ∘₁ rₘ  
#                 if opt.draw_marking 
#                     solid_mark[:,    iseg  ,iel] = xn .+ r[:,2]*opt.markrad 
#                     solid_mark[:,opt.nseg+2,iel].= NaN 
#                 end
#                 for isec = 1:opt.nsec
#                     solid_vertex[:,isec,iseg,iel] = xn .+ r[:,2]*opt.section[1,isec] + r[:,3]*opt.section[2,isec] 
#                     if iseg≤opt.nseg
#                         i1,i2,i3,i4 = idx(iel,iseg,isec),idx(iel,iseg  ,isec+1),idx(iel,iseg+1,isec  ),idx(iel,iseg+1,isec+1)
#                         solid_face[1,isec,iseg,iel,:] = SVector(i1,i2,i4)    
#                         solid_face[2,isec,iseg,iel,:] = SVector(i1,i4,i3)   
#                     end
#                 end
#             end  
#         end
#     end
#     return mut
# end

# function Muscade.display_drawing!(axis,::Type{Bar3D{Tmat,Udof}},obs,opt) where{Tmat,Udof}
#     scatter!(                                          axis, obs.node                         ,color = opt.line_color , marker=:circle,markersize=3)  
#     opt.style==:shape  &&                     lines!(  axis, obs.shape_x                      ,color = opt.line_color ,linewidth=.5                )
#     opt.style==:shape  && opt.draw_frame   && lines!(  axis, obs.shape_frame                  ,color = :grey          ,linewidth=.5                )    
#     opt.style==:solid  &&                     mesh!(   axis, obs.solid_vertex, obs.solid_face ,color = opt.solid_color                             )  
#     opt.style==:solid  && opt.draw_marking && lines!(  axis, obs.solid_mark                   ,color = opt.line_color                              )    
#     opt.Udof           &&                     lines!(  axis, obs.ucrest                       ,color = :red           ,linewidth=.5                )    
# end



