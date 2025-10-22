using StaticArrays, LinearAlgebra, Muscade, Muscade.Toolbox

struct Position3D{Nsensor,Nel} <: AbstractElement
    xₘ          :: SVector{3,𝕣}     # As-meshed position
    P           :: SMatrix{3,Nsensor,𝕣,Nel}
    D           :: SMatrix{3,Nsensor,𝕣,Nel}
end
"""
    .Toolbox.Position3D

An 3D single-node element, to be connected to a node with both translation and
rotation dofs.  The element makes zero contribution to residual or Lagrangian.  It
only provides requestables allowing to model accelerometers and optical position
measurements. For an inverse analysis, this is done by including the element in 
an `ElementCost`.

# Keyword arguments when adding elements:

- `P` `SMatrix{3,Nsensor,𝕣}`, giving the offset between the nodal position and the point(s)
    at which position(s) and/or accelerations will be measured.
- `D` `SMatrix{3,Nsensor,𝕣}`, the orientation of an accelerometer. If no accelerometer is 
    present at a given postion `P`, use NaNs.  If multiple accelerometers are present
    at the same position, one can repeat the columns of `P`.

# Requestables:

- `a` `a[isensor]` is the 
- `x`  `x[oder+1][:,isensor]` contains the position, velocity and acceleration (oder=0,1,2) of the sensor
       which position was described in the `isensor`-th column of `P`. 
- `rᵢ` `rₑ[oder+1]` is a vector containing a zero vector, the intrinsic rotation rate vector and its time derivative (oder=0,1,2), for the element's node.  
- `rₑ` `rₑ[oder+1]` is a vector containing the rotation vector, the extrinsic rotation rate vector and its time derivative (oder=0,1,2), for the element's node.  
- `R ` `R[oder+1]` is a matrix containing the rotation matrix, spin matrix and its time derivative (oder=0,1,2), for the element's node.  
- `xₙ` `x[oder+1]` is a vector containing the position, velocity and acceleration (oder=0,1,2) of the element's node.

#Keyword arguments when drawing elements:

When calling `draw!`, instruction to `Position3D` elements are given as in this example:

    draw!(axis,state>;Position3D=(L= @SVector [0.1,0,0],point_size=10,accelerometer_color=:orange))`

The opional keword arguments and their default values are    
- `L                   = 0.` (thus if not given, directions of accelerometers are not shown)
- `point_size          = 6`
- `point_color         = :black`
- `stalk_color         = :grey`
- `stalk_width         = 1`
- `accelerometer_color = :teal`
- `accelerometer_width = 2`
    
"""
function Position3D(node;P::SMatrix{3,Nsensor,𝕣},D::SMatrix{3,Nsensor,𝕣}) where{Nsensor}
    xₘ = SVector{3,𝕣}(coord(node)[1])
    return Position3D(xₘ,P,Muscade.colnormalize(D))
end
Muscade.doflist(::Type{<:Position3D}) = (inod  = (1  ,1  ,1  ,1  ,1  ,1  ), 
                                         class = (:X ,:X ,:X ,:X ,:X ,:X ), 
                                         field = (:t1,:t2,:t3,:r1,:r2,:r3) )
@espy function Muscade.residual(o::Position3D{Nsensor},   X,U,A,t,SP,dbg) where{Nsensor}
    P,ND   = constants(X),length(X)
    X_     = motion{P}(X)
    Δx_    = X_[SVector(1,2,3)]
    r_     = X_[SVector(4,5,6)]
    R_     = Rodrigues(r_)
    xₙ_    = o.xₘ + Δx_
    x_     = R_ ∘₁ o.P .+ xₙ_         # [3,Nsensor] of adiffs wrt time
    ☼x     = motion⁻¹{P,ND}(x_)       # a tuple of points position [3,Nsensor] matrix and its time derivatives 
    ♢x₀    = x[1]
    ♢xₙ    = motion⁻¹{P,ND}(xₙ_)
    ☼R     = motion⁻¹{P,ND}(R_)       # a tuple of node   rotation [3]         vec    and its time derivatives
    ♢rₑ    = motion⁻¹{P,ND}(r_)       # a tuple of node   rotation [3]         vec    and its time derivatives
    ♢rᵢ    = intrinsicrotationrates(R)
    ♢a     = o.D' ∘₁ ∂2(x)            # acceleration components in directions D
    return SVector{6,𝕣}(0 for i=1:6) ,noFB  
end

using GLMakie

function Muscade.allocate_drawing(axis,o::AbstractVector{Position3D{Nsensor,Nel}};kwargs...) where{Nsensor,Nel}
    args   = default{:Position3D}(kwargs,(;))
    opt    = default(args,(L                   = 0.,
                           point_size          = 6,
                           point_color         = :black,
                           stalk_color         = :grey,
                           stalk_width         = 1,
                           accelerometer_color = :teal,
                           accelerometer_width = 2)
                    )
    nel    = length(o)
    mut    = (pos = 𝕣2(undef,3,3*Nsensor*nel), 
              dir = 𝕣2(undef,3,3*Nsensor*nel),
              pts = 𝕣2(undef,3,  Nsensor*nel))
    return mut,opt
end

function Muscade.update_drawing(axis,o::AbstractVector{Position3D{Nsensor,Nel}},mut,opt, Λ,X,U,A,t,SP,dbg) where{Nsensor,Nel}
    X₀    = ∂0(X)
    U₀    = ∂0(U)
    nXdof,nUdof,nAdof = 6,0,0
    nel   = size(X₀,2)
    pos   = reshape(mut.pos,(3,3,Nsensor,nel))
    dir   = reshape(mut.dir,(3,3,Nsensor,nel))
    pts   = reshape(mut.pts,(3,  Nsensor,nel))
    req   = @request (x,xₙ,R)
    for (iel,oᵢ) ∈ enumerate(o) 
        X₀ᵢ = SVector{nXdof,𝕣}(X₀[idof,iel] for idof=1:nXdof)
        U₀ᵢ = SVector{nUdof,𝕣}(U₀[idof,iel] for idof=1:nUdof)
        Aᵢ  = SVector{nAdof,𝕣}(A[ idof,iel] for idof=1:nAdof)
        _,_,eleres = Muscade.residual(oᵢ,(X₀ᵢ,),(U₀ᵢ,),Aᵢ,t,SP,(caller=:update_drawing,dbg...),req)
        x   = ∂0(eleres.x)
        xₙ  = ∂0(eleres.xₙ)
        R   = ∂0(eleres.R)
        for isen = 1:Nsensor
            pos[:,1,isen,iel] .= xₙ
        end
        pos[:,2,:,iel] .= x  
        pos[:,3,:,iel] .= NaN  

        pts[:,:  ,iel] .= x

        dir[:,1,:,iel] .= x
        for isen = 1:Nsensor
            dir[:,2,isen,iel] .= x[:,isen] .+ opt.L * (R ∘₁ oᵢ.D[:,isen])    
        end
        dir[:,3,:,iel] .= NaN  
    end    
    return mut
end

function Muscade.display_drawing!(axis,::Type{<:Position3D},obs,opt) 
    scatter!(axis,obs.pts,markersize=opt.point_size,color=opt.point_color)
    lines!(  axis,obs.pos,color=opt.stalk_color,linewidth=opt.stalk_width)
    lines!(  axis,obs.dir,color=opt.accelerometer_color,linewidth=opt.accelerometer_width)
end

