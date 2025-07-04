#Assumes that BeamElements.jl has been included previously, and that "using GLMakie" has been invoked
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
- `Udof = true` wether to draw U-forces
- `draw_frame = false` wether to draw the local reference frame of each element
- `draw_marking = true` wether to draw "longitudinal marking" along the element.  Will only draw if style=:solid.
- `nseg = 1` number of segments to display the shape of a deformed element
- `solid_color = :yellow` color of the surface if `style=:solid`
- `line_color = :black` color of the line if `style=:sshape`
- `Uscale = 1.` How many meter is a Newton per meter?
"""
function Muscade.allocate_drawing(axis,o::AbstractVector{EulerBeam3D{Tmat,Udof}};kwargs...) where{Tmat,Udof}
    # define constant inputs to the drawing process
    args                 = default{:EulerBeam3D     }(kwargs,(;)     )  
    section              = default{:section         }(args,zeros(2,0))  
    nsec                 = size(section,2)                            
    opt=(
            nel          = length(o)                                  ,
            style        = default{:style           }(args,:shape    ),                       
            draw_frame   = default{:frame           }(args,false     ),                       
            draw_marking = default{:marking         }(args,true      ),                       
            nseg         = default{:nseg            }(args,1         ),                       
            solid_color  = default{:solid_          }(args,:yellow   ),                       
            line_color   = default{:line_color      }(args,:black    ),                       
            Uscale       = default{:Uscale          }(args,1.        ),   
            Udof         = Udof && default{:Udof    }(args,true      ),                                    
            nsec         = nsec                                       ,                    
            section      = section                                    ,
            markrad      = nsec==0 ? 0. : 1.01*maximum(section[1,:])      
        )
    opt.style==:solid && nsec<2 && muscadeerror("An section description must be provided for 'solid' plot")

    # we are going to allocate many arrays. The plotting process has options about what to draw and what to leave out.
    # To save memory, we set nel_something to zero if the corresponding arrays are not needed.
    nel_shape         = opt.style==:shape ? opt.nel   : 0
    nel_shape_frame   = opt.draw_frame    ? nel_shape : 0
    nel_solid         = opt.style==:solid ? opt.nel   : 0 
    nel_solid_marking = opt.draw_marking  ? nel_solid : 0
    nel_udof          = opt.Udof          ? opt.nel   : 0

    # This tuple contains all the "mutables", whose contents will change from step to step
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
    mut           = oldmut 
    X₀            = ∂0(X)
    U₀            = ∂0(U)
    it1,ir1,it2,ir2 = SVector{3}(1:3),SVector{3}(4:6),SVector{3}(7:9),SVector{3}(10:12)
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
            vₛₘ,rₛₘ,uₗ₂,vₗ₂,cₛₘ = corotated(oᵢ,X₀ₑ) 
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
            vₛₘ,rₛₘ,uₗ₂,vₗ₂,cₛₘ = corotated(oᵢ,X₀ₑ) 
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
            for (iseg,ζᵢ) ∈ enumerate(ζ)
                y  = SVector(yₐ(ζᵢ)*uₗ₂[1] , yᵤ(ζᵢ)*uₗ₂[2]+L*yᵥ(ζᵢ)*vₗ₂[3], yᵤ(ζᵢ)*uₗ₂[3]-L*yᵥ(ζᵢ)*vₗ₂[2])  
                xn = rₛₘ∘₁(tgₑ*ζᵢ+y)+cₛₘ # point on neutral axis
                v  = (iseg-1)/opt.nseg*Rodrigues⁻¹(rₛ₂ ∘₁ rₛ₁')
                r  = Rodrigues(v) ∘₁ rₛ₁ ∘₁ rₘ  
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


