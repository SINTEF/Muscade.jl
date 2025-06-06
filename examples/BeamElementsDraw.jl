#Assumes that BeamElements.jl has been included previously
using GLMakie

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
function Muscade.draw(axe,o::Vector{EulerBeam3D{Tmat,Udof}}, Λ,X,U,A,t,SP,dbg;kwargs...) where{Tmat,Udof}
    nel           = length(o)
    args          = default{:EulerBeam3D}(kwargs,(;)     )
    style         = default{:style      }(args,:simple   )
    draw_frame    = default{:frame      }(args,true      )
    draw_marking  = default{:marking    }(args,true      )
    nseg          = default{:nseg       }(args,10        )
    section       = default{:section    }(args,zeros(2,0))
    solid_color   = default{:color      }(args,:yellow   )
    line_color    = default{:color      }(args,:black    )
    Uscale        = default{:Uscale     }(args,1.        )
    nsec          = size(section,2)
    X₀            = ∂0(X)
    U₀            = ∂0(U)
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
        rx = reshape(x,(3,(nseg+2)*nel))
        if draw_frame  
            frame = Array{𝕣,4}(undef,3,3,3,nel)  # idim, point-point-lift, ivec, iel
            rframe = reshape(frame,(3,9nel)) 
        end  
        if Udof
            ucrest  = Array{𝕣,3}(undef,3,6,nel) # idim, 6point-lift,iel
            rucrest = reshape(ucrest,(3,6nel)) 
        end  
        for (iel,oᵢ) = enumerate(o)
            cₘ,rₘ,tgₘ,tgₑ,ζnod,ζgp,L  = oᵢ.cₘ,oᵢ.rₘ,oᵢ.tgₘ,oᵢ.tgₑ,oᵢ.ζnod,oᵢ.ζgp,oᵢ.L   
            X₀ₑ = view(X₀,:,iel)
            vₛₘ,rₛₘ,uₗ₂,vₗ₂,cₛₘ = corotated(oᵢ,X₀ₑ) 
            if draw_frame
                for ivec = 1:3
                    frame[:,1,ivec,iel] = cₛₘ
                    frame[:,2,ivec,iel] = cₛₘ + oᵢ.L/3*rₛₘ[:,ivec]
                    frame[:,3,ivec,iel].= NaN
                end
            end
            if Udof
                ucrest[:,1,iel] = cₛₘ - oᵢ.L/2*rₛₘ[:,1]
                ucrest[:,2,iel] = cₛₘ - oᵢ.L/2*rₛₘ[:,1] + rₛₘ ∘₁ view(U₀,:,iel) * Uscale
                ucrest[:,3,iel] = cₛₘ + oᵢ.L/2*rₛₘ[:,1] + rₛₘ ∘₁ view(U₀,:,iel) * Uscale
                ucrest[:,4,iel] = cₛₘ + oᵢ.L/2*rₛₘ[:,1]
                ucrest[:,5,iel] = cₛₘ - oᵢ.L/2*rₛₘ[:,1]
                ucrest[:,6,iel].= NaN
            end
            for (i,ζᵢ) ∈ enumerate(ζ)
                y          = SVector(yₐ(ζᵢ)*uₗ₂[1] , yᵤ(ζᵢ)*uₗ₂[2]+L*yᵥ(ζᵢ)*vₗ₂[3], yᵤ(ζᵢ)*uₗ₂[3]-L*yᵥ(ζᵢ)*vₗ₂[2])  
                x[:,i,iel] = rₛₘ∘₁(tgₑ*ζᵢ+y)+cₛₘ 
            end        
            x[:,nseg+2,iel] .= NaN
        end
        if Udof
            lines!(axe,rucrest,color = :red,linewidth=1)    
        end
        if draw_frame  
            lines!(axe,rframe,color = :grey,linewidth=1)    
        end
        lines!(  axe,rx,color = line_color)
        xnod  = view(x,:,[1,nseg+1],:)
        rxnod = reshape(xnod,(3,2*nel))
        scatter!(axe,rxnod,color = line_color, marker=:circle) 
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
        if Udof
            ucrest  = Array{𝕣,3}(undef,3,6,nel) # idim, 6point-lift,iel
            rucrest = reshape(ucrest,(3,6nel)) 
        end  
        for (iel,oᵢ) = enumerate(o)
            cₘ,rₘ,tgₘ,tgₑ,ζnod,ζgp,L  = oᵢ.cₘ,oᵢ.rₘ,oᵢ.tgₘ,oᵢ.tgₑ,oᵢ.ζnod,oᵢ.ζgp,oᵢ.L   
            X₀ₑ = view(X₀,:,iel)
            vₛₘ,rₛₘ,uₗ₂,vₗ₂,cₛₘ = corotated(oᵢ,X₀ₑ) 
            vᵧ₁,vᵧ₂          = vec3(X₀ₑ,4:6), vec3(X₀ₑ,10:12)
            rₛ₁              = Rodrigues(vᵧ₁)
            rₛ₂              = Rodrigues(vᵧ₂)
            if Udof
                ucrest[:,1,iel] = cₛₘ - oᵢ.L/2*rₛₘ[:,1]
                ucrest[:,2,iel] = cₛₘ - oᵢ.L/2*rₛₘ[:,1] + rₛₘ ∘₁ view(U₀,:,iel) * Uscale
                ucrest[:,3,iel] = cₛₘ + oᵢ.L/2*rₛₘ[:,1] + rₛₘ ∘₁ view(U₀,:,iel) * Uscale
                ucrest[:,4,iel] = cₛₘ + oᵢ.L/2*rₛₘ[:,1]
                ucrest[:,5,iel] = cₛₘ - oᵢ.L/2*rₛₘ[:,1]
                ucrest[:,6,iel].= NaN
            end
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
        if Udof
            lines!(axe,rucrest,color = :red,linewidth=1)    
        end
    end
end
