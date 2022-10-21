
using  Muscade,Muscade.Tools.Dialect,Muscade.Tools.Dots
#using  EspyInsideFunctions
using  StaticArrays,Printf,LinearAlgebra,StaticUnivariatePolynomials,GLMakie

export Turbine,AnchorLine

### Turbine

struct Turbine{Tsea,Tsky} <: AbstractElement
    xₘ      :: SVector{2,𝕣} # dx1,dx2
    z       :: 𝕣
    seadrag :: 𝕣
    sea     :: Tsea  # function
    skydrag :: 𝕣
    sky     :: Tsky  # function
end
Turbine(nod::Vector{Node};seadrag,sea,skydrag,sky) = Turbine([coords(nod)[1,1],coords(nod)[1,2],0.],coords(nod)[1,3],seadrag,sea,skydrag,sky)  
@espy function Muscade.lagrangian(o::Turbine, δX,X,U,A, χo,χcv, t,ε,dbg)
    :x, δx   = ∂0(X)+o.xₘ, ∂0(δX)  
    δW       = δx ∘₁ (o.sea(t,x)*(o.seadrag+A[1]) + o.sky(t,x)*(o.skydrag+A[2]))
    χn       = nostate
    return δW,χn
end
function Muscade.draw(axe,key,out, o::Turbine, δX,X,U,A, χo,χcv, t,ε,dbg)
    x    = ∂0(X)+o.xₘ  
    lines!(axe,SMatrix{2,3}(x[1],x[1],x[2],x[2],o.z-10,o.z+10)' ,color=:orange, linewidth=5)
end
Muscade.Xdofid(  ::Type{<:Turbine}) = (nod=[1,1,1],typ=[:dx1,:dx2])
Muscade.Adofid(  ::Type{<:Turbine}) = (nod=[2,2  ],typ=[:Δseadrag,:Δskydrag])
Muscade.espyable(::Type{<:Turbine}) = (x=(3,),)

### AnchorLine

struct AnchorLine <: AbstractElement
    xₘtop   :: SVector{3,𝕣}  # x1,x2,x3
    Δxₘtop  :: SVector{3,𝕣}  # as meshed, node to fairlead
    xₘbot   :: SVector{2,𝕣}  # x1,x2 (x3=0)
    L       :: 𝕣
    buoyancy:: 𝕣
end
AnchorLine(nod::Vector{Node};Δxₘtop,xₘbot,L,buoyancy) = AnchorLine(coords(nod)[1,:],Δxₘtop,xₘbot,L,buoyancy)
#AnchorLine(coords(nod)[1,:],SVector{3}(Δxₘtop),SVector{2}(xₘbot),L,buoyancy)
         
p = Polynomial(   2.82040487827,  -24.86027164695,   153.69500343165, -729.52107422849, 2458.11921356871,
              -5856.85610233072, 9769.49700812681,-11141.12651712473, 8260.66447746395,-3582.36704093187,
                687.83550335374)

@espy function Muscade.lagrangian(o::AnchorLine, δX,X,U,A, χo,χcv, t,ε,dbg)
    xₘtop,Δxₘtop,xₘbot,L,buoyancy = o.xₘtop,o.Δxₘtop,o.xₘbot,o.L,o.buoyancy      # a for anchor, t for TDP, f for fairlead
    x, δx    = ∂0(X), ∂0(δX)  
    :Xtop    = [x[1:2]...,0] + xₘtop
    α        =  x[3]                            # azimut from COG to fairlead
    c,s      = cos(α),sin(α)
    :ΔXtop   = [c -s 0;s c 0;0 0 1]*Δxₘtop       # arm of the fairlead
    :ΔXchain = Xtop[1:2]+ΔXtop[1:2]-xₘbot        # vector from anchor to fairlead
    :xaf     = norm(ΔXchain)                     # horizontal distance from anchor to fairlead
    :cr      = exp10(p((L-xaf)/Xtop[3]))*Xtop[3] # curvature radius at TDP
    :Fh      = -cr*buoyancy                      # horizontal force
    :ltf     = √(Xtop[3]^2+2Xtop[3]*cr)          # horizontal distance from fairlead to TDP
    δW       = view(δx,1:2) ∘₁ (ΔXchain/xaf.*Fh) + δx[3] * (ΔXtop[1]*r[2]-ΔXtop[2]*r[1])
    χn       = nostate
    return δW,χn
end
function Muscade.draw(axe,key,out, o::AnchorLine, δX,X,U,A, χo,χcv, t,ε,dbg)
    Muscade.lagrangian(out,key,o, δX,X,U,A, χo,χcv, t,ε,(dbg...,espy2draw=true))
    Laf,Xbot,Xtop,ΔXtop,ΔXchain,cr,xaf,Ltf = o.L, o.xₘbot, out[key.Xtop],out[key.ΔXtop],out[key.ΔXchain], out[key.cr], out[key.xaf], out[key.ltf]
    n     = ΔXchain./xaf  # horizontal normal vector from anchor to fairlead
    xat   = Laf-Ltf
    xtf   = xaf-xat
    Xtdp  = Xbot + n*xat
    x    = range(max(0,-xat),xtf,11)
    X    = n[1].*x.+Xtdp[1]
    Y    = n[2].*x.+Xtdp[2]
    Z    = cr.*(cosh.(x./cr).-1)
    lines!(axe,hcat(X,Y,Z)     ,color=:blue,  linewidth=2) # line
    scatter!(axe,Xtdp          ,markersize=20,color=:blue)
    scatter!(axe,Xbot          ,markersize=50,color=:red)
    if xat>0
        lines!(axe,hcat(Xbot,Xtdp)       ,color=:green, linewidth=2) # seafloor
    else
        x    = range(0,-xat,11)
        X    = n[1].*x.+Xtdp[1]
        Y    = n[2].*x.+Xtdp[2]
        Z    = cr.*(cosh.(x./cr).-1)
        lines!(axe,hcat(X,Y,Z)           ,color=:red,  linewidth=5) # line
    end
    lines!(axe,hcat(Xtop,Xtop+ΔXtop) ,color=:red , linewidth=2) # excentricity
end
Muscade.Xdofid(      ::Type{<:AnchorLine}) = (nod=[1,1,1],typ=[:dx1,:dx2,:rx3])
Muscade.Adofid(      ::Type{<:AnchorLine}) = (nod=[2,2  ],typ=[:ΔL,:Δbuoyancy])
Muscade.espyable(    ::Type{<:AnchorLine}) = (Xtop=(3,),ΔXtop=(3,),ΔXchain=(2,),xaf=scalar,cr=scalar,Fh=scalar,ltf=scalar)
Muscade.request2draw(::Type{<:AnchorLine}) = @request (Xtop,ΔXtop,ΔXchain,cr,xaf,ltf)



