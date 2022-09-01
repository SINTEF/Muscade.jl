
using  Muscade,Muscade.Tools.Dialect,Muscade.Tools.Dots
#using  EspyInsideFunctions
using  StaticArrays,Printf,LinearAlgebra,StaticUnivariatePolynomials,GLMakie

export Ballast
#export Turbine,AnchorLine,Ballast,HawserTop,Hawser


### Ballast
struct Ballast{Tsea} <: AbstractElement
    xₘ       :: SVector{3,𝕣} # dx1,dx2,dx3
    seadrag  :: 𝕣
    sea      :: Tsea  # function
    buoyancy :: 𝕣
end
Ballast(nod::Vector{Node};seadrag,sea,buoyancy) = Ballast(coords(nod)[1,:],seadrag,sea,buoyancy)  
@espy function Muscade.lagrangian(o::Ballast, δX,X,U,A, χo,χcv, t,ε,dbg)
    x, δx    = ∂0(X)+o.xₘ, ∂0(δX)  
    δW       = (view(δx,1:2) ∘₁ o.sea(t,view(x,1:2))) * (o.seadrag + A[1])
    δW      += δx[3  ] * -(o.buoyancy+A[2])
    χn       = nostate
    return δW,χn
end
Muscade.draw(axe,key,out, o::Ballast, δX,X,U,A, χo,χcv, t,ε,dbg) = scatter!(axe,∂0(X)+o.Xₘ,markersize=50,color=:blue)
Muscade.Xdofid(  ::Type{<:Ballast}) = (nod=[1,1,1],typ=[:dx1,:dx2,:dx3])
Muscade.Adofid(  ::Type{<:Ballast}) = (nod=[2,2  ],typ=[:Δseadrag,:Δbuoyancy])
Muscade.espyable(::Type{<:Ballast}) = (x=(3,),)

### Turbine
struct Turbine{Tsea,Tsky} <: AbstractElement
    xₘ      :: SVector{3,𝕣} # dx1,dx2,rx3
    z       :: 𝕣
    seadrag :: 𝕣
    sea     :: Tsea  # function
    skydrag :: 𝕣
    sky     :: Tsky  # function
end
Turbine(nod::Vector{Node};seadrag,sea,skydrag,sky) =
        Turbine(         [coords(nod)[1,1],coords(nod)[1,2],0.],coords(nod)[1,3],seadrag,sea,skydrag,sky)  
# Turbine(SVector(coords(nod)[1,1],coords(nod)[1,2],0.),coords(nod)[1,3],seadrag,sea,skydrag,sky)  
@espy function Muscade.lagrangian(o::Turbine, δX,X,U,A, χo,χcv, t,ε,dbg)
    x, δx    = ∂0(X)+o.xₘ, ∂0(δX)  
    δW       = view(δx,1:2) ∘₁ (o.sea(t,view(x,1:2))*(o.seadrag+A[1]) + o.sky(t,view(x,1:2))*(o.skydrag+A[2]))
    χn       = nostate
    return δW,χn
end
function Muscade.draw(axe,key,out, o::Turbine, δX,X,U,A, χo,χcv, t,ε,dbg)
    x    = ∂0(X)+o.xₘ  
    lines!(axe,SMatrix{2,3}(x[1],x[1],x[2],x[2],o.z-10,o.z+10)' ,color=:orange, linewidth=5)
end
Muscade.Xdofid(  ::Type{<:Turbine}) = (nod=[1,1,1],typ=[:dx1,:dx2,:rx3])
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


### HawserTop
struct HawserTop <: AbstractElement
    xₘtop   :: SVector{3,𝕣}  # x1,x2,x3
    xₘbot   :: SVector{3,𝕣}  # x1,x2,x3
    Δxₘtop  :: SVector{3,𝕣}  # as meshed, node to fairlead
    L₀      :: 𝕣
    EA      :: 𝕣
end
HawserTop(nod::Vector{Node};Δxₘtop,L₀,EA) = HawserTop(coords(nod)[1,:],coords(nod)[2,:],Δxₘtop,L₀,EA)
@espy function Muscade.lagrangian(o::HawserTop, δX,X,U,A, χo,χcv, t,ε,dbg)
    xₘtop,xₘbot,Δxₘtop,L₀,EA = o.xₘtop,o.xₘbot,o.Δxₘtop,o.L₀,o.EA
    :Xtop    = [∂0(X)[1:2]...,0] + xₘtop
    :Xbot    =  ∂0(X)[4:6]       + xₘbot
    α        =  ∂0(X)[3]           # azimut from COG to fairlead
    δx       =  ∂0(δX)  
    c,s      = cos(α),sin(α)
    :ΔXtop   = [c -s 0;s c 0;0 0 1]*Δxₘtop   # arm of the fairlead
    :ΔX      = Xtop+ΔXtop-Xbot               # vector from fairlead to anchor
    :L       = norm(ΔX)
    dir      = ΔX./L
    :T       = EA*(L/L₀-1)
    δW       = view(δx,1:2)∘₁(dir[1:2].*T) + δx[3]*(ΔXtop[1]*r[2]-ΔXtop[2]*r[1]) + view(δx,4:6)∘₁(-dir.*T)
    χn       = nostate
    return δW,χn
end
function Muscade.draw(axe,key,out, o::HawserTop, δX,X,U,A, χo,χcv, t,ε,dbg)
    Muscade.lagrangian(out,key,o, δX,X,U,A, χo,χcv, t,ε,(dbg...,espy2draw=true))
    Xtop,Xbot,ΔXtop = out[key.Xtop], out[key.Xbot], out[key.ΔXtop]
    lines!(axe,hcat(Xtop,Xtop+ΔXtop) ,color=:red,   linewidth=2)
    lines!(axe,hcat(Xbot,Xtop+ΔXtop) ,color=:black, linewidth=2)
end
Muscade.Xdofid(      ::Type{<:HawserTop})  = (nod=[1,1,1,2,2,2],typ=[:dx1,:dx2,:rx3,:dx1,:dx2,:dx3])
Muscade.Adofid(      ::Type{<:HawserTop})  = (nod=[3,3],typ=[:ΔL,:ΔEA])
Muscade.espyable(    ::Type{<:HawserTop})  = (Xtop=(3,),Xbot=(3,),ΔXtop=(3,),ΔX=(3,),L=scalar,T=scalar)
Muscade.request2draw(::Type{<:HawserTop})  = @request (Xtop,Xbot,ΔXtop)

### Hawser
struct Hawser <: AbstractElement
    xₘ1     :: SVector{3,𝕣}  # x1,x2,x3
    xₘ2     :: SVector{3,𝕣}  # x1,x2,x3
    L₀      :: 𝕣
    EA      :: 𝕣
end
Hawser(nod::Vector{Node};L₀,EA) = Hawser(coords(nod)[1,:],coords(nod)[2,:],L₀,EA)

@espy function Muscade.lagrangian(o::Hawser, δX,X,U,A, χo,χcv, t,ε,dbg)
    xₘ1,xₘ2,L₀,EA = o.xₘ1,o.xₘ2,o.L₀,o.EA
    :X1      = ∂0(X)[1:3] + xₘ1
    :X2      = ∂0(X)[4:6] + xₘ2
    δx       = ∂0(δX)
    :ΔX      = X2-X1               # vector from fairlead to anchor
    :L       = norm(ΔX)
    dir      = ΔX./L
    :T       = EA*(L/L₀-1)
    δW       = view(δx,4:6) ∘₁ (dir.*T) - view(δx,1:3) ∘₁ (-dir.*T)
    χn       = nostate
    return δW,χn
end
function draw(axe,key,out, o::Hawser, δX,X,U,A, χo,χcv, t,ε,dbg)
    X1      = ∂0(X)[1:3] + ol.xₘ1
    X2      = ∂0(X)[4:6] + o.xₘ2
    lines!(axe,hcat(X1,X2) ,color = :black, linewidth = 2)
end
Muscade.Xdofid(      ::Type{<:Hawser})  = (nod=[1,1,1,2,2,2],typ=[:dx1,:dx2,:dx3,:dx1,:dx2,:dx3])
Muscade.Adofid(      ::Type{<:Hawser})  = (nod=[3,3],typ=[:ΔL,:ΔEA])
Muscade.espyable(    ::Type{<:Hawser})  = (X1=(3,),X2=(3,),ΔX=(3,),L=scalar,T=scalar)
      



