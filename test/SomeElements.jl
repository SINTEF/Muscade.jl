using  Muscade
using  StaticArrays,LinearAlgebra #,GLMakie 

function horner(p::AbstractVector{Rp},x::Rx) where{Rp<:ℝ,Rx<:ℝ} # avoiding to use e.g. Polynomials.jl just for test code
    R = promote_type(Rp,Rx)
    y = zero(R) 
    for i ∈ reverse(p) 
        y = i + x*y
    end
    return y
end    

### Turbine

struct Turbine{Tsea,Tsky} <: AbstractElement
    xₘ      :: SVector{2,𝕣} # tx1,tx2
    z       :: 𝕣
    seadrag :: 𝕣
    sea     :: Tsea  # function
    skydrag :: 𝕣
    sky     :: Tsky  # function
end
Turbine(nod::Vector{Node};seadrag,sea,skydrag,sky) = Turbine(SVector(coord(nod)[1][1],coord(nod)[1][2]),coord(nod)[1][3],seadrag,sea,skydrag,sky)  
@espy function Muscade.residual(o::Turbine, X,U,A, t,SP,dbg)
    ☼x = ∂0(X)+o.xₘ  
    R  = -o.sea(t,x)*o.seadrag*(1+A[1]) - o.sky(t,x)*o.skydrag*(1+A[2])
    return R,noFB 
end
function Muscade.allocate_drawing(axe,o::AbstractVector{Teleobj};kwargs...) where{Teleobj<:Turbine}
    nel        = length(o)
    mut        = (;a=𝕣2(undef,3,3*nel))
    opt        = (Δz        = default{:height   }(kwargs,20     )/2,
                  color     = default{:color    }(kwargs,:orange),
                  linewidth = default{:linewidth}(kwargs,5      ) )
    return mut,opt
end
function Muscade.update_drawing(  axe,o::AbstractVector{Teleobj},mut,opt, Λ,X,U,A,t,SP,dbg) where{Teleobj<:Turbine}
    nel              = length(o)
    for iel          = 1:nel
        i            = 3*(iel-1)
        x            = ∂0(X)[:,iel]+o[iel].xₘ
        mut.a[1,i+1] = mut.a[1,i+2] = x[1] 
        mut.a[2,i+1] = mut.a[2,i+2] = x[2] 
        mut.a[3,i+1] = o[iel].z-opt.Δz
        mut.a[3,i+2] = o[iel].z+opt.Δz
        mut.a[:,i+3].= NaN
    end
    return mut
end
function Muscade.display_drawing!(            axe,::Type{<:Turbine},obs,opt) 
    lines!(axe,obs.a,color=opt.color, linewidth=opt.linewidth)
end

Muscade.doflist( ::Type{<:Turbine}) = (inod =(1   ,1   ,2        ,2        ),
                                       class=(:X  ,:X  ,:A       ,:A       ),
                                       field=(:tx1,:tx2,:Δseadrag,:Δskydrag))

### AnchorLine

struct AnchorLine <: AbstractElement
    xₘtop   :: SVector{3,𝕣}  # x1,x2,x3
    Δxₘtop  :: SVector{3,𝕣}  # as meshed, node to fairlead
    xₘbot   :: SVector{2,𝕣}  # x1,x2 (x3=0)
    L       :: 𝕣
    buoyancy:: 𝕣
end
AnchorLine(nod::Vector{Node};Δxₘtop,xₘbot,L,buoyancy) = AnchorLine(coord(nod)[1],Δxₘtop,xₘbot,L,buoyancy)
         
p = SVector(   2.82040487827,  -24.86027164695,   153.69500343165, -729.52107422849, 2458.11921356871,
              -5856.85610233072, 9769.49700812681,-11141.12651712473, 8260.66447746395,-3582.36704093187,
                687.83550335374)

Muscade.doflist(     ::Type{<:AnchorLine}) = (inod =(1   ,1   ,1   ,2  ,2         ),
                                              class=(:X  ,:X  ,:X  ,:A ,:A        ),
                                              field=(:tx1,:tx2,:rx3,:ΔL,:Δbuoyancy))

@espy function Muscade.lagrangian(o::AnchorLine, Λ,X,U,A,t,SP,dbg)
    xₘtop,Δxₘtop,xₘbot,L,buoyancy = o.xₘtop,o.Δxₘtop,o.xₘbot,o.L*(1+A[1]),o.buoyancy*(1+A[2])      # a for anchor, t for TDP, f for fairlead
    x        = ∂0(X)  
    ☼Xtop    = SVector(x[1],x[2],0.) + xₘtop
    α        =  x[3]                            # azimut from COG to fairlead
    c,s      = cos(α),sin(α)
    ☼ΔXtop   = SMatrix{3,3}(c,s,0,-s,c,0,0,0,1)*Δxₘtop       # arm of the fairlead
    ☼ΔXchain = Xtop[1:2]+ΔXtop[1:2]-xₘbot        # vector from anchor to fairlead
    ☼xaf     = norm(ΔXchain)                     # horizontal distance from anchor to fairlead
    ☼cr      = exp10(horner(p,(L-xaf)/Xtop[3]))*Xtop[3] # curvature radius at TDP
    ☼Fh      = -cr*buoyancy                      # horizontal force
    ☼ltf     = √(Xtop[3]^2+2Xtop[3]*cr)          # horizontal distance from fairlead to TDP
    Fd       = ΔXchain/xaf.*Fh
    m3       = ΔXtop[1]*Fd[2]-ΔXtop[2]*Fd[1]
    L        = Λ[1:2] ∘₁ Fd
    L       += Λ[3  ] *  m3 
    return L,noFB
end

function Muscade.allocate_drawing(axe,o::AbstractVector{E};kwargs...) where{E<:AnchorLine} 
    nel           = length(o)
    req           = @request (Xtop,ΔXtop,ΔXchain,cr,xaf,ltf)
    opt = (
        blue          = default{:blue         }(kwargs,:blue )  ,
        red           = default{:red          }(kwargs,:red  )  ,
        green         = default{:green        }(kwargs,:green)  ,
        linewidth     = default{:linewidth    }(kwargs,2     )  ,
        tdpmarkersize = default{:tdpmarkersize}(kwargs,10    )  ,
        botmarkersize = default{:botmarkersize}(kwargs,20    )  ,
        ncat          = default{:cat          }(kwargs,11    )  ,
        req           = req                                     ,
        ΔXtop         = 𝕣2(undef,3,nel)                         ,
        Xtop          = 𝕣2(undef,3,nel) 
         ) 
    mut = (
        Xbot   = 𝕣2(undef,3,nel)          ,
        Xtdp   = 𝕣2(undef,3,nel)          ,
        Xline  = 𝕣2(undef,3,nel*(opt.ncat+1)) ,
        Xexc   = 𝕣2(undef,3,nel*3)        ,
        Xfloor = 𝕣2(undef,3,nel*3)        
         ) 
    return mut,opt
end
function Muscade.update_drawing(  axe,o::AbstractVector{E},mut,opt, Λ,X,U,A,t,SP,dbg) where{E<:AnchorLine} 
    nel = length(o)
    for iel = 1:nel
        iline = (opt.ncat+1)*(iel-1)
        iexc  = (2   +1)*(iel-1)
        Λᵢ = Λ[    :,iel]
        Xᵢ = ∂0(X)[:,iel]
        Uᵢ = ∂0(U)[:,iel]
        Aᵢ = A[    :,iel]
        L,FB,out = Muscade.lagrangian(o[iel], Λᵢ,(Xᵢ,),(Uᵢ,),Aᵢ, t,SP,(dbg...,espy2draw=true),opt.req)
        Laf,mut.Xbot[1:2,iel],opt.Xtop[:,iel],opt.ΔXtop[:,iel],ΔXchain,cr,xaf,Ltf = o[iel].L, o[iel].xₘbot, out.Xtop,out.ΔXtop,out.ΔXchain, out.cr, out.xaf, out.ltf
        mut.Xbot[3,iel] = 0.
        n     = ΔXchain./xaf  # horizontal normal vector from anchor to fairlead
        xat   = Laf-Ltf
        xtf   = xaf-xat
        mut.Xtdp[1:2,iel]  = mut.Xbot[1:2,iel] + n*xat
        mut.Xtdp[3  ,iel]  = 0.
        ζ                     = range(max(0,-xat),xtf,opt.ncat)
        mut.Xline[1,iline.+(1:opt.ncat)] = n[1].*ζ.+mut.Xtdp[1,iel]
        mut.Xline[2,iline.+(1:opt.ncat)] = n[2].*ζ.+mut.Xtdp[2,iel]
        mut.Xline[3,iline.+(1:opt.ncat)] = cr.*(cosh.(ζ./cr).-1)
        mut.Xline[:,iline +1+ opt.ncat ].= NaN
        mut.Xexc[  :,iexc+1]         = opt.Xtop[:,iel]
        mut.Xexc[  :,iexc+2]         = opt.Xtop[:,iel] + opt.ΔXtop[:,iel]
        mut.Xexc[  :,iexc+3]        .= NaN
        mut.Xfloor[:,iexc+1]         = mut.Xbot[:,iel]
        mut.Xfloor[:,iexc+2]         = mut.Xtdp[:,iel] 
        mut.Xfloor[:,iexc+3]        .= NaN
    end

    return mut
end

function Muscade.display_drawing!(            axe,::Type{<:AnchorLine},obs,opt) 
    lines!(  axe,obs.Xline  ,color=opt.blue , linewidth  = opt.linewidth) # line
    lines!(  axe,obs.Xexc   ,color=opt.red  , linewidth  = opt.linewidth) # excentricity
    lines!(  axe,obs.Xfloor ,color=opt.green, linewidth  = opt.linewidth) # seafloor
    scatter!(axe,obs.Xtdp   ,color=opt.blue , markersize = opt.tdpmarkersize)
    scatter!(axe,obs.Xbot   ,color=opt.red  , markersize = opt.botmarkersize)
end


#### Spring

struct Spring{D} <: AbstractElement
    x₁     :: SVector{D,𝕣}  # x1,x2,x3
    x₂     :: SVector{D,𝕣} 
    EA     :: 𝕣
    L      :: 𝕣
end
Spring{D}(nod::Vector{Node};EA) where{D}= Spring{D}(coord(nod)[1],coord(nod)[2],EA,norm(coord(nod)[1]-coord(nod)[2]))
@espy function Muscade.residual(o::Spring{D}, X,U,A, t,SP,dbg) where{D}
    ☼L₀      = o.L *exp10(A[1]) 
    ☼EA      = o.EA*exp10(A[2]) 
    x₁       = ∂0(X)[SVector{D}(i   for i∈1:D)]+o.x₁
    x₂       = ∂0(X)[SVector{D}(i+D for i∈1:D)]+o.x₂
    Δx       = x₁-x₂
    ☼L       = norm(Δx)
    ☼T       = EA*(L-L₀)/L₀   
    F₁       = Δx/L*T # external force on node 1
    R        = vcat(F₁,-F₁)
    return R,noFB
end
Muscade.doflist(     ::Type{Spring{D}}) where{D}=(
    inod  = (( 1 for i=1: D)...,(2 for i=1:D)...,3,3),
    class = ((:X for i=1:2D)...,:A,:A),
    field = ((Symbol(:tx,i) for i=1: D)...,(Symbol(:tx,i) for i=1: D)...,:ΞL₀,:ΞEI)) # \Xi

Muscade.no_second_order(::Type{<:Spring}) = Val(false)


### SdofOscillator

struct SdofOscillator <: AbstractElement
    K₁ :: 𝕣
    K₂ :: 𝕣
    C₁ :: 𝕣
    C₂ :: 𝕣
    M₁ :: 𝕣
    M₂ :: 𝕣
end
SdofOscillator(nod::Vector{Node};K₁=0.::𝕣,K₂::𝕣=0.,C₁=0.::𝕣,C₂::𝕣=0.,M₁=0.::𝕣,M₂::𝕣=0.) = SdofOscillator(K₁,K₂,C₁,C₂,M₁,M₂)
@espy function Muscade.residual(o::SdofOscillator, X,U,A, t,SP,dbg) 
    x,x′,x″,u = ∂0(X)[1], ∂1(X)[1], ∂2(X)[1], ∂0(U)[1]
    R         = SVector(-u +o.K₁*x +o.K₂*x^2  +o.C₁*x′ +o.C₂*x′^2 +o.M₁*x″ +o.M₂*x″^2)
    return R,noFB
end
Muscade.doflist( ::Type{SdofOscillator})  = (inod =(1 ,1 ), class=(:X,:U), field=(:tx1,:tu1))

### AdjustableSdofOscillator

struct AdjustableSdofOscillator <: AbstractElement
    K :: 𝕣
    C :: 𝕣
    M :: 𝕣
end
AdjustableSdofOscillator(nod::Vector{Node};K=1.::𝕣,C=0.::𝕣,M=0.::𝕣) = AdjustableSdofOscillator(K,C,M)
@espy function Muscade.residual(o::AdjustableSdofOscillator, X,U,A, t,SP,dbg) 
    x,x′,x″,u = ∂0(X)[1], ∂1(X)[1], ∂2(X)[1], ∂0(U)[1]
    # ☼K        = o.K *exp10(A[1]) 
    # R         = SVector(-u +K*x +o.C*x′ +o.M*x″)
    ☼C        = o.C *exp10(A[1]) 
    R         = SVector(-u +o.K*x +C*x′ +o.M*x″)
    return R,noFB
end
#Muscade.doflist( ::Type{AdjustableSdofOscillator})  = (inod =(1 ,1, 1, 1), class=(:X,:U,:A,:A), field=(:tx1,:tu1,:ΞK,:ΞC))
Muscade.doflist( ::Type{AdjustableSdofOscillator})  = (inod =(1 ,1, 1), class=(:X,:U,:A), field=(:tx1,:tu1,:ΞK))
