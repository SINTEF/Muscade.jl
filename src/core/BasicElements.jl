struct XdofCost{Tcost,Field,Derivative} <: AbstractElement
    cost :: Tcost # Function 
end
XdofCost(nod::Vector{Node};field::Symbol,cost::Tcost,derivative=0::𝕫) where{Tcost<:Function} = XdofCost{Tcost,field,derivative}(cost)
Muscade.doflist(::Type{XdofCost{Tcost,Field,Derivative}}) where{Tcost,Field,Derivative} = (inod =(1,), class=(:X,), field=(Field,))
@espy function Muscade.lagrangian(o::XdofCost{Tcost,Field,Derivative}, δX,X,U,A, t,ε,dbg) where{Tcost,Field,Derivative}
    :J = o.cost(∂n(X,Derivative)[1])
    return J
end
Muscade.espyable(::Type{<:XdofCost}) = (J=scalar,)

#-------------------------------------------------
# TODO extend this element (or create a new one "UdofLoad", that creates a U-dof acting on a X-dof and associates a cost to the Udof)
struct UdofCost{Tcost,Field,Derivative} <: AbstractElement
    cost :: Tcost # Function 
end
UdofCost(nod::Vector{Node};field::Symbol,cost::Tcost,derivative=0::𝕫) where{Tcost<:Function} = UdofCost{Tcost,field,derivative}(cost)
Muscade.doflist(::Type{UdofCost{Tcost,Field,Derivative}}) where{Tcost,Field,Derivative} = (inod =(1,), class=(:U,), field=(Field,))
@espy function Muscade.lagrangian(o::UdofCost{Tcost,Field,Derivative}, δX,X,U,A, t,ε,dbg) where{Tcost,Field,Derivative}
    :J = o.cost(∂n(XU,Derivative)[1])
    return J
end
Muscade.espyable(::Type{<:UdofCost}) = (J=scalar,)

#-------------------------------------------------

struct AdofCost{Tcost,Field} <: AbstractElement
    cost :: Tcost # Function 
end
AdofCost(nod::Vector{Node};field::Symbol,cost::Tcost) where{Tcost<:Function} = AdofCost{Tcost,field}(cost)
Muscade.doflist(::Type{AdofCost{Tcost,Field}}) where{Tcost,Field} = (inod=(1,), class=(:A,), field=(Field,))
@espy function Muscade.lagrangian(o::AdofCost{Tcost,Field}, δX,X,U,A, t,ε,dbg) where{Tcost,Field}
    :J = o.cost(A[1])
    return J
end
Muscade.espyable(::Type{<:AdofCost}) = (J=scalar,)

#-------------------------------------------------

struct DofLoad{Tvalue,Field} <: AbstractElement
    value      :: Tvalue # Function
end
DofLoad(nod::Vector{Node};field::Symbol,value::Tvalue) where{Tvalue<:Function} = DofLoad{Tvalue,field}(value)
Muscade.doflist(::Type{DofLoad{Tvalue,Field}}) where{Tvalue,Field}=(inod=(1,), class=(:X,), field=(Field,))
@espy function Muscade.residual(o::DofLoad, X,U,A, t,ε,dbg) 
    :F = o.value(t)
    return SVector{1}(-F)
end
Muscade.espyable(::Type{<:DofLoad}) = (F=scalar,)

#-------------------------------------------------

slack(g,λ,γ) = (g+λ)/2-hypot(γ,(g-λ)/2) # Modified interior point method's take on KKT's-complementary slackness 

KKT(      g,λ,γ) = g*λ # A pseudo potential
KKT(g::∂ℝ{P,N,R},λ::∂ℝ{P,N,R},γ::𝕣) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(KKT(g.x,λ.x,γ) , λ.x*g.dx + slack(g.x,λ.x,γ)*λ.dx)
KKT(g::∂ℝ{P,N,R},λ:: ℝ       ,γ::𝕣) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(KKT(g.x,λ  ,γ) , λ  *g.dx                            )
KKT(g:: ℝ       ,λ::∂ℝ{P,N,R},γ::𝕣) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(KKT(g  ,λ.x,γ) ,            slack(g  ,λ.x,γ)*λ.dx)
function KKT(a::∂ℝ{Pa,Na,Ra},b::∂ℝ{Pb,Nb,Rb},γ::𝕣) where{Pa,Pb,Na,Nb,Ra<:ℝ,Rb<:ℝ}
    if Pa==Pb
        R = promote_type(Ra,Rb)
        return ∂ℝ{Pa,Na}(convert(R,KKT(g.x,λ.x,γ)),convert.(R,λ.x*g.dx + slack(g.x,λ.x,γ)*λ.dx))
    elseif Pa> Pb
        R = promote_type(Ra,typeof(b))
        return ∂ℝ{Pa,Na}(convert(R,KKT(g.x,λ  ,γ)),convert.(R,λ  *g.dx                            ))
    else
        R = promote_type(typeof(a),Rb)
        return ∂ℝ{Pb,Nb}(convert(R,KKT(g  ,λ.x,γ)),convert.(R,            slack(g  ,λ.x,γ)*λ.dx))
    end
end

#-------------------------------------------------

struct HolonomicConstraint{N,xinod,xfield,λinod,λfield,Tg} <: AbstractElement
    g :: Tg # Function
end
HolonomicConstraint(nod::Vector{Node};xinod::NTuple{N,𝕫},xfield::NTuple{N,Symbol},
                                      λinod::         𝕫 ,λfield::         Symbol ,
                                      g::Function) where{N} =
    HolonomicConstraint{N,xinod,xfield,λinod,λfield,typeof(g)}(g)
Muscade.doflist(::Type{<:HolonomicConstraint{N,xinod,xfield,λinod,λfield}}) where{N,xinod,xfield,λinod,λfield} = 
   (inod=(xinod...,λinod), class=ntuple(i->:X,N+1), field=(xfield...,λfield)) 
@espy function Muscade.lagrangian(o::HolonomicConstraint{N}, δX,X,U,A, t,ε,dbg) where{N}
    P          = constants(δX,∂0(X))
    X∂         = directional{P}(∂0(X),δX) 
    x ,λ       = X∂[SVector{N}(1:N)], -X∂[N+1] 
    return ∂{P}(o.g(x)*λ)    # = δ(gλ) = δg*λ+δλ*g = δx∘∇ₓg*λ+δλ*g
end

@espy function Muscade.residual(o::HolonomicConstraint{N}, X,U,A, t,ε,dbg) where{N}
    P          = constants(∂0(X))
    x,λ        = ∂0(X)[SVector{N}(1:N)], -∂0(X)[N+1]
    x∂         = variate{P,N}(x) 
    g,∇ₓg      = value_∂{P,N}(o.g(x∂)) 
    return  SVector{N+1}(∇ₓg*λ...,g)
end

# directional{3, 2}(::Tuple{SVector{2, ∂ℝ{2, 4, ∂ℝ{1, 4, Float64}}}}, 
#                   ::      SVector{2, ∂ℝ{2, 4, ∂ℝ{1, 4, Float64}}})


id1(v) = v[1]
struct DofHold <: AbstractElement end  
DofHold(nod::Vector{Node};field::Symbol,λfield::Symbol=Symbol(:λ,field)) = HolonomicConstraint{1,(1,),(field,),1,λfield,typeof(id1)}(id1)

#-------------------------------------------------

# struct DofHold{Field,λfield} <: AbstractElement
# end
# DofHold(nod::Vector{Node};field::Symbol,λfield::Symbol=Symbol(:λ,field))  = DofHold{field,λfield}()
# Muscade.doflist(::Type{DofHold{Field,λfield}}) where{Field,λfield}=(inod=(1,1), class=(:X,:X), field=(Field,λfield))
# @espy function Muscade.residual(o::DofHold, X,U,A, t,ε,dbg) 
#     x,:λ       = ∂0(X)[1],∂0(X)[2] # it's +∂0(X)[2]: "internal force" λ will be negative for x in the positive direction
#     return SVector{2}(-λ,x)
# end
# Muscade.espyable(::Type{<:DofHold}) = (λ=scalar,)


#-------------------------------------------------

struct Spring{D} <: AbstractElement
    x₁     :: SVector{D,𝕣}  # x1,x2,x3
    x₂     :: SVector{D,𝕣} 
    EI     :: 𝕣
    L      :: 𝕣
end
Spring{D}(nod::Vector{Node};EI) where{D}= Spring{D}(coord(nod)[1],coord(nod)[2],EI,norm(coord(nod)[1]-coord(nod)[2]))
@espy function Muscade.residual(o::Spring{D}, X,U,A, t,ε,dbg) where{D}
    x₁       = ∂0(X)[SVector{D}(i   for i∈1:D)]+o.x₁
    x₂       = ∂0(X)[SVector{D}(i+D for i∈1:D)]+o.x₂
    :L₀      = o.L *exp10(A[1]) 
    :EI      = o.EI*exp10(A[2]) 
    Δx       = x₁-x₂
    :L       = norm(Δx)
    :T       = EI*(L-L₀)
    F₁       = Δx/L*T # external force on node 1
    R        = vcat(F₁,-F₁)
    return R
end
Muscade.doflist(     ::Type{Spring{D}}) where{D}=(
    inod  = (( 1 for i=1: D)...,(2 for i=1:D)...,3,3),
    class = ((:X for i=1:2D)...,:A,:A),
    field = ((Symbol(:tx,i) for i=1: D)...,(Symbol(:tx,i) for i=1: D)...,:ΞL₀,:ΞEI)) # \Xi
Muscade.espyable(    ::Type{<:Spring}) = (EI=scalar,L₀=scalar,L=scalar,T=scalar)


