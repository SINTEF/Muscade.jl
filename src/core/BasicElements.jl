struct Xclass end
struct Uclass end
struct Aclass end
Base.Symbol(::Type{Xclass}) = :X
Base.Symbol(::Type{Uclass}) = :U
Base.Symbol(::Type{Aclass}) = :A

#-------------------------------------------------

struct DofCost{Derivative,Class,Field,Tcost} <: AbstractElement
    cost :: Tcost # Function 
end
DofCost(nod::Vector{Node};class::DataType,field::Symbol,cost::Tcost,derivative=0::𝕫) where{Tcost<:Function} = DofCost{derivative,class,field,Tcost}(cost)
doflist(::Type{<:DofCost{Derivative,Class,Field}}) where{Derivative,Class,Field} = (inod =(1,), class=(Symbol(Class),), field=(Field,))
espyable(::Type{<:DofCost}) = (J=scalar,)
@espy function lagrangian(o::DofCost{Derivative,Xclass}, δX,X,U,A, t,γ,dbg) where{Derivative}
    :J = o.cost(∂n(X,Derivative)[1],t)
    return J
end
@espy function lagrangian(o::DofCost{Derivative,Uclass}, δX,X,U,A, t,γ,dbg) where{Derivative}
    :J = o.cost(∂n(U,Derivative)[1],t)
    return J
end
@espy function lagrangian(o::DofCost{Derivative,Aclass}, δX,X,U,A, t,γ,dbg) where{Derivative}
    :J = o.cost(A[1])
    return J
end

#-------------------------------------------------

struct DofLoad{Tvalue,Field} <: AbstractElement
    value      :: Tvalue # Function
end
DofLoad(nod::Vector{Node};field::Symbol,value::Tvalue) where{Tvalue<:Function} = DofLoad{Tvalue,field}(value)
doflist(::Type{DofLoad{Tvalue,Field}}) where{Tvalue,Field}=(inod=(1,), class=(:X,), field=(Field,))
@espy function residual(o::DofLoad, X,U,A, t,γ,dbg) 
    :F = o.value(t)
    return SVector{1}(-F)
end
espyable(::Type{<:DofLoad}) = (F=scalar,)

#-------------------------------------------------

slack(g,λ,γ) = (g+λ)/2-hypot(γ,(g-λ)/2) # Modified interior point method's take on KKT's-complementary slackness 

KKT(      g,λ,γ) = g*λ # A pseudo-potential with strange derivatives
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

struct Constraint{λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield,Tg,Tkind} <: AbstractElement
    g        :: Tg    # g(x,t) for Xconstraints, or g(x,u,a,t) otherwise
    kind     :: Tkind # kind(t)->symbol, or Symbol for Aconstraints
end
Constraint{    λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield                       }(g,kind) where{λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield} =
    Constraint{λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield,typeof(g),typeof(kind)}(g,kind)
Constraint(nod::Vector{Node};xinod::NTuple{Nx,𝕫},xfield::NTuple{Nx,Symbol},
                                         uinod::NTuple{Nu,𝕫},ufield::NTuple{Nu,Symbol},
                                         ainod::NTuple{Na,𝕫},afield::NTuple{Na,Symbol},
                                         λinod::𝕫, λclass::Symbol, λfield,
                                         g::Function ,kind::Function) where{Nx,Nu,Na} =
                 Constraint{λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield}(g,kind)
doflist(::Type{<:Constraint{λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield}}) where{λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield} = 
   (inod =(xinod...           ,uinod...           ,ainod...           ,λinod         ), 
    class=(ntuple(i->:X,Nx)...,ntuple(i->:U,Nu)...,ntuple(i->:A,Na)...,Symbol(λclass)), 
    field=(xfield...          ,ufield...          ,afield...          ,λfield        )) 

off,equal,inequal = :off,:equal,:inequal # because @espy has its own ways with symbols...
@espy function residual(o::Constraint{Xclass,Nx,0,0}, X,U,A, t,γ,dbg) where{Nx}
    P          = constants(∂0(X))
    x,λ        = ∂0(X)[SVector{Nx}(1:Nx)], -∂0(X)[Nx+1]
    x∂         = variate{P,Nx}(x) 
    g,∇ₓg      = value_∂{P,Nx}(o.g(x∂,t)) 
    z = if o.kind(t)==off;     λ # what to set to zero
    elseif o.kind(t)==equal;   g
    elseif o.kind(t)==inequal; slack(g,λ,γ) 
    else MuscadeException("kind(t) must have value :off, :equal or :inequal",dbg)
    end
    return  SVector{Nx+1}((∇ₓg*λ)...,z)
end
@espy function lagrangian(o::Constraint{Xclass,Nx,0,0}, δX,X,U,A, t,γ,dbg) where{Nx}  
    P          = constants(δX,∂0(X))
    X∂         = directional{P}(∂0(X),δX) 
    x,λ        = X∂[SVector{Nx}(1:Nx)], -X∂[Nx+1] 
    g          = o.g(x,t)
    m = if o.kind(t)==off;     λ^2 # what to minimize
    elseif o.kind(t)==equal;   g*λ
    elseif o.kind(t)==inequal; KKT(g,λ,γ) 
    else MuscadeException("kind(t) must have value :off, :equal or :inequal",dbg)
    end
    return ∂{P}(m)    # = δ(gλ) = δg*λ+δλ*g = δx∘∇ₓg*λ+δλ*g   
end
@espy function lagrangian(o::Constraint{Uclass,Nx,Nu,Na}, δX,X,U,A, t,γ,dbg) where{Nx,Nu,Na}
    x,u,a,λ    = ∂0(X),∂0(U)[SVector{Nu}(1:Nu)],A,-∂0(U)[Nu+1] 
    g          = o.g(x,u,a,t)
    m = if o.kind(t)==off;     λ^2 # what to minimize
    elseif o.kind(t)==equal;   g*λ
    elseif o.kind(t)==inequal; KKT(g,λ,γ) 
    else MuscadeException("kind(t) must have value :off, :equal or :inequal",dbg)
    end
    return m 
end
@espy function lagrangian(o::Constraint{Aclass,Nx,Nu,Na}, δX,X,U,A, t,γ,dbg) where{Nx,Nu,Na}
    x,u,a,λ    = ∂0(X),∂0(U),A[SVector{Nu}(1:Nu)],-A[Nu+1] 
    g          = o.g(x,u,a,t)
    m = if o.kind==off;     λ^2 # what to minimize
    elseif o.kind==equal;   g*λ
    elseif o.kind==inequal; KKT(g,λ,γ) 
    else MuscadeException("kind must have value :off, :equal or :inequal",dbg)
    end
    return m 
end

#-------------------------------------------------

struct Hold <: AbstractElement end  
# id1(v,t) = v[1]
# eq(t)    = :equal
Hold(nod::Vector{Node};field::Symbol,λfield::Symbol=Symbol(:λ,field)) = 
    Constraint{Xclass,1, 0, 0, (1,),(field,),(),   (),    (),   (),    1,    λfield}((v,t)->v[1] , t->:equal)
#   Constraint{λclass,Nx,Nu,Na,xinod,xfield, uinod,ufield,ainod,afield,λinod,λfield}

#-------------------------------------------------
