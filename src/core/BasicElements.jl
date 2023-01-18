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

S(  λ,g,γ) = (g+λ    -hypot(g-λ,2γ))/2 # Modified interior point method's take on KKT's-complementary slackness 
S∂g(λ,g,γ) = (1-(g-λ)/hypot(g-λ,2γ))/2

#KKT(λ        ,g         ,γ::𝕣,λₛ,gₛ)                 = λ*g # A pseudo-potential with strange derivatives
KKT(λ::𝕣        ,g::𝕣         ,γ::𝕣,λₛ,gₛ)                 = λ*g # A pseudo-potential with strange derivatives
KKT(λ::∂ℝ{P,N,R},g::∂ℝ{P,N,R},γ::𝕣,λₛ,gₛ) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(0, S∂g(λ.x/λₛ,g.x/gₛ,γ)*λ.x*g.dx + gₛ*S(λ.x/λₛ,g.x/gₛ,γ)*λ.dx)
KKT(λ:: ℝ       ,g::∂ℝ{P,N,R},γ::𝕣,λₛ,gₛ) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(0, S∂g(λ.x/λₛ,g.x/gₛ,γ)*λ.x*g.dx                           )
KKT(λ::∂ℝ{P,N,R},g:: ℝ       ,γ::𝕣,λₛ,gₛ) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(0,                                gₛ*S(λ.x/λₛ,g.x/gₛ,γ)*λ.dx)
function KKT(λ::∂ℝ{Pλ,Nλ,Rλ},g::∂ℝ{Pg,Ng,Rg},γ::𝕣,λₛ,gₛ) where{Pλ,Pg,Nλ,Ng,Rλ<:ℝ,Rg<:ℝ}
    if Pλ==Pg
        R = promote_type(Rλ,Rg)
        return ∂ℝ{Pλ,Nλ}(convert(R,KKT(λ.x,g.x,γ,λₛ,gₛ)),convert.(R,     S∂g(λ.x/λₛ,g.x/gₛ,γ)*λ.x*g.dx + gₛ*S(λ.x/λₛ,g.x/gₛ,γ)*λ.dx))
    elseif Pλ> Pg
        R = promote_type(Rλ,typeof(b))
        return ∂ℝ{Pλ,Nλ}(convert(R,KKT(λ  ,g.x,γ,λₛ,gₛ)),convert.(R,     S∂g(λ.x/λₛ,g.x/gₛ,γ)*λ.x*g.dx                            ))
    else
        R = promote_type(typeof(a),Rg)
        return ∂ℝ{Pg,Ng}(convert(R,KKT(λ.x,g  ,γ,λₛ,gₛ)),convert.(R,                                    gₛ*S(λ.x/λₛ,g.x/gₛ,γ)*λ.dx))
    end
end

#-------------------------------------------------

struct Constraint{λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield,Tg,Tkind} <: AbstractElement
    g        :: Tg    # g(x,t) for Xconstraints, or g(x,u,a,t) otherwise
    kind     :: Tkind # kind(t)->symbol, or Symbol for Aconstraints
    gₛ        :: 𝕣
    λₛ        :: 𝕣  
end
Constraint{    λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield                       }(g,kind,gₛ,λₛ) where
              {λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield} =
    Constraint{λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield,typeof(g),typeof(kind)}(g,kind,gₛ,λₛ)
Constraint(nod::Vector{Node};xinod::NTuple{Nx,𝕫},xfield::NTuple{Nx,Symbol},
                                         uinod::NTuple{Nu,𝕫},ufield::NTuple{Nu,Symbol},
                                         ainod::NTuple{Na,𝕫},afield::NTuple{Na,Symbol},
                                         λinod::𝕫, λclass::Symbol, λfield,
                                         gₛ::𝕣=1.,λₛ::𝕣=1.,
                                         g::Function ,kind::Function) where{Nx,Nu,Na} =
                 Constraint{λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield}(g,kind,gₛ,λₛ)
doflist(::Type{<:Constraint{λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield}}) where
                           {λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield} = 
   (inod =(xinod...           ,uinod...           ,ainod...           ,λinod         ), 
    class=(ntuple(i->:X,Nx)...,ntuple(i->:U,Nu)...,ntuple(i->:A,Na)...,Symbol(λclass)), 
    field=(xfield...          ,ufield...          ,afield...          ,λfield        )) 

off,equal,inequal = :off,:equal,:inequal # because @espy has its own ways with symbols... TODO improve @espy
@espy function residual(o::Constraint{Xclass,Nx,0,0}, X,U,A, t,γ,dbg) where{Nx}
    P,gₛ,λₛ     = constants(∂0(X)),o.gₛ,o.λₛ
    x,λ        = ∂0(X)[SVector{Nx}(1:Nx)], ∂0(X)[Nx+1]
    x∂         = variate{P,Nx}(x) 
    g,g∂x      = value_∂{P,Nx}(o.g(x∂,t)) 
    return if o.kind(t)==off;     SVector{Nx+1}(         ntuple(i->0,Nx)...,-gₛ/λₛ*λ         ) 
    elseif    o.kind(t)==equal;   SVector{Nx+1}((                -g∂x*λ)...,-     g         )
    elseif    o.kind(t)==inequal; SVector{Nx+1}((-S∂g(λ/λₛ,g/gₛ,γ)*g∂x*λ)...,-gₛ*S(λ/λₛ,g/gₛ,γ)) 
    else MuscadeException("kind(t) must have value :off, :equal or :inequal",dbg)
    end
end
@espy function lagrangian(o::Constraint{class,Nx,Nu,Na}, δX,X,U,A, t,γ,dbg) where{class<:Union{Uclass,Aclass},Nx,Nu,Na}
    if class==Uclass; x,u,a,λ = ∂0(X),∂0(U)[SVector{Nu}(1:Nu)],A                   ,∂0(U)[Nu+1] end
    if class==Aclass; x,u,a,λ = ∂0(X),∂0(U)                   ,A[SVector{Na}(1:Na)],A[    Na+1] end
    g = o.g(x,u,a,t)
    return if o.kind(t)==off;     -o.gₛ/(2o.λₛ)*λ^2 
    elseif    o.kind(t)==equal;   -g*λ
    elseif    o.kind(t)==inequal; -KKT(λ,g,γ,o.λₛ,o.gₛ) 
    else MuscadeException("kind(t) must have value :off, :equal or :inequal",dbg)
    end
end

#-------------------------------------------------

struct Hold <: AbstractElement end  
# id1(v,t) = v[1]
# eq(t)    = :equal
Hold(nod::Vector{Node};field::Symbol,λfield::Symbol=Symbol(:λ,field)) = 
    Constraint{Xclass,1, 0, 0, (1,),(field,),(),   (),    (),   (),    1,    λfield}((v,t)->v[1] , t->:equal,1.,1.)
#   Constraint{λclass,Nx,Nu,Na,xinod,xfield, uinod,ufield,ainod,afield,λinod,λfield}

#-------------------------------------------------
