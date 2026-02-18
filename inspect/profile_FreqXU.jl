#module TestFreqXU

# cd("C:\\Users\\philippem\\.julia\\dev\\Muscade")
# using Pkg 
# Pkg.activate(".")



using Test
using Muscade
using StaticArrays,SparseArrays

###

using StaticArrays
struct El1 <: AbstractElement
    K :: 𝕣
    C :: 𝕣
    M :: 𝕣
end
El1(nod::Vector{Node};K::𝕣,C::𝕣,M::𝕣) = El1(K,C,M)
@espy function Muscade.residual(o::El1, X,U,A, t,SP,dbg) 
    x,x′,x″,u,ΞC,ΞM = ∂0(X)[1], ∂1(X)[1], ∂2(X)[1], ∂0(U)[1], A[1], A[2]
    r         = -u + o.K*x + o.C*exp10(ΞC)*x′ + o.M*exp10(ΞM)*x″
    return SVector(r),noFB
end
Muscade.doflist( ::Type{El1})  = (inod =(1 ,1 ,1 ,1), class=(:X,:U,:A,:A), field=(:tx1,:utx1,:ΞC,:ΞM))

include("../test/SomeElements.jl")

model           = Model(:TrueModel)
n1              = addnode!(model,𝕣[0])  
n2              = addnode!(model,𝕣[1])  
n3              = addnode!(model,𝕣[]) # anode for spring
e1              = addelement!(model,El1,[n1], K=1.,C=0.05,M=1.)
e2              = addelement!(model,El1,[n2], K=0.,C=0.0 ,M=1.)
e3              = addelement!(model,Spring{1},[n1,n2,n3], EI=1.1)
@functor with() f(x) = x^2
@functor with() fu(u,t) = u^2
@functor with() l1(tx1,t) = (tx1-0.1*sin(t))^2
@functor with() l2(tx1,t) = (tx1-0.1*cos(t))^2
e4              = addelement!(model,SingleDofCost,[n3];class=:A,field=:ΞL₀,     cost=f)
e5              = addelement!(model,SingleDofCost,[n3];class=:A,field=:ΞEI,     cost=f)
e6              = addelement!(model,SingleDofCost,[n1];class=:A,field=:ΞC ,     cost=f)
e7              = addelement!(model,SingleDofCost,[n1];class=:A,field=:ΞM ,     cost=f)
e8              = addelement!(model,SingleDofCost,[n2];class=:A,field=:ΞC ,     cost=f)
e9              = addelement!(model,SingleDofCost,[n2];class=:A,field=:ΞM ,     cost=f)
e10             = addelement!(model,SingleUdof   ,[n1];Xfield=:tx1,Ufield=:utx1,cost=fu)
e11             = addelement!(model,SingleUdof   ,[n2];Xfield=:tx1,Ufield=:utx1,cost=fu)
e12             = addelement!(model,SingleDofCost,[n1];class=:X,field=:tx1,     cost=l1)
e13             = addelement!(model,SingleDofCost,[n2];class=:X,field=:tx1,     cost=l2)

state0          = initialize!(model)   

OX               = 2
OU               = 0


using Profile,ProfileView,BenchmarkTools
mission = :profile
if mission == :report
    stateXUA         = solve(FreqXU{OX,OU};Δt=1., p=8, t₀=0.,initialstate=state0)
elseif mission == :time
    @btime state    = stateXUA         = solve(FreqXU{OX,OU};Δt=1., p=8, t₀=0.,initialstate=state0,verbose=false)
elseif mission == :profile
    local stateXUA         = solve(FreqXU{OX,OU};Δt=1., p=16, t₀=0.,initialstate=state0,verbose=false)
    Profile.clear()
    Profile.@profile for i=1:10
        local stateXUA         = solve(FreqXU{OX,OU};Δt=1., p=16, t₀=0.,initialstate=state0,verbose=false)
    end
    ProfileView.view(fontsize=30);
end
# After clicking on a bar in the flame diagram, you can type warntype_last() and see the result of 
# code_warntype for the call represented by that bar.

;
