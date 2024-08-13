#module TestDirectXUA

using Test
using Muscade

###

using StaticArrays
struct SdofOscillator <: AbstractElement
    K :: 𝕣
    C :: 𝕣
    M :: 𝕣
end
SdofOscillator(nod::Vector{Node};K::𝕣,C::𝕣,M::𝕣) = SdofOscillator(K,C,M)
@espy function Muscade.residual(o::SdofOscillator, X,U,A, t,SP,dbg) 
    x,x′,x″,u,ΞC,ΞM = ∂0(X)[1], ∂1(X)[1], ∂2(X)[1], ∂0(U)[1], A[1], A[2]
    r         = -u + o.K*x + o.C*exp10(ΞC)*x′ + o.M*exp10(ΞM)*x″
    return SVector(r),noFB
end
Muscade.doflist( ::Type{SdofOscillator})  = (inod =(1 ,1 ,1 ,1), class=(:X,:U,:A,:A), field=(:x,:u,:ΞC,:ΞM))

###

model1          = Model(:TrueModel)
n               = addnode!(model1,𝕣[ 0, 0])  
e               = addelement!(model1,SdofOscillator,[n], K=1.,C=0.05,M=1.)
state0          = initialize!(model1;nXder=2,time=0.)  # make space for 1st-order derivatives, 
# setdof!(state0,[1.];field=:x,nodID=[n],order=1)
# time            = 1.:.1:100
# state1          = solve(SweepX{2};  initialstate=state0,time,verbose=false)
# x               = getdof(state1;field=:x,nodID=[n],order=0 )
# x               = reshape(x,length(x))
# using GLMakie
# fig             = Figure(size = (1000,800))
# axe             = Axis(fig[1,1],title="Test",xlabel="time",ylabel="x")
# oedge           = lines!(  axe,time,x , linewidth = 1)

###
include("../src/DirectXUA.jl")
dis             = state0.dis
out,asm,Ydofgr,Adofgr = prepare(AssemblyDirect    ,model1,dis,3)
#out2,asm2             = prepare(AssemblyDirectLine)


stateXUA           = solve(DirectXUA;initialstate=state0)

#end 

;