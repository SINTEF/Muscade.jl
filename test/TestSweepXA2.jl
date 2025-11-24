#module TestSweepXA2
using Test
using Muscade


include("SomeElements.jl")
include("../examples/DryFriction.jl")


K,C,M           = 1.,.4,.3

model           = Model(:TestModel)
node            = addnode!(model,𝕣[])
osc             = addelement!(model,AdjustableSdofOscillator,[node]; K,C,M)

@functor with(    ) acost(a,σ)=(a/σ)^2
@functor with(σ=.1) xcost(x,t)=(x/σ)^2
cK              = addelement!(model,SingleAcost  ,[node];field=:ΞK,costargs=(2.02,),cost=acost)
cX              = addelement!(model,SingleDofCost,[node];class=:X ,field=:tx1    ,cost=xcost)

initialstate    = initialize!(model;time=0.)
x,x′,x″         = 0.,1.,0.   
initialstate    = setdof!(initialstate,[x′];field=:tx1,nodID=[node],order=2)  # initial speed



Δt    = 0.01
t     = Δt:Δt:40*Δt
#t     = Δt:Δt:5*Δt
stateSweep = solve(SweepXA{2};  initialstate,time= t,verbose=true,catcherror=true,maxAiter=1,maxΔa=1e-10)

# stateDirect = solve(DirectXUA{2,0,1};initialstate=[initialstate],time= [t])
# fig      = Figure(size = (1000,800))
# axeX      = Axis(fig[1,1])
# #axeA      = Axis(fig[2,1])
# axeΛ      = Axis(fig[2,1])
# x = [s.X[1][1] for s∈stateDirect[1]]

;
#end