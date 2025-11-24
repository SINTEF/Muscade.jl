#module TestSweepXA2
using Test
using Muscade


include("SomeElements.jl")
include("../examples/DryFriction.jl")


K,C,M           = 1.,.4,.3

model           = Model(:TestModel)
node            = addnode!(model,𝕣[])
osc             = addelement!(model,AdjustableSdofOscillator,[node]; K,C,M)

@functor with(σ=2.) acost(a  )=(a/σ)^2
@functor with(σ=.1) xcost(x,t)=(x/σ)^2
cK              = addelement!(model,SingleAcost  ,[node];field=:ΞK,               cost=acost)
cX              = addelement!(model,SingleDofCost,[node];class=:X ,field=:tx1    ,cost=xcost)

initialstate    = initialize!(model;time=0.)
x,x′,x″         = 0.,1.,0.   
initialstate    = setdof!(initialstate,[x′];field=:tx1,nodID=[node],order=1)  # initial speed


Δt    = 0.1
t     = Δt:Δt:100*Δt
state0 = solve(SweepX{2};  initialstate,time= t,verbose=true,catcherror=true)
state = solve(SweepXA{2};  initialstate,time= t,verbose=true,catcherror=true,maxAiter=20,maxΔa=1e-10)

@show state[1].A[1]

using GLMakie
fig      = Figure(size = (1000,800))
axeX      = Axis(fig[1,1])
x0 = [s.X[1][1] for s∈state0]
x = [s.X[1][1] for s∈state]
lines!(axeX,t,x,color=:black)
lines!(axeX,t,x0,color=:red)
display(fig)



;
#end