#module TestSweepXA2
using Test
using Muscade


include("SomeElements.jl")
include("../examples/DryFriction.jl")

Δt              = 0.01
t               = Δt:Δt:1000*Δt

### Target

K,C,M           = 1.,.4,.3

model           = Model(:Target)
node            = addnode!(model,𝕣[])
osc             = addelement!(model,AdjustableSdofOscillator,[node]; K,C,M)

initialstate    = initialize!(model;time=0.)
x,x′,x″         = 0.,1.,0.   
initialstate    = setdof!(initialstate,[x′];field=:tx1,nodID=[node],order=1)  # initial speed

state0          = solve(SweepX{ 2};  initialstate,time= t,verbose=false,catcherror=true)
target          = Muscade.FunctionFromVector(t,getdof(state0;field=:tx1,nodID=[node]))

### Optimize

K,C,M           = 1.,.6,.3

model           = Model(:Optim)
node            = addnode!(model,𝕣[])
osc             = addelement!(model,AdjustableSdofOscillator,[node]; K,C,M)

@functor with(σ=2.        ) acost(a  )=(a/σ)^2
@functor with(σ=.01,target) xcost(x,t)=((x-target(t))/σ)^2
cK              = addelement!(model,SingleAcost  ,[node];field=:ΞK,               cost=acost)
cX              = addelement!(model,SingleDofCost,[node];class=:X ,field=:tx1    ,cost=xcost)

initialstate    = initialize!(model;time=0.)
x,x′,x″         = 0.,1.,0.   
initialstate    = setdof!(initialstate,[x′];field=:tx1,nodID=[node],order=1)  # initial speed

state  = solve(SweepXA{2,2,2};  initialstate,time= t,verbose=true,catcherror=true,maxAiter=20,maxΔa=1e-10)

# @show state[1].A[1]

# using GLMakie
# fig      = Figure(size = (1000,800))
# axeX      = Axis(fig[1,1])
# x0 = [s.X[1][1] for s∈state0]
# x  = [s.X[1][1] for s∈state]
# lines!(axeX,t,x0,color=:red)
# lines!(axeX,t,x.+.01,color=:black)
# display(fig)


;
#end