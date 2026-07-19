module TestSweepXA2
using Test
using Muscade


include("SomeElements.jl")
include("../examples/DryFriction.jl")

# Δt              = 0.01
# t               = Δt:Δt:1000*Δt
Δt              = 1.
t               = Δt:Δt:10*Δt

### Target

K,C,M           = 1.,.4,.3

model           = Model(:Target)
node            = addnode!(model,𝕣[])
osc             = addelement!(model,AdjustableSdofOscillator,[node]; K,C,M)

initialstate    = initialize!(model;time=0.)
x,x′,x″         = 0.,1.,0.   
initialstate    = setdof!(initialstate,[x′];field=:tx1,nodID=[node],order=1)  # initial speed

state0          = solve(SweepX{ 2};  initialstate,time= t,verbose=false,catcherror=true)
xtarget         = getdof(state0;field=:tx1,nodID=[node])
target          = Muscade.FunctionFromVector(t,xtarget)

### Optimize

K,C,M           = 1.,.6,.3

model           = Model(:Optim)
node            = addnode!(model,𝕣[])
osc             = addelement!(model,AdjustableSdofOscillator,[node]; K,C,M)

@functor with(σ=10.        ) acost(a  )=( a           /σ)^2
@functor with(σ=.0001,target) xcost(x,t)=((x-target(t))/σ)^2
cK              = addelement!(model,SingleAcost  ,[node];field=:ΞC,               cost=acost)
cX              = addelement!(model,SingleDofCost,[node];class=:X ,field=:tx1    ,cost=xcost)

initialstate    = initialize!(model;time=0.)
x,x′,x″         = 0.,1.,0.   
initialstate    = setdof!(initialstate,[x′];field=:tx1,nodID=[node],order=1)  # initial speed

state           = solve(SweepXA{2};  initialstate,time= t,verbose=false,catcherror=true,maxAiter=50,maxΔa=1e-10)
statewrong      = solve(SweepX{ 2};  initialstate,time= t,verbose=false,catcherror=true)

A              = getdof(state[1];class=:A,field=:ΞC,nodID=[node])
xfitted        = getdof(state;field=:tx1,nodID=[node])
xwrong         = getdof(statewrong;field=:tx1,nodID=[node])

# # Adjust the damping in the green model to match the red measurements
# using GLMakie
# fig       = Figure(size = (1000,800))
# axeX      = Axis(fig[1,1])
# scatter!(axeX,t[1:20:end],xtarget[1,1:20:end],color=:red)
# lines!(  axeX,t,xfitted[1,:],color=:black)
# lines!(  axeX,t,xwrong[1,:],color=:lightgreen)
# display(fig)

@testset "StaticX" begin
    @test  A ≈ [log10(0.4)-log10(0.6)] rtol= 1e-6
    @test  xfitted≈xtarget rtol= 1e-5
end


end