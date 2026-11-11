#module TestSweepXA2
using Test
using Muscade


include("SomeElements.jl")
include("../examples/DryFriction.jl")

# Δt              = 0.01
Δt              = 1.
t               = Δt:Δt:10

### Target

K,C,M           = 1.,.4,.3

model           = Model(:Target)
node            = addnode!(model,𝕣[])
osc             = addelement!(model,AdjustableSdofOscillator,[node]; K,C,M)
initialstate    = initialize!(model;time=0.)

t1              = t
x,x′,x″         = 0.,1.,0.   
initialstate1   = setdof!(initialstate,[x′];field=:tx1,nodID=[node],order=1)  # initial speed
state1          = solve(SweepX{ 2};  initialstate=initialstate1,time= t1,verbose=false,catcherror=true)

t2              = t.+20.
x,x′,x″         = 1.,0.,0.   
initialstate2   = setdof!(initialstate,[x];field=:tx1,nodID=[node],order=0)  # initial speed
state2          = solve(SweepX{ 2};  initialstate=initialstate2,time= t2,verbose=false,catcherror=true)

motion1         = getdof(state1;field=:tx1,nodID=[node])[1,:]
motion2         = getdof(state2;field=:tx1,nodID=[node])[1,:]
target          = interpolator([t1;t2],[motion1;motion2])


### Optimize

K,C,M           = 1.,.6,.3

model           = Model(:Optim)
node            = addnode!(model,𝕣[])
osc             = addelement!(model,AdjustableSdofOscillator,[node]; K,C,M)
@functor with(σ=10.        )  acost(a  )=( a           /σ)^2
@functor with(σ=.0001,target) xcost(x,t)=((x-target(t))/σ)^2
cC              = addelement!(model,SingleAcost  ,[node];field=:ΞC,               cost=acost)
cX              = addelement!(model,SingleDofCost,[node];class=:X ,field=:tx1    ,cost=xcost)


initialstate    = initialize!(model;time=0.)

x,x′,x″         = 0.,1.,0.   
initialstate1   = setdof!(initialstate,[x′];field=:tx1,nodID=[node],order=1)  # initial speed
sweep1          = Sweep{2}(isweep=1,time=t1,dependson=-1)

x,x′,x″         = 1.,0.,0.   
initialstate2   = setdof!(initialstate,[x];field=:tx1,nodID=[node],order=0)  # initial position
sweep2          = Sweep{2}(isweep=2,time=t2,dependson=-2)
@assert ~(initialstate1 ≡ initialstate2)

state           = solve(SweepXA; initialstates=[initialstate1,initialstate2],sweeps=[sweep1,sweep2],
                        verbose=true,catcherror=true,maxAiter=2,maxΔa=1e-10,multithread=false)
statewrong1     = solve(SweepX{ 2};  initialstate=initialstate1,time= t1,verbose=false,catcherror=true)
statewrong2     = solve(SweepX{ 2};  initialstate=initialstate2,time= t2,verbose=false,catcherror=true)

A               = getdof(state[1][1];class=:A,field=:ΞC,nodID=[node])
xfitted1        = getdof(state[1];field=:tx1,nodID=[node])[1,:]
xfitted2        = getdof(state[2];field=:tx1,nodID=[node])[1,:]
xwrong1         = getdof(statewrong1;field=:tx1,nodID=[node])[1,:]
xwrong2         = getdof(statewrong2;field=:tx1,nodID=[node])[1,:]

# # Adjust the damping in the green model to match the red measurements
# using GLMakie
# fig       = Figure(size = (1000,800))
# axeX      = Axis(fig[1,1])
# scatter!(axeX,t[1:20:end],xtarget[1,1:20:end],color=:red)
# lines!(  axeX,t,xfitted[1,:],color=:black)
# lines!(  axeX,t,xwrong[1,:],color=:lightgreen)
# display(fig)

# @testset "StaticX" begin
#     @test  A ≈ [log10(0.4)-log10(0.6)] rtol= 1e-6
#     @test  xfitted≈xtarget rtol= 1e-5
# end


#end