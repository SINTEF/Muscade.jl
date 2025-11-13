#module TestSweepXA0
using Test
using Muscade


include("SomeElements.jl")
#include("../examples/DryFriction.jl")


K,C,M           = 1.,.3,0.

model           = Model(:TestModel)
node            = addnode!(model,𝕣[])
osc             = addelement!(model,AdjustableSdofOscillator,[node]; K,C,M)

@functor with(    )  acost(a,σ)=(a/σ)^2
@functor with(σ=.1)  xcost(x,t)=(x/σ)^2
@functor with(    )  xload(  t)= .5cos(3t)
# cK              = addelement!(model,SingleAcost  ,[node];          field=:ΞK,costargs=(.2,),cost =acost)
#cX              = addelement!(model,SingleDofCost,[node];class=:X ,field=:tx1,              cost =xcost)
#cF               = addelement!(model,DofLoad      ,[node];          field=:tx1,              value=xload)

initialstate    = initialize!(model;time=0.)
x,x′,x″         = 0.3,1.,0.   
initialstate    = setdof!(initialstate,[x′];field=:tx1,nodID=[node],order=1)  # initial speed

t               = .1:.1:4π
state           = solve(SweepX{2};  initialstate,time=t)
#state           = solve(SweepXA{1};  initialstate,time= t,verbose=true,catcherror=true,maxAiter=100)

x = [state[i].X[1][1] for i = 1:length(t)]
#λ = [state[i].Λ[1][1] for i = 1:length(t)]


using GLMakie

fig      = Figure(size = (600,450))
display(fig) # open interactive window (gets closed down by "save")
axe      = Axis(fig[1,1],title="Test",xlabel="t")
ox    = lines!(  axe,t,x    ,color = :black, linewidth = 1.)
#oλ    = lines!(  axe,t,λ/100,color = :red, linewidth = 1.)
#@show state[1].A[1]
;#end
