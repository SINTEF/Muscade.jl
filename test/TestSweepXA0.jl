#module TestSweepXA2
using Test
using Muscade


include("SomeElements.jl")
include("../examples/DryFriction.jl")


K,C,M           = .3,1.,1.

model           = Model(:TestModel)
node            = addnode!(model,𝕣[])
osc             = addelement!(model,AdjustableSdofOscillator,[node]; K,C,M)

@functor with(     ) acost(a,σ)=(a/σ)^2
@functor with(σ=.1) xcost(x,t)=(x/σ)^2
@functor with()  xload(t)= 1.#cos(t)
cK              = addelement!(model,SingleAcost  ,[node];          field=:ΞK,costargs=(.2,),cost =acost)
cX              = addelement!(model,SingleDofCost,[node];class=:X ,field=:tx1,              cost =xcost)
cF              = addelement!(model,DofLoad      ,[node];          field=:tx1,              value=xload)

initialstate    = initialize!(model;time=0.)
x,x′,x″         = 0.,1.,0.   
initialstate    = setdof!(initialstate,[x′];field=:tx1,nodID=[node],order=1)  # initial speed

t               = 2.:2.
state           = solve(SweepXA{0};  initialstate,time= t,verbose=true,catcherror=true,maxAiter=100)


;
