#module TestStaticXUA

using Test,StaticArrays,SparseArrays
using Muscade

include("SomeElements.jl")


model           = Model(:TestModel)
n1              = addnode!(model,𝕣[0,0,+100]) # turbine
n2              = addnode!(model,𝕣[])  # Anod for turbine 
n3              = addnode!(model,𝕣[])  # Anod for anchor
sea(t,x)        = SVector(1.,0.)*t
sky(t,x)        = SVector(0.,10.)
α(i)            = SVector(cos(i*2π/3),sin(i*2π/3))
e1              = addelement!(model,Turbine   ,[n1,n2], seadrag=1e6, sea=sea, skydrag=1e5, sky=sky)
e2              =[addelement!(model,AnchorLine,[n1,n3], Δxₘtop=vcat(5*α(i),[0.]), xₘbot=250*α(i), L=290., buoyancy=-5e3) for i∈0:2]
e3              = addelement!(model,XdofCost  ,[n1], field=:tx1      ,cost=x->0.0001x^2)
#e4              =  addelement!(model,XdofCost  ,[n1], field=:tx2      ,cost=x->0.0001x^2)
e5              = addelement!(model,AdofCost  ,[n2], field=:Δseadrag ,cost=a->10a^2)
e6              = addelement!(model,AdofCost  ,[n2], field=:Δskydrag ,cost=a->10a^2)
e7              = addelement!(model,AdofCost  ,[n3], field=:ΔL       ,cost=a->10a^2+a)
e8              = addelement!(model,AdofCost  ,[n3], field=:Δbuoyancy,cost=a->10a^2)
state           = solve(staticXUA;model,time=[0.],verbose=false)
;
#end
