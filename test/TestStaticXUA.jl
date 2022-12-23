#module TestStaticXUA

using Test,StaticArrays,SparseArrays
using Muscade

include("SomeElements.jl")

model           = Model(:TestModel)
n1              = addnode!(model,𝕣[0,0])  # moving node
n2              = addnode!(model,𝕣[10,0]) # anchor 1 
n3              = addnode!(model,𝕣[0,10]) # anchor 2 
n4              = addnode!(model,𝕣[])     # Anod for springs
e1              = addelement!(model,Spring{2},[n1,n2,n4], EI=1)
e2              = addelement!(model,Spring{2},[n1,n3,n4], EI=1)
e3              = addelement!(model,DofLoad  ,[n1], field=:tx1      ,value=t->      t)
e4              = addelement!(model,DofHold  ,[n2], field=:tx1)
e5              = addelement!(model,DofHold  ,[n2], field=:tx2)
e6              = addelement!(model,DofHold  ,[n3], field=:tx1)
e7              = addelement!(model,DofHold  ,[n3], field=:tx2)
e8              = addelement!(model,XdofCost ,[n1], field=:tx1      ,cost=x->    1x^2)
e9              = addelement!(model,XdofCost ,[n1], field=:tx2      ,cost=x->    1x^2)
e10             = addelement!(model,AdofCost ,[n4], field=:ΞL₀      ,cost=a->  0.1a^2)
e10             = addelement!(model,AdofCost ,[n4], field=:ΞEI      ,cost=a->  0.1a^2)
state           = solve(staticX;model,time=[0.,1.],verbose=true)
describe(state[2])
state           = solve(staticXUA;model,time=[1.],verbose=true)
describe(state[1])
;
#end 

