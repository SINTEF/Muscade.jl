module TestStaticXUA

using Test
using Muscade


model           = Model(:TestModel)
n1              = addnode!(model,𝕣[ 0, 0])  # moving node
n2              = addnode!(model,𝕣[10, 0])  # anchor 1 
n3              = addnode!(model,𝕣[ 0,10])  # anchor 2 
n4              = addnode!(model,𝕣[     ])  # A-nod for springs
e1              = addelement!(model,Spring{2},[n1,n2,n4], EI=1)
e2              = addelement!(model,Spring{2},[n1,n3,n4], EI=1)
@once f1 f1(t)  = t
e3              = addelement!(model,DofLoad  ,[n1], field=:tx1      ,value=f1)
e4              = addelement!(model,Hold  ,[n2], field=:tx1)
e5              = addelement!(model,Hold  ,[n2], field=:tx2)
e6              = addelement!(model,Hold  ,[n3], field=:tx1)
e7              = addelement!(model,Hold  ,[n3], field=:tx2)
@once f2 f2(x)  = 1x^2
@once f3 f3(a)  = 0.1a^2
e8              = addelement!(model,XdofCost ,[n1], field=:tx1      ,cost=f2)
e9              = addelement!(model,XdofCost ,[n1], field=:tx2      ,cost=f2)
e10             = addelement!(model,AdofCost ,[n4], field=:ΞL₀      ,cost=f3)
e10             = addelement!(model,AdofCost ,[n4], field=:ΞEI      ,cost=f3)
@testset "StaticX" begin
    stateX           = solve(staticX;model,time=[0.,1.],verbose=false)
    @test stateX[2].X[1] ≈ [ 1.000830542358214,    0.056562064402879385,    0.0,    0.0,    0.0,    0.0,   -1.0006330261310143,    0.006289232571302405,    0.0006330261310144671,   -0.006289232571302405]
end
@testset "StaticXUA" begin
    stateX             = solve(staticX;  model,time=[.5,1.],verbose=false)
    stateXUA           = solve(staticXUA;model,initial=stateX,maxYiter= 50,verbose=false)
    @test stateXUA[2].X[1] ≈ [  0.16947517267111387,    -0.09872147216175686,     0.0,     0.0,     0.0,     0.0,    -0.9998314994105624,    -0.01004064780561606,    -0.00016850058943765545,     0.01004064780561606]
    @test stateXUA[2].A    ≈ [0.004212461115295247,    0.5743380076037062]
    @test stateXUA[2].A == stateXUA[1].A
end
end 

