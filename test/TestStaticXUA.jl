module TestStaticXUA

using Test
using Muscade

include("SomeElements.jl")

model           = Model(:TestModel)
n1              = addnode!(model,𝕣[ 0, 0])  # moving node
n2              = addnode!(model,𝕣[10, 0])  # anchor 1 
n3              = addnode!(model,𝕣[ 0,10])  # anchor 2 
n4              = addnode!(model,𝕣[     ])  # A-nod for springs
e1              = addelement!(model,Spring{2},[n1,n2,n4], EI=1)
e2              = addelement!(model,Spring{2},[n1,n3,n4], EI=1)
@once f1(t)     = t
e3              = addelement!(model,DofLoad  ,[n1], field=:tx1      ,value=f1)
e4              = addelement!(model,Hold  ,[n2], field=:tx1)
e5              = addelement!(model,Hold  ,[n2], field=:tx2)
e6              = addelement!(model,Hold  ,[n3], field=:tx1)
e7              = addelement!(model,Hold  ,[n3], field=:tx2)
@once f2(x,t)   = 1x^2
@once f3(a)     = 0.1a^2
e8              = addelement!(model,DofCost ,class=:X, field=:tx1,[n1]      ,cost=f2)
e9              = addelement!(model,DofCost ,class=:X, field=:tx2,[n1]      ,cost=f2)
e10             = addelement!(model,DofCost ,class=:A, field=:ΞL₀,[n4]      ,cost=f3)
e11             = addelement!(model,DofCost ,class=:A, field=:ΞEI,[n4]      ,cost=f3)
initialstate    = initialize!(model)
@testset "StaticX" begin
    stateX           = solve(StaticX;initialstate,time=[0.,1.],verbose=false)
    @test stateX[2].X[1] ≈ [ 1.000830542358214,    0.056562064402879385,    0.0,    0.0,    0.0,    0.0,   -1.0006330261310143,    0.006289232571302405,    0.0006330261310144671,   -0.006289232571302405]
end
@testset "StaticXUA" begin
    stateX             = solve(StaticX;  initialstate,time=[.5,1.],verbose=false)
    stateXUA           = solve(StaticXUA;initialstate=stateX,maxYiter= 50,verbose=false)
    @test stateXUA[2].X[1] ≈ [  0.16947517267111387,    -0.09872147216175686,     0.0,     0.0,     0.0,     0.0,    -0.9998314994105624,    -0.01004064780561606,    -0.00016850058943765545,     0.01004064780561606]
    @test stateXUA[2].A    ≈ [0.004212461115295247,    0.5743380076037062]
    @test stateXUA[2].A == stateXUA[1].A
end
end 

