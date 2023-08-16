module TestStaticXUAstepwise

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
e8              = addelement!(model,SingleDofCost ,class=:X, field=:tx1,[n1]      ,cost=f2)
e9              = addelement!(model,SingleDofCost ,class=:X, field=:tx2,[n1]      ,cost=f2)
e10             = addelement!(model,SingleDofCost ,class=:A, field=:ΞL₀,[n4]      ,cost=f3)
e11             = addelement!(model,SingleDofCost ,class=:A, field=:ΞEI,[n4]      ,cost=f3)
initialstate    = initialize!(model)
@testset "StaticX" begin
    stateX           = solve(StaticX;initialstate,time=[0.,1.],verbose=false)
    @test stateX[2].X[1] ≈ [1.0008305423582147, 0.056562064402881376, 0.0, 0.0, 0.0, 0.0, 
    -1.0006330261310146, 0.006289232571302629, 
    0.0006330261310144907, -0.006289232571302629]
end
stateX             = solve(StaticX;  initialstate,time=[.5,1.],verbose=false)
stateXUA           = solve(StaticXUAstepwise;initialstate=stateX,maxYiter= 50,verbose=false)
@testset "StaticXUAstepwise" begin
    @test stateXUA[2].X[1] ≈ [  0.16947517267111387,    -0.09872147216175686,     0.0,     0.0,     0.0,     0.0,    -0.9998314994105624,    -0.01004064780561606,    -0.00016850058943765545,     0.01004064780561606]
    @test stateXUA[2].A    ≈ [0.004212461115295247,    0.5743380076037062]
    @test stateXUA[2].A == stateXUA[1].A
end
stateXUA           = solve(StaticXUAstepwise;initialstate=stateX,saveiter=true,verbose=false)
@testset "StaticXUAstepwise-saveiter" begin
    @test stateXUA[6][2].X[1] ≈ [  0.16947517267111387,    -0.09872147216175686,     0.0,     0.0,     0.0,     0.0,    -0.9998314994105624,    -0.01004064780561606,    -0.00016850058943765545,     0.01004064780561606]
    @test stateXUA[6][2].A    ≈ [0.004212461115295247,    0.5743380076037062]
    @test stateXUA[6][2].A == stateXUA[6][1].A
    @test stateXUA[5][2].X[1] ≈ [0.16947513467261088, -0.09872142687522917, 0.0, 0.0, 0.0, 0.0, -0.9998315045843422, -0.010040642705804102, -0.00016849541565779142, 0.010040642705804102]
    @test stateXUA[5][2].A    ≈ [0.004212461757011719, 0.5743379495209457]
    @test !(stateXUA[5][1].A ≡ stateXUA[6][1].A)
end

end 

