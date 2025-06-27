#module TestDirectXUA001

using Test
using Muscade
include("SomeElements.jl")

model           = Model(:TestModel)
n1              = addnode!(model,𝕣[ 0, 0])  # moving node
n2              = addnode!(model,𝕣[10, 0])  # anchor 1 
n3              = addnode!(model,𝕣[ 0,10])  # anchor 2 
n4              = addnode!(model,𝕣[     ])  # A-nod for springs
e1              = addelement!(model,Spring{2},[n1,n2,n4], EI=10)
e2              = addelement!(model,Spring{2},[n1,n3,n4], EI=10)
@once f1 f1(t)  = t
e3              = addelement!(model,DofLoad,[n1], field=:tx1      ,value=f1)
e4              = addelement!(model,Hold   ,[n2], field=:tx1)
e5              = addelement!(model,Hold   ,[n2], field=:tx2)
e6              = addelement!(model,Hold   ,[n3], field=:tx1)
e7              = addelement!(model,Hold   ,[n3], field=:tx2)
@once f2 f2(x,t)   = 1x^2
@once f3 f3(a)     = 0.1a^2
e8              = addelement!(model,SingleDofCost ,class=:X, field=:tx1,[n1]      ,cost=f2)
e9              = addelement!(model,SingleDofCost ,class=:X, field=:tx2,[n1]      ,cost=f2)
e10             = addelement!(model,SingleDofCost ,class=:A, field=:ΞL₀,[n4]      ,cost=f3)
e11             = addelement!(model,SingleDofCost ,class=:A, field=:ΞEI,[n4]      ,cost=f3)
initialstate    = initialize!(model)
stateX          = solve(SweepX{0};  initialstate,time=[0.],verbose=false)
stateXUA        = solve(DirectXUA{0,0,1};initialstate=stateX[1],time = 0:.1:1,verbose=false)
@testset "solution" begin
    @test stateXUA[2].X[1] ≈ [-0.04401449821525448, -0.07557699120101641, 0.0, 0.0, 0.0, 0.0, -0.10000328717156994, -0.0007524827404401577, 3.287171569929111e-6, 0.0007524827404401577]
    @test stateXUA[2].A    ≈ [0.003263970117473947, 0.509972206219955]
    @test stateXUA[2].A ≡ stateXUA[1].A
end
stateXUAcv           = solve(DirectXUA{0,0,1};initialstate=stateX[1],time = 0:.1:1,saveiter=true,verbose=false)
@testset "saveiter" begin
    @test stateXUAcv[7][2].X[1] ≈ stateXUA[2].X[1]
    @test stateXUAcv[7][2].A    ≈ [0.003263970117473947, 0.509972206219955]
    @test stateXUAcv[7][2].A == stateXUAcv[7][1].A
    @test stateXUAcv[6][2].X[1] ≈  [-0.04401449852921919, -0.07557700317881072, 0.0, 0.0, 0.0, 0.0, -0.10000328718956761, -0.0007524829249451943, 3.287189567597725e-6, 0.0007524829249451943]
    @test stateXUAcv[6][2].A    ≈ [0.003263970369544787, 0.5099724089221498]
    @test !(stateXUAcv[7][1].A == stateXUAcv[6][1].A)
end
#end 

