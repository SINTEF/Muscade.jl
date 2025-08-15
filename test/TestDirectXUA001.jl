module TestDirectXUA001

using Test
using Muscade
include("SomeElements.jl")

model           = Model(:TestModel)
n1              = addnode!(model,𝕣[ 0, 0])  # moving node
n2              = addnode!(model,𝕣[10, 0])  # anchor 1 
n3              = addnode!(model,𝕣[ 0,10])  # anchor 2 
n4              = addnode!(model,𝕣[     ])  # A-nod for springs
e1              = addelement!(model,Spring{2},[n1,n2,n4], EA=10) # springs share Adofs
e2              = addelement!(model,Spring{2},[n1,n3,n4], EA=10)
@once load load(t)  = 0.1*t
e3              = addelement!(model,DofLoad,[n1], field=:tx1      ,value=load)
e3b             = addelement!(model,DofLoad,[n1], field=:tx2      ,value=load)
e4              = addelement!(model,Hold   ,[n2], field=:tx1)
e5              = addelement!(model,Hold   ,[n2], field=:tx2)
e6              = addelement!(model,Hold   ,[n3], field=:tx1)
e7              = addelement!(model,Hold   ,[n3], field=:tx2)
@once positionMeas positionMeas(x,t)   = 0.5*((x-0.12t)/0.01)^2
@once acost acost(a)     = 0.5*(a/.1)^2
e8              = addelement!(model,SingleDofCost ,class=:X, field=:tx1,[n1]      ,cost=positionMeas)
e9              = addelement!(model,SingleDofCost ,class=:X, field=:tx2,[n1]      ,cost=positionMeas)
e10             = addelement!(model,SingleDofCost ,class=:A, field=:ΞL₀,[n4]      ,cost=acost)
e11             = addelement!(model,SingleDofCost ,class=:A, field=:ΞEI,[n4]      ,cost=acost)
initialstate    = initialize!(model)
stateX          = solve(SweepX{0};  initialstate,time=[0.],verbose=false)
stateXUA        = solve(DirectXUA{0,0,1};initialstate=[stateX[1]],time = [0:.1:1],maxΔλ=.5,maxΔa=1e-4,maxΔx=1e-4,verbose=false,maxiter=50)
iexp = 1
@testset "solution" begin
    @test stateXUA[iexp][2].X[1]' ≈ [0.0153749  0.0153749  0.0  0.0  0.0  0.0  -0.0100154  1.54223e-5  1.54223e-5  -0.0100154] rtol=1e-4
    @test stateXUA[iexp][2].A'    ≈ [ -0.000189233  -0.0411784] rtol=1e-4
    @test stateXUA[iexp][2].A ≡ stateXUA[1][1].A
end
stateXUAcv           = solve(DirectXUA{0,0,1};initialstate=[stateX[1]],time = [0:.1:1],saveiter=true,maxΔλ=.5,maxΔa=1e-4,maxΔx=1e-4,verbose=false)
jiter = findlastassigned(stateXUAcv)
@testset "saveiter" begin
    @test stateXUAcv[jiter][iexp][2].X[1] ≈ stateXUA[iexp][2].X[1]
    @test stateXUAcv[jiter][iexp][2].A    ≈ [ -0.0001892328954715017, -0.04117844040354315]
    @test stateXUAcv[1][iexp][2].A === stateXUAcv[1][iexp][1].A
    @test stateXUAcv[jiter][iexp][2].X[1]' ≈  [0.0153749  0.0153749  0.0  0.0  0.0  0.0  -0.0100154  1.54223e-5  1.54223e-5  -0.0100154] rtol=1e-4
    @test !(stateXUAcv[jiter][iexp][1].A == stateXUAcv[1][iexp][1].A)
end
end 

