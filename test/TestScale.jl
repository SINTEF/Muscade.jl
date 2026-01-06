module TestScale

using Test,StaticArrays,SparseArrays
using Muscade

include("SomeElements.jl")

model1          = Model(:TestModel)
n1              = addnode!(model1,𝕣[ 0, 0])  # moving node
n2              = addnode!(model1,𝕣[10, 0])  # anchor 1 
n3              = addnode!(model1,𝕣[ 0,10])  # anchor 2 
n4              = addnode!(model1,𝕣[     ])  # A-nod for springs
e1              = addelement!(model1,Spring{2},[n1,n2,n4], EA=10) # springs share Adofs
e2              = addelement!(model1,Spring{2},[n1,n3,n4], EA=10)
@functor with() load(t)  = 0.1*t
e3              = addelement!(model1,DofLoad,[n1], field=:tx1      ,value=load)
e3b             = addelement!(model1,DofLoad,[n1], field=:tx2      ,value=load)
e4              = addelement!(model1,Hold   ,[n2], field=:tx1)
e5              = addelement!(model1,Hold   ,[n2], field=:tx2)
e6              = addelement!(model1,Hold   ,[n3], field=:tx1)
e7              = addelement!(model1,Hold   ,[n3], field=:tx2)
@functor with() positionMeas(x,t)   = 0.5*((x-0.12t)/0.01)^2
@functor with() acost(a)     = 0.5*(a/.1)^2
e8              = addelement!(model1,SingleDofCost ,class=:X, field=:tx1,[n1]      ,cost=positionMeas)
e9              = addelement!(model1,SingleDofCost ,class=:X, field=:tx2,[n1]      ,cost=positionMeas)
e10             = addelement!(model1,SingleAcost   ,          field=:ΞL₀,[n4]      ,cost=acost)
e11             = addelement!(model1,SingleAcost   ,          field=:ΞEI,[n4]      ,cost=acost)

model2          = deepcopy(model1)

initialstate1    = initialize!(model1)
stateX1          = solve(SweepX{0};  initialstate=initialstate1,time=[0.,1.],verbose=false)
stateXUA1        = solve(DirectXUA{0,0,1};initialstate=[stateX1[1]],time = [0:.1:1],maxΔλ=.5,maxΔa=1e-4,maxΔx=1e-4,verbose=false)

setscale!(model2;scale=(X=(tx1=1.,tx2=10.,rx3=2.),A=(Δseadrag=3.,Δskydrag=4.,ΔL=5)),Λscale=2)  
initialstate2    = initialize!(model2)
stateX2          = solve(SweepX{0};  initialstate=initialstate2,time=[0.],verbose=false)
stateXUA2        = solve(DirectXUA{0,0,1};initialstate=[stateX2[1]],time = [0:.1:1],maxΔλ=.5,maxΔa=1e-4,maxΔx=1e-4,verbose=false)

step = 1
iexp = 1
@testset "ScaleSweepX" begin
    @test  stateX1[step].X[1] ≈ stateX2[step].X[1]
end
@testset "ScaleDirectXUAstepwise" begin
    @test  stateXUA1[iexp][step].X[1] ≈ stateXUA2[iexp][step].X[1]
    @test  stateXUA1[iexp][step].A    ≈ stateXUA2[iexp][step].A
end

scale=Muscade.study_scale(stateX1[end],verbose=false) # scaling suggestions from looking at a state
@testset "study_scale" begin
    @test scale.Λ.tx1 ≈ 81000.0
    @test scale.Λ.tx2 ≈ 81000.0
    @test scale.X.λtx2 ≈ 1.2e-5
    @test scale.U == (;)
    @test scale.A.ΞL₀ ≈  6.3e-7
end    
end
