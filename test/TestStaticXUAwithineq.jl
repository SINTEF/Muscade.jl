module TestStaticXUAwithineq

using Test
using Muscade
include("SomeElements.jl")

model           = Model(:TestModel)
n1              = addnode!(model,𝕣[ 0, 0])  # moving node
n2              = addnode!(model,𝕣[10, 0])  # anchor 1 
n3              = addnode!(model,𝕣[ 0,10])  # anchor 2 
n4              = addnode!(model,𝕣[     ])  # "A-nod" for springs
e1              = addelement!(model,Spring{2},[n1,n2,n4], EI=1)
e2              = addelement!(model,Spring{2},[n1,n3,n4], EI=1)
@once f1(t)     = t
e3              = addelement!(model,DofLoad,[n1], field=:tx1      ,value=f1)
e4              = addelement!(model,Hold   ,[n2], field=:tx1)
e5              = addelement!(model,Hold   ,[n2], field=:tx2)
e6              = addelement!(model,Hold   ,[n3], field=:tx1)
e7              = addelement!(model,Hold   ,[n3], field=:tx2)
@once f2(x,t)   = 1x^2
@once f3(a)     = 0.1a^2
e8              = addelement!(model, SingleDofCost, [n1], class=:X, field=:tx1, cost=f2)
e9              = addelement!(model, SingleDofCost, [n1], class=:X, field=:tx2, cost=f2)
e10             = addelement!(model, SingleDofCost, [n4], class=:A, field=:ΞL₀, cost=f3)
e11             = addelement!(model, SingleDofCost, [n4], class=:A, field=:ΞEI, cost=f3)
@once gap(x,u,a,t) = 1-sum(x.^2)
e12             = addelement!(model,DofConstraint , [n1], xinod=(1,1),xfield=(:tx1,:tx2), λinod=1,λclass=:U,λfield=:λcsr, gap=gap,mode=positive)
initialstate    = initialize!(model)
setdof!(initialstate,1.;class=:U,field=:λcsr)
stateX          = solve(SweepX{0};  initialstate,time=[.5,1.],maxΔx=1e-2,verbose=false)
stateXUA        = solve(StaticXUA;initialstate=stateX,γfac=0.1,verbose=false)
@testset "StaticXUA" begin
    @test stateXUA[2].X[1] ≈ [  0.16947515777519506, -0.09872146391604254, 0.0, 0.0, 0.0, 0.0,-0.9998314994395636,-0.01004064695204067,-0.00016850056043641565, 0.01004064695204067]
    @test stateXUA[2].A    ≈ [0.0042124607696232995,    0.5743380448543424]
    @test stateXUA[2].A == stateXUA[1].A
end

end 
;

