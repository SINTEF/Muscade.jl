module TestModelDescription

using Test
using Muscade

include("SomeElements.jl")

model           = Model(:TestModel)
n1              = addnode!(model,𝕣[0,0,-10])
n2              = addnode!(model,𝕣[])
n3              = addnode!(model,𝕣[])
@testset "Nodes" begin
    @test n1 == 1
    @test n2 == 2
end

# the arguments to element constructors are irrelevant - elements are not tested here
sea(t,x)        = SVector(1.,0.)
sky(t,x)        = SVector(0.,1.)
e1              = addelement!(model,Turbine   ,[n1,n2], seadrag=2., sea=sea, skydrag=3., sky=sky)
e2              = addelement!(model,AnchorLine,[n1,n3], Δxₘtop=SVector(0,2.,0), xₘbot=SVector(94.,0.), L=170., buoyancy=-1.)
@testset "Receipts" begin
    @test n1 == 1
    @test n2 == 2
    @test e1 == 1
    @test e2 == 2
end
@testset "model.nod" begin
    @test model.nod[1].ID    == 1
    @test model.nod[1].coord ≈  [0.0, 0.0, -10.0]
    @test model.nod[1].dofID == [1, 2, 5]
    @test model.nod[1].eleID == [1, 2]

    @test model.nod[2].ID    == 2
    @test model.nod[2].coord == 𝕣[]
    @test model.nod[2].dofID == [3,4]
    @test model.nod[2].eleID == [1]

    @test model.nod[3].ID    == 3
    @test model.nod[3].coord == 𝕣[]
    @test model.nod[3].dofID == [6,7]
    @test model.nod[3].eleID == [2]
end 
@testset "model.ele" begin
    @test model.ele[1].ID       == 1
    @test model.ele[1].eletypID == 1
    @test model.ele[1].iele     == 1
    @test model.ele[1].nodID    == [1,2]
    @test model.ele[1].dofID    == [1,2,3,4]

    @test model.ele[2].ID       == 2
    @test model.ele[2].eletypID == 2
    @test model.ele[2].iele     == 1
    @test model.ele[2].nodID    == [1,3]
    @test model.ele[2].dofID    == [1,2,5,6,7]
end
@testset "model.dof" begin
    @test model.dof[1].ID       == 1
    @test model.dof[1].nodID    == 1
    @test model.dof[1].doftypID == 1
    @test model.dof[1].eleID    == [1,2]

    @test model.dof[3].ID       == 3
    @test model.dof[3].nodID    == 2
    @test model.dof[3].doftypID == 3
    @test model.dof[3].eleID    == [1]

    @test model.dof[5].ID       == 5
    @test model.dof[5].nodID    == 1
    @test model.dof[5].doftypID == 5
    @test model.dof[5].eleID    == [2]

    @test model.dof[6].ID       == 6
    @test model.dof[6].nodID    == 3
    @test model.dof[6].doftypID == 6
    @test model.dof[6].eleID    == [2]
end
@testset "model.eletyp" begin
    @test model.eletyp[1].ID    == 1
    @test model.eletyp[1].type  == Turbine{typeof(sea), typeof(sky)}

    @test model.eletyp[2].ID    == 2
    @test model.eletyp[2].type  == AnchorLine
end
@testset "model.doftyp" begin
    @test model.doftyp[1].ID    == 1
    @test model.doftyp[1].class == :X
    @test model.doftyp[1].field == :tx1
    @test model.doftyp[1].scale ≈ 1.

    @test model.doftyp[3].ID    == 3
    @test model.doftyp[3].class == :A
    @test model.doftyp[3].field == :Δseadrag
    @test model.doftyp[3].scale ≈ 1.

    @test model.doftyp[5].ID    == 5
    @test model.doftyp[5].class == :X
    @test model.doftyp[5].field == :rx3
    @test model.doftyp[5].scale ≈ 1.
end

end # module
