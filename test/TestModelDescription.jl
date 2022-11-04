module TestModelDescription

using Test
using Muscade

include("SomeElements.jl")

model           = Model(:TestModel)
n1              = addnode!(model,𝕣[0,0,-10])
n2              = addnode!(model,𝕣[])
n3              = addnode!(model,𝕣[])
@testset "Nodes" begin
    @test n1 == Muscade.NodID(1)
    @test n2 == Muscade.NodID(2)
end

# the arguments to element constructors are irrelevant - elements are not tested here
sea(t,x)        = SVector(1.,0.)
sky(t,x)        = SVector(0.,1.)
e1              = addelement!(model,Turbine   ,[n1,n2], seadrag=2., sea=sea, skydrag=3., sky=sky)
e2              = addelement!(model,AnchorLine,[n1,n3], Δxₘtop=SVector(0,2.,0), xₘbot=SVector(94.,0.), L=170., buoyancy=-1.)
@testset "Receipts" begin
    @test n1 == NodID(1)
    @test n2 == NodID(2)
    @test e1 == EleID(1,1)
    @test e2 == EleID(2,1)
end
@testset "model.nod" begin
    @test model.nod[1].ID    == NodID(1)
    @test model.nod[1].coord ≈  [0.0, 0.0, -10.0]
    @test model.nod[1].dofID == [DofID(:X,1), DofID(:X,2), DofID(:X,3)]
    @test model.nod[1].eleID == [EleID(1,1), EleID(2,1)]

    @test model.nod[2].ID    == NodID(2)
    @test model.nod[2].coord == 𝕣[]
    @test model.nod[2].dofID == [DofID(:A,1), DofID(:A,2)]#[3,4]
    @test model.nod[2].eleID == [EleID(1,1)]

    @test model.nod[3].ID    == NodID(3)
    @test model.nod[3].coord == 𝕣[]
    @test model.nod[3].dofID == [DofID(:A,3), DofID(:A,4)]
    @test model.nod[3].eleID == [EleID(2,1)]
end 
@testset "model.ele" begin
    @test model.ele[1][1].ID       == EleID(1,1)
    @test model.ele[1][1].ieletyp  == 1
    @test model.ele[1][1].iele     == 1
    @test model.ele[1][1].nodID    == [NodID(1),NodID(2)]
    @test model.ele[1][1].dofID    == [DofID(:X, 1), DofID(:X, 2), DofID(:A, 1), DofID(:A, 2)]

    @test model.ele[2][1].ID       == EleID(2,1)
    @test model.ele[2][1].ieletyp  == 2
    @test model.ele[2][1].iele     == 1
    @test model.ele[2][1].nodID    == [NodID(1),NodID(3)]
    @test model.ele[2][1].dofID    == [DofID(:X, 1), DofID(:X, 2), DofID(:X, 3), DofID(:A, 3), DofID(:A, 4)]
end
@testset "model.dof" begin
    @test model.dof[DofID(:X,1)].ID       == DofID(:X,1)
    @test model.dof[DofID(:X,1)].nodID    == NodID(1)
    @test model.dof[DofID(:X,1)].idoftyp  == 1
    @test model.dof[DofID(:X,1)].eleID    == [EleID(1,1),EleID(2,1)]

    @test model.dof[DofID(:A,1)].ID       == DofID(:A,1)
    @test model.dof[DofID(:A,1)].nodID    == NodID(2)
    @test model.dof[DofID(:A,1)].idoftyp  == 3
    @test model.dof[DofID(:A,1)].eleID    == [EleID(1,1)]

    @test model.dof[DofID(:X,3)].ID       == DofID(:X,3)
    @test model.dof[DofID(:X,3)].nodID    == NodID(1)
    @test model.dof[DofID(:X,3)].idoftyp  == 5
    @test model.dof[DofID(:X,3)].eleID    == [EleID(2,1)]

    @test model.dof[DofID(:A,3)].ID       == DofID(:A,3)
    @test model.dof[DofID(:A,3)].nodID    == NodID(3)
    @test model.dof[DofID(:A,3)].idoftyp  == 6
    @test model.dof[DofID(:A,3)].eleID    == [EleID(2,1)]
end
@testset "model.doftyp" begin
    @test model.doftyp[1].class == :X
    @test model.doftyp[1].field == :tx1
    @test model.doftyp[1].scale ≈ 1.

    @test model.doftyp[3].class == :A
    @test model.doftyp[3].field == :Δseadrag
    @test model.doftyp[3].scale ≈ 1.

    @test model.doftyp[5].class == :X
    @test model.doftyp[5].field == :rx3
    @test model.doftyp[5].scale ≈ 1.
end

end # module
