#module TestModelDescription

using Test
using Muscade
using Muscade: DofID,EleID,NodID

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
setscale!(model;scale=(X=(tx1=1.,tx2=1.,rx3=2.),A=(Δseadrag=3.,Δskydrag=4.,ΔL=5)),Λscale=2)  # scale = (X=(tx=10,rx=1),A=(drag=3.))
@testset "model.doftyp" begin
    @test model.doftyp[1].class == :X
    @test model.doftyp[1].field == :tx1
    @test model.doftyp[1].scale ≈ 1.
    @test model.doftyp[1].dofID == [DofID(:X, 1)]

    @test model.doftyp[3].class == :A
    @test model.doftyp[3].field == :Δseadrag
    @test model.doftyp[3].scale ≈ 3.
    @test model.doftyp[3].dofID == [DofID(:A, 1)]

    @test model.doftyp[5].class == :X
    @test model.doftyp[5].field == :rx3
    @test model.doftyp[5].scale ≈ 2.
    @test model.doftyp[5].dofID == [DofID(:X, 3)]

    @test model.scaleΛ          ≈ 2.
end

@testset "model inspection" begin
    @test model.eleobj[e1] == Turbine{typeof(sea), typeof(sky)}([0.0, 0.0], -10.0, 2.0, sea, 3.0, sky) 
    @test model.nod[n1].coord ≈ [0.0, 0.0, -10.0]
    getndof(model,:X) == 3
    getndof(model,(:X,:A)) == (3,2)
    Muscade.getnele(model) == 2
    Muscade.getnele(model,1) == 1
end

@testset "get" begin
    @test get(model,e1) == (eletyp = Turbine{typeof(sea), typeof(sky)}, nods = ((inod = 1, coord = [0.0, 0.0, -10.0]), (inod = 2, coord = Float64[])), dofs = ((idof = 1, inod = 1, class = :X, field = :tx1, scale = 1.0), (idof = 2, inod = 1, class = :X, field = :tx2, scale = 1.0), (idof = 1, inod = 2, class = :A, field = :Δseadrag, scale = 3.0), (idof = 2, inod = 2, class = :A, field = :Δskydrag, scale = 4.0)))
    @test get(model,n1) == (coord = [0.0, 0.0, -10.0], eles = Any[(eleID = EleID(1, 1), eletyp = Turbine{typeof(sea), typeof(sky)}), (eleID = EleID(2, 1), eletyp = AnchorLine)], dofs = Any[(dofID = DofID(:X, 1), class = :X, field = :tx1, scale = 1.0), (dofID = DofID(:X, 2), class = :X, field = :tx2, scale = 1.0), (dofID = DofID(:X, 3), class = :X, field = :rx3, scale = 2.0)])
    @test get(model,DofID(:X, 1)) == (class = :X, field = :tx1, nodID = NodID(1), scale = 1.0)
end

#end # module
