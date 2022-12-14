module TestAssemble

using Test,StaticArrays,SparseArrays
using Muscade
using Muscade.ElTest

include("SomeElements.jl")

### Turbine
sea(t,x) = SVector(1.,0.)
sky(t,x) = SVector(0.,1.)
turbine  = Turbine(SVector(0.,0.),-10., 2.,sea, 3.,sky)
δX       = @SVector [1.,1.]
X        = @SVector [1.,2.]
U        = @SVector 𝕣[]
A        = @SVector [0.,0.]  # [Δseadrag,Δskydrag]


L,Lδx,Lx,Lu,La   = Muscade.gradient(turbine,δX,[X],[U],A, 0.,0.,())

@testset "Turbine gradient" begin
    @test Lδx           ≈ [-2, -3]
    @test Lx            ≈ [0, 0]
    @test length(Lu)    == 0
    @test La            ≈ [-2, -3]
end

Lδx,Lx,Lu,La   = test_static_element(turbine;δX,X,U,A,verbose=false)

@testset "test_static_element" begin
    @test Lδx           ≈ [-2, -3]
    @test Lx            ≈ [0, 0]
    @test length(Lu)    == 0
    @test La            ≈ [-2, -3]
end

# ###  AnchorLine

anchorline      = AnchorLine(SVector(0.,0.,100.), SVector(0,2.,0), SVector(94.,0.), 170., -1.)

δX       = @SVector [1.,1.,1.]
X        = @SVector [0.,0.,0.]
U        = @SVector 𝕣[]
A        = @SVector [0.,0.]  # [Δseadrag,Δskydrag]
L,Lδx,Lx,Lu,La   = Muscade.gradient(anchorline,δX,[X],[U],A, 0.,0.,())
@testset "anchorline1" begin
    @test Lδx           ≈ [-12.25628901693551, 0.2607721067433087, 24.51257803387102]
    @test Lx            ≈ [-0.91509745608786, 0.14708204066349, 1.3086506986891027]
    @test length(Lu)    == 0
    @test La            ≈ [-156.06324599170992, 12.517061123678818]
end

model           = Model(:TestModel)
n1              = addnode!(model,𝕣[0,0,+100])
n2              = addnode!(model,𝕣[])
n3              = addnode!(model,𝕣[])
sea(t,x)        = SVector(1.,0.)
sky(t,x)        = SVector(0.,1.)
e1              = addelement!(model,Turbine   ,[n1,n2], seadrag=2., sea=sea, skydrag=3., sky=sky)
e2              = addelement!(model,AnchorLine,[n1,n3], Δxₘtop=SVector(5.,0.,0), xₘbot=SVector(150.,0.), L=180., buoyancy=-1e3)
setscale!(model;scale=(X=(tx1=1.,tx2=1.,rx3=2.),A=(Δseadrag=3.,Δskydrag=4.,ΔL=5)),Λscale=2)  # scale = (X=(tx=10,rx=1),A=(drag=3.))
dis = Muscade.Disassembler(model)


@testset "Disassembler" begin
    @test  dis.dis[1].index[1].X == [1,2]
    @test  dis.dis[1].index[1].U == []
    @test  dis.dis[1].index[1].A == [1,2]
    @test  dis.dis[2].index[1].X == [1,2,3]
    @test  dis.dis[2].index[1].U == []
    @test  dis.dis[2].index[1].A == [3,4]
    @test  dis.scaleX ≈  [1,1,2]
    @test  dis.scaleU ≈  𝕫[]
    @test  dis.scaleA ≈  [3,4,5,1]
end

dofgr       = Muscade.allXdofs(model,dis)
Λ,X,U,A     = Muscade.indexedstate(dofgr)
nΛ,nX,nU,nA = Muscade.gradientstructure(dofgr,dis.dis[1]) # number of dofs of each class in the gradient returned by an element
iΛ,iX,iU,iA = Muscade.gradientpartition(nΛ,nX,nU,nA)  # indices into said gradient
@testset "dofgr" begin
    @test dofgr.nX == 3
    @test dofgr.nU == 0
    @test dofgr.nA == 4
    @test dofgr.iΛ == Int64[]
    @test dofgr.iX == 1:3
    @test dofgr.iU == Int64[]
    @test dofgr.iA == Int64[]
    @test dofgr.jΛ == 1:0
    @test dofgr.jX == 1:3
    @test dofgr.jU == 4:3
    @test dofgr.jA == 4:3
    @test dofgr.scaleΛ == Float64[]
    @test dofgr.scaleX ≈  [1.0, 1.0, 2.0]
    @test dofgr.scaleU == Float64[]
    @test dofgr.scaleA == Float64[]
end
@testset "state" begin
    @test Λ == [0, 0 ,0]
    @test X == [1, 2, 3]
    @test U == Int64[]
    @test A == [0, 0, 0, 0]
    @test nΛ == 0
    @test nX == 2
    @test nU == 0
    @test nA == 0
    @test iΛ == 1:0
    @test iX == 1:2
    @test iU == 3:2
    @test iA == 3:2
end

neletyp     = 2
asmvec      = Vector{𝕫2}(undef,neletyp)  
Lλ          = Muscade.asmvec!(asmvec,dofgr,dis) 
@testset "asmvec" begin
    @test typeof(Lλ)== Vector{Float64}
    @test length(Lλ)== 3
    @test asmvec[1] == [1; 2;;]
    @test asmvec[2] == [1; 2; 3;;]
end
out,asm,dofgr = Muscade.prepare(Muscade.OUTstaticX,model,dis)
Muscade.zero!(out)
@testset "prepare" begin
    @test  out.Lλ ≈ [0,0,0]
    @test  out.Lλx ≈ sparse([1,2,3,2,3], [1,2,2,3,3], [0,0,0,0,0], 3, 3)
    @test  asm[1,1] == [1; 2;;]
    @test  asm[1,2] == [1; 2; 3;;]
    @test  asm[2,1] == [1; 2; 4; 5;;]
    @test  asm[2,2] == [1; 2; 3; 4; 5; 6; 7; 8; 9;;]
end
@testset "dofgr again" begin
    @test dofgr.nX == 3
    @test dofgr.nU == 0
    @test dofgr.nA == 4
    @test dofgr.iΛ == Int64[]
    @test dofgr.iX == 1:3
    @test dofgr.iU == Int64[]
    @test dofgr.iA == Int64[]
    @test dofgr.jΛ == 1:0
    @test dofgr.jX == 1:3
    @test dofgr.jU == 4:3
    @test dofgr.jA == 4:3
    @test dofgr.scaleΛ == Float64[]
    @test dofgr.scaleX ≈  [1.0, 1.0, 2.0]
    @test dofgr.scaleU == Float64[]
    @test dofgr.scaleA == Float64[]
end

state = Muscade.State(model,dis)
Muscade.assemble!(out,asm,dis,model,state, 0.,())

@testset "assemble" begin
    @test  out.Lλ ≈ [-304261.42399716884, -6.0, 0.0]
    @test  out.Lλx ≈ sparse([1,2,3,2,3], [1,2,2,3,3], [20646.13919595113, 2098.3270620494404, 20983.270620494404, 20983.2706204944, 6.294981186148321e6], 3, 3)
end

end
