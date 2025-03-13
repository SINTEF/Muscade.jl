module TestBeamElements

using Test, Muscade, StaticArrays, LinearAlgebra
include("../examples/BeamElements.jl")

a = SA[1,0,0]
b = SA[0,1,1]
r = adjust(a,b)
R = Rodrigues(r)
u = R*a

v1      = variate{1,3}(SA[.1,.2,.3])
M1      = Rodrigues(v1)
w1,w∂v1 = value_∂{1,3}(Rodrigues⁻¹(M1))

v2      = variate{1,3}(SA[1e-7,2e-7,1e-8])
M2      = Rodrigues(v2)
w2,w∂v2 = value_∂{1,3}(Rodrigues⁻¹(M2))


@testset "rotations" begin
    @test r ≈ [0.0, -1.1107207345395913, 1.1107207345395913]
    @test u ≈ [2.220446049250313e-16, 0.7071067811865476, 0.7071067811865476]
    @test v1 ≈ w1
    @test w∂v1 ≈ I#[1 0 0;0 1 0;0 0 1]
    @test v2 ≈ w2
    @test w∂v2 ≈ I#[1 0 0;0 1 0;0 0 1]
end

###

model           = Model(:TestModel)
node1           = addnode!(model,𝕣[0,0,0])
node2           = addnode!(model,𝕣[4,3,0])
elnod           = [model.nod[n.inod] for n∈[node1,node2]]
mat             = BeamCrossSection(EA=10.,EI=3.,GJ=4.)

beam            = EulerBeam3D(elnod;mat,orient2=SVector(0.,1.,0.))

@testset "constructor" begin
    @test beam.cₘ    ≈ [2.0, 1.5, 0.0]
    @test beam.rₘ    ≈ [0.8 -0.6 0.0; 0.6 0.8 -0.0; 0.0 0.0 1.0]
    @test beam.ζgp   ≈ [-0.2886751345948129, 0.2886751345948129]
    @test beam.ζnod  ≈ [-0.5, 0.5]
    @test beam.tgₘ   ≈ [4.0, 3.0, 0.0]
    @test beam.tgₑ   ≈ [5.0, 0.0, 0.0]
    @test beam.Nε[1] ≈ [-.2, 0, 0, 0, 0, 0, .2, 0, 0, 0, 0, 0]
    @test beam.Nκ[1][2,2] ≈ -0.1385640646055102
    @test beam.Nκ[1][3,5] ≈ 0.5464101615137755
    @test beam.Nu[1][1,1] ≈ 0.7886751345948129
    @test beam.dL    ≈ [2.5, 2.5]
end


##

model           = Model(:TestModel)
node1           = addnode!(model,𝕣[0,0,0])
node2           = addnode!(model,𝕣[1,0,0])
elnod           = [model.nod[n.inod] for n∈[node1,node2]]
mat             = BeamCrossSection(EA=10.,EI=3.,GJ=4.)
beam            = EulerBeam3D(elnod;mat,orient2=SVector(0.,1.,0.))
t,SP,dbg  = 0.,(;),(status=:testing,)
U = (SVector{0,𝕣}(),)
A = SVector{0,𝕣}()

x = SVector(0.,0.,0.,0.,0.,0.,0.1,0.0,0.,0.,0.,0.); X = (x,)
R,FB=Muscade.residual(beam,   X,U,A,t,SP,dbg) 
@testset "residual tension" begin
    @test R        ≈  [-0.9999999999999998, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9999999999999998, 0.0, 0.0, 0.0, 0.0, 0.0]
    @test FB === nothing
end

x = SVector(0.,0.,0.,0.,0.,0.,0.,0.1,0.,0.,0.,0.); X = (x,)
R,FB=Muscade.residual(beam,   X,U,A,t,SP,dbg) 
@testset "residual flex" begin
    @test R        ≈  [0.305626505038752, -3.557508839178609, 0.0, 0.0, 0.0, -1.7940357448412423, -0.305626505038752, 3.557508839178609, 0.0, 0.0, 0.0, -1.7940357448412423]
    @test FB === nothing
end

x = SVector(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.1); X = (x,)
R,FB=Muscade.residual(beam,   X,U,A,t,SP,dbg) 
@testset "residual flex" begin
    @test R        ≈  [0.0, 1.8000000000002365, 0.0, 0.0, 0.0, 0.6000000000000143, 0.0, -1.8000000000002365, 0.0, 0.0, 0.0, 1.2000000000002227]
    @test FB === nothing
end

x = SVector(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.1,0.,0.); X = (x,)
R,FB=Muscade.residual(beam,   X,U,A,t,SP,dbg) 
@testset "residual torsion" begin
    @test R        ≈ [0.0, 0.0, 0.0, -0.40000000000000124, 0.0, 0.0, 0.0, 0.0, 0.0, 0.40000000000000135, 0.0, 0.0]
    @test FB === nothing
end


end