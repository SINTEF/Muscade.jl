module TestEulerBeam3D

using Test, Muscade, StaticArrays

a = SA[1,0,0]
b = SA[0,1,1]
r = Elements.adjust(a,b)
R = Elements.Rodrigues(r)
u = R*a

v1      = variate{1,3}(SA[.1,.2,.3])
M1      = Elements.Rodrigues(v1)
w1,w∂v1 = value_∂{1,3}(Elements.Rodrigues⁻¹(M1))

v2      = variate{1,3}(SA[1e-7,2e-7,1e-8])
M2      = Elements.Rodrigues(v2)
w2,w∂v2 = value_∂{1,3}(Elements.Rodrigues⁻¹(M2))

v3      = variate{1,3}(SA[1e-7,2e-7,1e-8])
M3      = Elements.Rodrigues(v3)
w3,w∂v3 = value_∂{1,3}(Elements.Rodrigues⁻¹(M3))


@testset "rotations" begin
    @test r ≈ [0.0, -1.1107207345395913, 1.1107207345395913]
    @test u ≈ [2.220446049250313e-16, 0.7071067811865476, 0.7071067811865476]
    @test v1 ≈ w1
    @test w∂v1 ≈ I#[1 0 0;0 1 0;0 0 1]
    @test v2 ≈ w2
    @test w∂v2 ≈ I#[1 0 0;0 1 0;0 0 1]
    @test v3 ≈ w3
    @test w∂v3 ≈ I#[1 0 0;0 1 0;0 0 1]
end



###

model           = Model(:TestModel)
node1           = addnode!(model,𝕣[0,0,0])
node2           = addnode!(model,𝕣[4,3,0])
elnod           = [model.nod[n.inod] for n∈[node1,node2]]
mat             = Elements.BeamCrossSection(EA=10.,EI=3.,GJ=4.)
beam            = Elements.EulerBeam3D(elnod;mat,orient2=SVector(0.,1.,0.))

@testset "constructor" begin
    @test beam.cₘ    ≈ [2.0, 1.5, 0.0]
    @test beam.rₘ    ≈ [0.8 -0.6 0.0; 0.6 0.8 -0.0; 0.0 0.0 1.0]
    @test beam.ζgp   ≈ [-0.2886751345948129, 0.2886751345948129]
    @test beam.ζnod  ≈ [-0.5, 0.5]
    @test beam.tgₘ   ≈ [4.0, 3.0, 0.0]
    @test beam.tgₑ   ≈ [5.0, 0.0, 0.0]
    @test beam.Nε[1] ≈ [-.2, 0, 0, 0, 0, 0, .2, 0, 0, 0, 0, 0]
    @test beam.Nκ[1][2,2] ≈ -0.2078460969082653
    @test beam.Nκ[1][3,5] ≈ 0.2878460969082653
    @test beam.Nu[1][1,1] ≈ 0.2886751345948129
    @test beam.dL    ≈ [2.5, 2.5]
end

t,SP,dbg  = 0.,(;),(status=:testing,)
x = SVector(1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.)
X = (x,)
U = (SVector{0,𝕣}(),)
A = SVector{0,𝕣}()
R,FB=residual(beam,   X,U,A,t,SP,dbg) 
@testset "residual flex" begin
    @test R        ≈ [1.1323673394544749, 1.0097682842764775, 0.0, 0.0, 0.0, -0.18389858276699608, -1.1323673394544749, -1.0097682842764775, 0.0, 0.0, 0.0, -0.1838985827669961]
    @test FB === nothing
end

###

model           = Model(:TestModel)
node1           = addnode!(model,𝕣[0,0,0])
node2           = addnode!(model,𝕣[1,0,0])
elnod           = [model.nod[n.inod] for n∈[node1,node2]]
mat             = Elements.BeamCrossSection(EA=10.,EI=3.,GJ=4.)
beam            = Elements.EulerBeam3D(elnod;mat,orient2=SVector(0.,1.,0.))
t,SP,dbg  = 0.,(;),(status=:testing,)
x = SVector(0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.)
X = (x,)
U = (SVector{0,𝕣}(),)
A = SVector{0,𝕣}()
R,FB=residual(beam,   X,U,A,t,SP,dbg) 
@testset "residual torsion" begin
    @test R        ≈ [0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0]
    @test FB === nothing
end


end


