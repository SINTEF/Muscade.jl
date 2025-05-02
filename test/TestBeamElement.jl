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
    @test w∂v1 ≈ LinearAlgebra.I#[1 0 0;0 1 0;0 0 1]
    @test v2 ≈ w2
    @test w∂v2 ≈ LinearAlgebra.I#[1 0 0;0 1 0;0 0 1]
end

###
L               = 5
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

    @test beam.yₐ    ≈ [-1/√3,1/√3]
    @test beam.yᵤ    ≈ [-0.7698003589195012,0.7698003589195012]
    @test beam.yᵥ    ≈ [-1/6,-1/6]
    @test beam.κₐ    ≈ [2/L,2/L]
    @test beam.κᵤ    ≈ [0.2771281292110204,-0.2771281292110204]
    @test beam.κᵥ    ≈ [2/L,2/L]

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
    @test R        ≈  [-1,0,0,0,0,0,1,0,0,0,0,0]
    @test FB === nothing
end

x = SVector(0.,0.,0.,0.,0.,0.,0.,0.1,0.,0.,0.,0.); X = (x,)
R,FB=Muscade.residual(beam,   X,U,A,t,SP,dbg) 
@testset "residual flex" begin
    @test R        ≈ [-0.04962809790010789,    -3.604962809790013,     0.0,     0.0,     0.0,    -1.8,     0.04962809790010789,     3.604962809790013,     0.0,     0.0,     0.0,    -1.8]
    @test FB === nothing
end

x = SVector(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.1); X = (x,)
R,FB=Muscade.residual(beam,   X,U,A,t,SP,dbg) 
@testset "residual flex" begin
    @test R        ≈  [ -0.08992502499553628,    1.7970014996429078,    0.0,    0.0,    0.0,    0.5985007498214541,    0.08992502499553628,   -1.7970014996429078,    0.0,    0.0,    0.0,    1.1985007498214544  ]
    @test FB === nothing
end

x = SVector(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.1,0.,0.); X = (x,)
R,FB=Muscade.residual(beam,   X,U,A,t,SP,dbg) 
@testset "residual torsion" begin
    @test R        ≈ [0.0, 0.0, 0.0, -0.40000000000000124, 0.0, 0.0, 0.0, 0.0, 0.0, 0.40000000000000135, 0.0, 0.0]
    @test FB === nothing
end

using Printf
X = (x,x,x)
out = diffed_residual(beam; X,U,A,t,SP)
iλ,ix,iu,ia = 1,2,3,4

R = out.R
K = out.∇R[ix][1]
C = out.∇R[ix][2]
M = out.∇R[ix][3]
H = out.∇R[iu][1]

# using Profile,ProfileView,BenchmarkTools
# mission = :profile
# if  mission == :time
#     @btime out = diffed_residual(beam; X,U,A,t,SP)
# elseif mission == :profile
#     Profile.clear()
#     Profile.@profile for i=1:10000
#         out = diffed_residual(beam; X,U,A,t,SP)
#     end
#     ProfileView.view(fontsize=30);
#     # After clicking on a bar in the flame diagram, you can type warntype_last() and see the result of 
#     # code_warntype for the call represented by that bar.
# end
;

# @printf "\nR\n"
# print_element_array(beam,:X,out.R)    #  R
# @printf "\nK=∂R/∂X₀\n"
# print_element_array(beam,:X,out.∇R[2][1])  # K

# X = (x,x,x)
# out = diffed_residual(beam; X,U,A,t,SP)
# @printf "\nC=∂R/∂X₁\n"
# print_element_array(beam,:X,out.∇R[2][2])  # C
# @printf "\nM=∂R/∂X₂\n"
# print_element_array(beam,:X,out.∇R[2][3])  # M


end



