module TestBeamElement

using Test, Muscade, StaticArrays, LinearAlgebra
include("../examples/BeamElement.jl")

a = SA[1,0,0]
b = SA[0,1,1]
r = adjust(a,b)
R = Rodrigues(r)
u = R*a

# v1      = variate{1,3}(SA[.1,.2,.3])
# M1      = Rodrigues(v1)
# w1,w∂v1 = value_∂{1,3}(Rodrigues⁻¹(M1))

# v2      = variate{1,3}(SA[1e-7,2e-7,1e-8])
# M2      = Rodrigues(v2)
# w2,w∂v2 = value_∂{1,3}(Rodrigues⁻¹(M2))


# @testset "rotations" begin
#     @test r ≈ [0.0, -1.1107207345395913, 1.1107207345395913]
#     @test u ≈ [2.220446049250313e-16, 0.7071067811865476, 0.7071067811865476]
#     @test v1 ≈ w1
#     @test w∂v1 ≈ LinearAlgebra.I#[1 0 0;0 1 0;0 0 1]
#     @test v2 ≈ w2
#     @test w∂v2 ≈ LinearAlgebra.I#[1 0 0;0 1 0;0 0 1]
# end

###
L               = 5
model           = Model(:TestModel)
node1           = addnode!(model,𝕣[0,0,0])
node2           = addnode!(model,𝕣[4,3,0])
elnod           = [model.nod[n.inod] for n∈[node1,node2]]
mat             = BeamCrossSection(EA=10.,EI₂=3.,EI₃=3.,GJ=4.,μ=1.,ι₁=1.0)

beam            = EulerBeam3D(elnod;mat,orient2=SVector(0.,1.,0.))

@testset "constructor" begin
    @test beam.cₘ    ≈ [2.0, 1.5, 0.0]
    @test beam.rₘ    ≈ [0.8 -0.6 0.0; 0.6 0.8 -0.0; 0.0 0.0 1.0]
    # @test beam.ζgp   ≈ [-0.2886751345948129, 0.2886751345948129]
    @test beam.ζnod  ≈ [-0.5, 0.5]
    @test beam.tgₘ   ≈ [4.0, 3.0, 0.0]
    @test beam.tgₑ   ≈ [5.0, 0.0, 0.0]

    # @test beam.yₐ    ≈ [-1/√3,1/√3]
    # @test beam.yᵤ    ≈ [-0.7698003589195012,0.7698003589195012]
    # @test beam.yᵥ    ≈ [-1/6,-1/6]
    # @test beam.κₐ    ≈ [2/L,2/L]
    # @test beam.κᵤ    ≈ [0.2771281292110204,-0.2771281292110204]
    # @test beam.κᵥ    ≈ [2/L,2/L]

    # @test beam.dL    ≈ [2.5, 2.5]
end


## Testing residual
EA = 10.
EI₂ = 3.
EI₃ = 2.
GJ = 4. 
L =  2.
μ = 1. 
ι₁ = 5.
model           = Model(:TestModel)
node1           = addnode!(model,𝕣[0,0,0])
node2           = addnode!(model,𝕣[L,0,0])
elnod           = [model.nod[n.inod] for n∈[node1,node2]]
mat             = BeamCrossSection(EA=EA,EI₂=EI₂,EI₃=EI₃,GJ=GJ,μ=μ,ι₁=ι₁)
beam            = EulerBeam3D(elnod;mat,orient2=SVector(0.,1.,0.))
t,SP,dbg  = 0.,(;),(status=:testing,)
U = (SVector{0,𝕣}(),)
A = SVector{0,𝕣}()

t1n2 = 0.1;
x = SVector(0.,     0.,     0.,     0.,     0.,     0.,            t1n2,    0.0,    0.,     0.,     0.,     0.); X = (x,)
R,FB=Muscade.residual(beam,   X,U,A,t,SP,dbg) 
@testset "residual tension" begin
    @test R        ≈  [ -EA/L*t1n2, 0.0, 0.0, 0.0, 0.0, 0.0, EA/L*t1n2, 0.0, 0.0, 0.0, 0.0, 0.0 ]
    @test FB === nothing
end

t2n2 = 0.01
x = SVector(0.,     0.,     0.,     0.,     0.,     0.,            0.,     t2n2,    0.,     0.,     0.,     0.); X = (x,)
R,FB=Muscade.residual(beam,   X,U,A,t,SP,dbg) 
@testset "residual flex2t" begin
    @test R ≈ [ 0,     -12*EI₃/L^3*t2n2,    0,     0,   0,    -6*EI₃/L^2*t2n2,  0,     12*EI₃/L^3*t2n2,    0,     0,   0,     -6*EI₃/L^2*t2n2]  atol=1e-2;
    @test FB === nothing
end

r3n2 = 0.01
x = SVector(0.,     0.,     0.,     0.,     0.,     0.,            0.,     0.,     0.,     0.,     0.,     r3n2); X = (x,)
R,FB=Muscade.residual(beam,   X,U,A,t,SP,dbg) 
@testset "residual flex2r" begin
    @test R ≈ [ 0,     6*EI₃/L^2*r3n2,    0,     0,   0,    2*EI₃/L*r3n2,  0,     -6*EI₃/L^2*r3n2,    0,     0,   0,     4*EI₃/L*r3n2] atol=1e-2
    @test FB === nothing
end

r1n2 = 0.1
x = SVector(0.,     0.,     0.,     0.,     0.,     0.,            0.,     0.,     0.,     r1n2,    0.,     0.); X = (x,)
R,FB=Muscade.residual(beam,   X,U,A,t,SP,dbg) 
@testset "residual torsion" begin
    @test R ≈ [0.0, 0.0, 0.0, -GJ/L*r1n2, 0.0, 0.0, 0.0, 0.0, 0.0, GJ/L*r1n2, 0.0, 0.0]
    @test FB === nothing
end

displacement =  SVector(0.,0.,0.,0.,0.,0.,  0.,0.,0.,0.,0.,0.); 
velocity     =  SVector(0.,0.,0.,0.,0.,0.,  0.,0.,0.,0.,0.,0.); 
acceleration =  SVector(0.,0.,0.,0.,0.,0.,  0.,0.,0.,0.,0.,0.); 
X            = (displacement,velocity,acceleration)

out = Muscade.diffed_residual(beam; X,U,A,t,SP)
iλ,ix,iu,ia = 1,2,3,4
R = out.R
K = out.∇R[ix][1]
C = out.∇R[ix][2]
M = out.∇R[ix][3]
H = out.∇R[iu][1]

# using Printf
# @printf "\nR\n"
# print_element_array(beam,:X,out.R)    #  R
# @printf "\nK=∂R/∂X₀\n"
# print_element_array(beam,:X,out.∇R[2][1])  # K
# @printf "\nC=∂R/∂X₁\n"
# print_element_array(beam,:X,out.∇R[2][2])  # C
# @printf "\nM=∂R/∂X₂\n"
# print_element_array(beam,:X,out.∇R[2][3])  # M

@testset "axial stiffness" begin
    # axial force induced by inline displacement of same node
    @test K[1,1]        ≈  EA/L 
    @test K[7,7]        ≈  EA/L 
    # axial force induced by inline displacement of opposite node
    @test K[1,7]        ≈ -EA/L 
    @test K[7,1]        ≈ -EA/L 
end
@testset "bending stiffness" begin
    # transverse force induced by translation of same node
    @test K[2,2]        ≈ 12*EI₃/L^3  
    @test K[8,8]        ≈ 12*EI₃/L^3  
    @test K[3,3]        ≈ 12*EI₂/L^3  
    @test K[9,9]        ≈ 12*EI₂/L^3  
    # transverse force induced by translation of the opposite node
    @test K[3,9]        ≈ -12*EI₂/L^3  
    @test K[9,3]        ≈ -12*EI₂/L^3  
    @test K[2,8]        ≈ -12*EI₃/L^3  
    @test K[8,2]        ≈ -12*EI₃/L^3  
    # transverse force induced by rotation of the opposite node
    @test K[8,6]        ≈ -6*EI₃/L^2  
    @test K[2,12]       ≈  6*EI₃/L^2   
    @test K[9,5]        ≈  6*EI₂/L^2  
    @test K[3,11]       ≈ -6*EI₂/L^2  
    # bending moment induced by rotation of same node
    @test K[5,5]        ≈ 4*EI₂/L  
    @test K[11,11]      ≈ 4*EI₂/L  
    @test K[6,6]        ≈ 4*EI₃/L  
    @test K[12,12]      ≈ 4*EI₃/L  
    # bending moment induced by rotation of opposite node
    @test K[5,11]       ≈ 2*EI₂/L  
    @test K[11,5]       ≈ 2*EI₂/L  
    @test K[6,12]       ≈ 2*EI₃/L  
    @test K[12,6]       ≈ 2*EI₃/L  
    # bending moment induced by translation of opposite node
    @test K[5,9]        ≈  6*EI₂/L^2  
    @test K[11,3]       ≈ -6*EI₂/L^2  
    @test K[6,8]        ≈ -6*EI₃/L^2  
    @test K[12,2]       ≈  6*EI₃/L^2  
end
@testset "torsional stiffness" begin
    # torsional moment induced by own rotation about element axis
    @test K[4,4]        ≈ GJ/L  
    @test K[10,10]      ≈ GJ/L  
    # torsional moment induced by rotation of opposite node about element axis
    @test K[4,10]       ≈ -GJ/L  
    @test K[10,4]       ≈ -GJ/L  
end
@testset "spurious stiffness" begin
    # no axial force from anything else than displacements about element axis
    @test norm(K[1, [2,3,4,5,6,8,9,10,11,12]])  ≈ 0.
    @test norm(K[7, [2,3,4,5,6,8,9,10,11,12]])  ≈ 0.
    # no torsion from anything else than rotations about element axis
    @test norm(K[4, [1,2,3,5,6,7,8,9,11,12]])   ≈ 0.
    @test norm(K[10,[1,2,3,5,6,7,8,9,11,12]])   ≈ 0.
    # no transverse force from anything else than translations in same plane or rotation in orthogonal plane
    @test norm(K[2, [1,3,4,5,7,9,10,11]])       ≈ 0.
    @test norm(K[3, [1,2,4,6,7,8,10,12]])       ≈ 0.
    @test norm(K[8, [1,3,4,5,7,9,10,11]])       ≈ 0.
    @test norm(K[9, [1,2,4,6,7,8,10,12]])       ≈ 0.
    # no bending moment force from anything else than translations in same plane or rotation in orthogonal plane
    @test norm(K[5, [1,2,4,6,7,8,10,12]])       ≈ 0.
    @test norm(K[6, [1,3,4,5,7,9,10,11]])       ≈ 0.
    @test norm(K[11,[1,2,4,6,7,8,10,12]])       ≈ 0.
    @test norm(K[12,[1,3,4,5,7,9,10,11]])       ≈ 0.
end


@testset "axial inertia" begin
    # axial force induced by inline displacement of same node
    @test M[1,1]        ≈  μ*L/3   
    @test M[7,7]        ≈  μ*L/3   
    # axial force induced by inline displacement of opposite node
    @test M[1,7]        ≈ μ*L/6  
    @test M[7,1]        ≈ μ*L/6  
end
@testset "bending inertia" begin
    # transverse force induced by translation of same node
    @test M[8,8]        ≈ 156*μ*L/420    
    @test M[2,2]        ≈ 156*μ*L/420   
    @test M[3,3]        ≈ 156*μ*L/420   
    @test M[9,9]        ≈ 156*μ*L/420   
    # transverse force induced by translation of the opposite node
    @test M[3,9]        ≈ 54*μ*L/420   
    @test M[9,3]        ≈ 54*μ*L/420   
    @test M[2,8]        ≈ 54*μ*L/420   
    @test M[8,2]        ≈ 54*μ*L/420     
    # transverse force induced by rotation of the opposite node
    @test M[8,6]        ≈  13*μ*L^2/420 
    @test M[2,12]       ≈ -13*μ*L^2/420 
    @test M[9,5]        ≈ -13*μ*L^2/420 
    @test M[3,11]       ≈  13*μ*L^2/420 
    # # bending moment induced by rotation of same node
    @test M[5,5]        ≈ 4*μ*L^3/420 
    @test M[11,11]      ≈ 4*μ*L^3/420 
    @test M[6,6]        ≈ 4*μ*L^3/420 
    @test M[12,12]      ≈ 4*μ*L^3/420 
    # # bending moment induced by rotation of opposite node
    @test M[5,11]       ≈ -3*μ*L^3/420  
    @test M[11,5]       ≈ -3*μ*L^3/420  
    @test M[6,12]       ≈ -3*μ*L^3/420  
    @test M[12,6]       ≈ -3*μ*L^3/420  
    # # bending moment induced by translation of opposite node
    @test M[5,9]        ≈ -13*μ*L^2/420 
    @test M[11,3]       ≈  13*μ*L^2/420 
    @test M[6,8]        ≈  13*μ*L^2/420 
    @test M[12,2]       ≈ -13*μ*L^2/420 
end
@testset "torsional inertia" begin
    # shape function for local roll acceleration not used yet
    @test M[4,4]        ≈ ι₁*L/4
    @test M[10,10]      ≈ ι₁*L/4
    @test M[4,10]       ≈ ι₁*L/4
    @test M[10,4]       ≈ ι₁*L/4
end
@testset "spurious inertia" begin
    # no axial force from anything else than displacements about element axis
    @test norm(M[1, [2,3,4,5,6,8,9,10,11,12]])  ≈ 0.
    @test norm(M[7, [2,3,4,5,6,8,9,10,11,12]])  ≈ 0.
    # no torsion from anything else than rotations about element axis
    @test norm(M[4, [1,2,3,5,6,7,8,9,11,12]])   ≈ 0.
    @test norm(M[10,[1,2,3,5,6,7,8,9,11,12]])   ≈ 0.
    # no transverse force from anything else than translations in same plane or rotation in orthogonal plane
    @test norm(M[2, [1,3,4,5,7,9,10,11]])       ≈ 0.
    @test norm(M[3, [1,2,4,6,7,8,10,12]])       ≈ 0.
    @test norm(M[8, [1,3,4,5,7,9,10,11]])       ≈ 0.
    @test norm(M[9, [1,2,4,6,7,8,10,12]])       ≈ 0.
    # no bending moment force from anything else than translations in same plane or rotation in orthogonal plane
    @test norm(M[5, [1,2,4,6,7,8,10,12]])       ≈ 0.
    @test norm(M[6, [1,3,4,5,7,9,10,11]])       ≈ 0.
    @test norm(M[11,[1,2,4,6,7,8,10,12]])       ≈ 0.
    @test norm(M[12,[1,3,4,5,7,9,10,11]])       ≈ 0.
end
;
# # using Profile,ProfileView,BenchmarkTools
# # mission = :profile
# # if  mission == :time
# #     @btime out = Muscade.diffed_residual(beam; X,U,A,t,SP)
# # elseif mission == :profile
# #     Profile.clear()
# #     Profile.@profile for i=1:10000
# #         out = Muscade.diffed_residual(beam; X,U,A,t,SP)
# #     end
# #     ProfileView.view(fontsize=30);
# #     # After clicking on a bar in the flame diagram, you can type warntype_last() and see the result of 
# #     # code_warntype for the call represented by that bar.
# # end
# ;
end