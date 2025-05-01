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


## Testing residual
EA = 10.
EI = 3.
GJ = 4. 
L =  1.
model           = Model(:TestModel)
node1           = addnode!(model,𝕣[0,0,0])
node2           = addnode!(model,𝕣[L,0,0])
elnod           = [model.nod[n.inod] for n∈[node1,node2]]
mat             = BeamCrossSection(EA=EA,EI=EI,GJ=GJ)
beam            = EulerBeam3D(elnod;mat,orient2=SVector(0.,1.,0.))
t,SP,dbg  = 0.,(;),(status=:testing,)
U = (SVector{0,𝕣}(),)
A = SVector{0,𝕣}()

t1n2 = 0.1;
x = SVector(0.,     0.,     0.,     0.,     0.,     0.,            t1n2,    0.0,    0.,     0.,     0.,     0.); X = (x,)
R,FB=Muscade.residual(beam,   X,U,A,t,SP,dbg) 
@testset "residual tension" begin
    @test R        ≈  [ -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
    @test R        ≈  [ -EA/L*t1n2, 0.0, 0.0, 0.0, 0.0, 0.0, EA/L*t1n2, 0.0, 0.0, 0.0, 0.0, 0.0 ]
    @test FB === nothing
end

t2n2 = 0.01
x = SVector(0.,     0.,     0.,     0.,     0.,     0.,            0.,     t2n2,    0.,     0.,     0.,     0.); X = (x,)
R,FB=Muscade.residual(beam,   X,U,A,t,SP,dbg) 
@testset "residual flex1" begin
    #corect value with new beam @test R        ≈  [0.305626505038752, -3.557508839178609, 0.0, 0.0, 0.0, -1.7940357448412423, -0.305626505038752, 3.557508839178609, 0.0, 0.0, 0.0, -1.7940357448412423]
    @test norm(R - [ 0,     -12*EI/L^3*t2n2,    0,     0,   0,    -6*EI/L^2*t2n2,  0,     12*EI/L^3*t2n2,    0,     0,   0,     -6*EI/L^2*t2n2],Inf)  < 1e-3;
    @test FB === nothing
end

r3n2 = 0.01
x = SVector(0.,     0.,     0.,     0.,     0.,     0.,            0.,     0.,     0.,     0.,     0.,     r3n2); X = (x,)
R,FB=Muscade.residual(beam,   X,U,A,t,SP,dbg) 
@testset "residual flex2" begin
    #correct value with new beam
    #  @test R        ≈  [0.0, 1.8000000000002365, 0.0, 0.0, 0.0, 0.6000000000000143, 0.0, -1.8000000000002365, 0.0, 0.0, 0.0, 1.2000000000002227]
    @test norm(R - [ 0,     6*EI/L^3*r3n2,    0,     0,   0,    2*EI/L*r3n2,  0,     -6*EI/L^3*r3n2,    0,     0,   0,     4*EI/L*r3n2],Inf)  < 1e-3;
    @test FB === nothing
end

r1n2 = 0.1
x = SVector(0.,     0.,     0.,     0.,     0.,     0.,            0.,     0.,     0.,     r1n2,    0.,     0.); X = (x,)
R,FB=Muscade.residual(beam,   X,U,A,t,SP,dbg) 
@testset "residual torsion" begin
    # @test R        ≈ [0.0, 0.0, 0.0, -0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.0, 0.0]
    @test R        ≈ [0.0, 0.0, 0.0, -GJ/L*r1n2, 0.0, 0.0, 0.0, 0.0, 0.0, GJ/L*r1n2, 0.0, 0.0]
    @test FB === nothing
end




displacement =  SVector(0.,0.,0.,0.,0.,0.,  0.,0.,0.,0.,0.,0.); 
velocity =      SVector(0.,0.,0.,0.,0.,0.,  0.,0.,0.,0.,0.,0.); 
acceleration =  SVector(0.,0.,0.,0.,0.,0.,  0.,0.,0.,0.,0.,0.); 
X = (displacement,velocity,acceleration)
out = diffed_residual(beam; X,U,A,t,SP)
iλ,ix,iu,ia = 1,2,3,4
typeof(out.∇R[ix][1])

K = out.∇R[ix][1]
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
    @test K[2,2]        ≈ 12*EI/L^3
    @test K[8,8]        ≈ 12*EI/L^3
    @test K[3,3]        ≈ 12*EI/L^3
    @test K[9,9]        ≈ 12*EI/L^3
    # transverse force induced by tranlation of the opposite node
    @test K[3,9]        ≈ -12*EI/L^3
    @test K[9,3]        ≈ -12*EI/L^3
    @test K[2,8]        ≈ -12*EI/L^3
    @test K[8,2]        ≈ -12*EI/L^3
    # transverse force induced by rotation of the opposite node
    @test K[8,6]        ≈ -6*EI/L^2
    @test K[2,12]       ≈ 6*EI/L^2
    @test K[9,5]        ≈ 6*EI/L^2
    @test K[3,11]       ≈ -6*EI/L^2
    # bending moment induced by rotation of same node
    @test K[5,5]        ≈ 4*EI/L
    @test K[11,11]      ≈ 4*EI/L
    @test K[6,6]        ≈ 4*EI/L
    @test K[12,12]      ≈ 4*EI/L
    # bending moment induced by rotation of opposite node
    @test K[5,11]       ≈ 2*EI/L^2
    @test K[11,5]       ≈ 2*EI/L^2
    @test K[6,12]       ≈ 2*EI/L^2
    @test K[12,6]       ≈ 2*EI/L^2
    # bending moment induced by translation of opposite node
    @test K[5,9]        ≈ 6*EI/L^2
    @test K[11,3]       ≈ -6*EI/L^2
    @test K[6,8]        ≈ -6*EI/L^2
    @test K[12,2]       ≈ 6*EI/L^2
end
@testset "torsional stiffness" begin
    # torsional moment induced by own rotation about element axis
    @test K[4,4]        ≈ GJ/L
    @test K[10,10]      ≈ GJ/L
    # torsional moment induced by rotation of opposite node about element axis
    @test K[4,10]       ≈ -GJ/L
    @test K[10,4]       ≈ -GJ/L
end
@testset "no spurious coupling" begin
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

# R = out.R
# C = out.∇R[ix][2]
# M = out.∇R[ix][3]
# H = out.∇R[iu][1]

# using Printf
# @printf "\nR\n"
# print_element_array(beam,:X,out.R)    #  R
# @printf "\nK=∂R/∂X₀\n"
# print_element_array(beam,:X,out.∇R[2][1])  # K
# @printf "\nC=∂R/∂X₁\n"
# print_element_array(beam,:X,out.∇R[2][2])  # K
# @printf "\nM=∂R/∂X₂\n"
# print_element_array(beam,:X,out.∇R[2][3])  # M


# Cantilever bend with out-of-plane load
# Leading to a three-dimensional response which mobilizes axial force, bending moment and torque.
# Comparison to solutions by Longva (2015) and Crisfield (1990)
R = 100.0;  # Radius of the bend [m]
EI = 833.33;  # Bending stiffness [Nm²]
EA = 1e9;  # Axial stiffness [N]
GJ = 705.;  # Torsional stiffness [Nm²]
Fy = 300.; # then 450 and 600 [N]

nel         = 8
nnodes      = nel+1   
nodeCoord   = hcat(R*cos.(3π/2 .+ ((1:nnodes).-1)/(nnodes-1)*π/2),zeros(Float64,nnod,1),R*(1 .+ sin.(3π/2 .+ ((1:nnodes).-1)/(nnodes-1)*π/2)))
mat         = BeamCrossSection(EA=EA,EI=EI,GJ=GJ)
model       = Model(:TestModel)
nodid       = addnode!(model,nodeCoord)
mesh        = hcat(nodid[1:nnodes-1],nodid[2:nnodes])
eleid       = addelement!(model,EulerBeam3D,mesh;mat=mat,orient2=SVector(0.,1.,0.))
[addelement!(model,Hold,[nodid[1]]  ;field) for field∈[:t1,:t2,:t3,:r1,:r2,:r3]]                                    # Clamp at one end
addelement!(model,DofLoad,[nodid[nnodes]];field=:t2,value=t->-min(1,t)*Fy) ;                                        # Out-of-plane load at other
initialstate    = initialize!(model);
state           = solve(SweepX{0};initialstate,time=[0.,0.])
x_ = getdof(state[2];field=:t1,nodID=[nodid[nnodes]]) #Compare to 58.56 (Longva,2015) or 58.53 (Crisfield, 1990)
y_ = getdof(state[2];field=:t2,nodID=[nodid[nnodes]]) #Compare to 40.47 (Longva,2015) or 40.53 (Crisfield, 1990)
z_ = getdof(state[2];field=:t3,nodID=[nodid[nnodes]]) #Compare to 22.16 (Longva,2015) or 22.14 (Crisfield, 1990)
# Load                300 N                       450 N                       600 N
#                     x,y,z                       x,y,z                       x,y,z
# Disp Longva         58.56, 40.47, 22.18         51.99, 48.72, 18.45         46.91, 53.64, 15.65 
# Disp Crisfield      58.53, 40.53, 22.16         51.93, 48.79, 18.43         46.84, 53.71, 15.61

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
#end
