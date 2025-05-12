#module TestBeamElements

using Test, Muscade, StaticArrays, LinearAlgebra
include("../examples/BeamElements.jl")

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
mat             = BeamCrossSection(EA=10.,EI=3.,GJ=4.)

beam            = EulerBeam3D(elnod;mat,orient2=SVector(0.,1.,0.))

# @testset "constructor" begin
#     @test beam.cₘ    ≈ [2.0, 1.5, 0.0]
#     @test beam.rₘ    ≈ [0.8 -0.6 0.0; 0.6 0.8 -0.0; 0.0 0.0 1.0]
#     @test beam.ζgp   ≈ [-0.2886751345948129, 0.2886751345948129]
#     @test beam.ζnod  ≈ [-0.5, 0.5]
#     @test beam.tgₘ   ≈ [4.0, 3.0, 0.0]
#     @test beam.tgₑ   ≈ [5.0, 0.0, 0.0]

#     @test beam.yₐ    ≈ [-1/√3,1/√3]
#     @test beam.yᵤ    ≈ [-0.7698003589195012,0.7698003589195012]
#     @test beam.yᵥ    ≈ [-1/6,-1/6]
#     @test beam.κₐ    ≈ [2/L,2/L]
#     @test beam.κᵤ    ≈ [0.2771281292110204,-0.2771281292110204]
#     @test beam.κᵥ    ≈ [2/L,2/L]

#     @test beam.dL    ≈ [2.5, 2.5]
# end


## Testing residual
EA = 10.
EI = 3.
GJ = 4. 
L =  2.
μ = 1. # currently hard-coded in beam element [kg/m]
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
    # @test R        ≈  [ -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
    @test R        ≈  [ -EA/L*t1n2, 0.0, 0.0, 0.0, 0.0, 0.0, EA/L*t1n2, 0.0, 0.0, 0.0, 0.0, 0.0 ]
    @test FB === nothing
end

t2n2 = 0.01
x = SVector(0.,     0.,     0.,     0.,     0.,     0.,            0.,     t2n2,    0.,     0.,     0.,     0.); X = (x,)
R,FB=Muscade.residual(beam,   X,U,A,t,SP,dbg) 
@testset "residual flex1" begin
    #corect value with new beam @test R        ≈  [0.305626505038752, -3.557508839178609, 0.0, 0.0, 0.0, -1.7940357448412423, -0.305626505038752, 3.557508839178609, 0.0, 0.0, 0.0, -1.7940357448412423]
    @test R ≈ [ 0,     -12*EI/L^3*t2n2,    0,     0,   0,    -6*EI/L^2*t2n2,  0,     12*EI/L^3*t2n2,    0,     0,   0,     -6*EI/L^2*t2n2]  atol=1e-2;
    @test FB === nothing
end

r3n2 = 0.01
x = SVector(0.,     0.,     0.,     0.,     0.,     0.,            0.,     0.,     0.,     0.,     0.,     r3n2); X = (x,)
R,FB=Muscade.residual(beam,   X,U,A,t,SP,dbg) 
@testset "residual flex2" begin
    #correct value with new beam
    #  @test R        ≈  [0.0, 1.8000000000002365, 0.0, 0.0, 0.0, 0.6000000000000143, 0.0, -1.8000000000002365, 0.0, 0.0, 0.0, 1.2000000000002227]
    @test R ≈ [ 0,     6*EI/L^2*r3n2,    0,     0,   0,    2*EI/L*r3n2,  0,     -6*EI/L^2*r3n2,    0,     0,   0,     4*EI/L*r3n2] atol=1e-2
    @test FB === nothing
end

r1n2 = 0.1
x = SVector(0.,     0.,     0.,     0.,     0.,     0.,            0.,     0.,     0.,     r1n2,    0.,     0.); X = (x,)
R,FB=Muscade.residual(beam,   X,U,A,t,SP,dbg) 
@testset "residual torsion" begin
    # @test R        ≈ [0.0, 0.0, 0.0, -0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.0, 0.0]
    @test R ≈ [0.0, 0.0, 0.0, -GJ/L*r1n2, 0.0, 0.0, 0.0, 0.0, 0.0, GJ/L*r1n2, 0.0, 0.0]
    @test FB === nothing
end

displacement =  SVector(0.,0.,0.,0.,0.,0.,  0.,0.,0.,0.,0.,0.); 
velocity     =  SVector(0.,0.,0.,0.,0.,0.,  0.,0.,0.,0.,0.,0.); 
acceleration =  SVector(0.,0.,0.,0.,0.,0.,  0.,0.,0.,0.,0.,0.); 
X = (displacement,velocity,acceleration)

out = diffed_residual(beam; X,U,A,t,SP)
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
    @test K[2,2]        ≈ 12*EI/L^3  
    @test K[8,8]        ≈ 12*EI/L^3  
    @test K[3,3]        ≈ 12*EI/L^3  
    @test K[9,9]        ≈ 12*EI/L^3  
    # transverse force induced by translation of the opposite node
    @test K[3,9]        ≈ -12*EI/L^3  
    @test K[9,3]        ≈ -12*EI/L^3  
    @test K[2,8]        ≈ -12*EI/L^3  
    @test K[8,2]        ≈ -12*EI/L^3  
    # transverse force induced by rotation of the opposite node
    @test K[8,6]        ≈ -6*EI/L^2  
    @test K[2,12]       ≈  6*EI/L^2   
    @test K[9,5]        ≈  6*EI/L^2  
    @test K[3,11]       ≈ -6*EI/L^2  
    # bending moment induced by rotation of same node
    @test K[5,5]        ≈ 4*EI/L  
    @test K[11,11]      ≈ 4*EI/L  
    @test K[6,6]        ≈ 4*EI/L  
    @test K[12,12]      ≈ 4*EI/L  
    # bending moment induced by rotation of opposite node
    @test K[5,11]       ≈ 2*EI/L  
    @test K[11,5]       ≈ 2*EI/L  
    @test K[6,12]       ≈ 2*EI/L  
    @test K[12,6]       ≈ 2*EI/L  
    # bending moment induced by translation of opposite node
    @test K[5,9]        ≈  6*EI/L^2  
    @test K[11,3]       ≈ -6*EI/L^2  
    @test K[6,8]        ≈ -6*EI/L^2  
    @test K[12,2]       ≈  6*EI/L^2  
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


@testset "axial stiffness" begin
    # axial force induced by inline displacement of same node
    @test K[1,1]        ≈  EA/L   
    @test K[7,7]        ≈  EA/L   
    # axial force induced by inline displacement of opposite node
    @test K[1,7]        ≈ -EA/L  
    @test K[7,1]        ≈ -EA/L  
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
@testset "torsional stiffness" begin
    # torsional moment induced by own rotation about element axis
    @test K[4,4]        ≈ GJ/L 
    @test K[10,10]      ≈ GJ/L 
    # torsional moment induced by rotation of opposite node about element axis
    @test K[4,10]       ≈ -GJ/L 
    @test K[10,4]       ≈ -GJ/L 
end
;
# # Cantilever bend, with out-of-plane load leading to a three-dimensional response mobilizing axial force, bending moment and torque.
# R = 100.0;  # Radius of the bend [m]
# EI = 833.33e3;  # Bending stiffness [Nm²]
# EA = 1e9;  # Axial stiffness [N]
# GJ = 705e3;  # Torsional stiffness [Nm²]
# Fy = 300.; # 300 then 450 and 600 [N]
# nel         = 8
# nnodes      = nel+1   
# nodeCoord   = hcat( 0 .+ R*cos.(3π/2 .+ ((1:nnodes).-1)/(nnodes-1)*π/4),
#                     0 .+ zeros(Float64,nnodes,1),
#                     R .+ R*sin.(3π/2 .+ ((1:nnodes).-1)/(nnodes-1)*π/4))
# mat         = BeamCrossSection(EA=EA,EI=EI,GJ=GJ)
# model       = Model(:TestModel)
# nodid       = addnode!(model,nodeCoord)
# mesh        = hcat(nodid[1:nnodes-1],nodid[2:nnodes])
# eleid       = addelement!(model,EulerBeam3D,mesh;mat=mat,orient2=SVector(0.,1.,0.))
# [addelement!(model,Hold,[nodid[1]]  ;field) for field∈[:t1,:t2,:t3,:r1,:r2,:r3]]                                   # Clamp at one end
# function load(t)
#     t<=1. ? load=t*300. : 
#     t>1. && t<=2. ? load=300. +(t-1)*150. :
#     load=450. +(t-2)*150. 
# end
# addelement!(model,DofLoad,[nodid[nnodes]];field=:t2,value=t->load(t)) ;                                        # Force along axis2 at other
# initialstate    = initialize!(model);
# loadSteps = [0.,1.,2.,3.];
# nLoadSteps = length(loadSteps)
# state           = solve(SweepX{0};initialstate,time=loadSteps,verbose=true,maxΔx=1e-9)
# # Plot of the beam configuration 
# x_ = [getdof(state[idxLoad];field=:t1,nodID=nodid[1:nnodes]) for idxLoad ∈ 1:nLoadSteps] 
# y_ = [getdof(state[idxLoad];field=:t2,nodID=nodid[1:nnodes]) for idxLoad ∈ 1:nLoadSteps] 
# z_ = [getdof(state[idxLoad];field=:t3,nodID=nodid[1:nnodes]) for idxLoad ∈ 1:nLoadSteps] 

# using GLMakie
# currentDir = @__DIR__
# fig      = Figure(size = (1000,1000))
# ax = Axis3(fig[1,1],xlabel="x [m]", ylabel="y [m]", zlabel="z [m]",aspect=:equal)
# for idxLoad ∈ 1:nLoadSteps
#     lines!(ax,nodeCoord[:,1]+x_[idxLoad][:], nodeCoord[:,2]+y_[idxLoad][:] , nodeCoord[:,3]+z_[idxLoad][:]                  , label="F="*string(load(loadSteps[idxLoad]))*" N");
# end
# xlims!(ax, 0,71); ylims!(ax, 0,71); zlims!(ax, 0,71); axislegend()
# save(normpath(joinpath(currentDir,"beamTest1.png")),fig)
# Comparison to solutions by Longva (2015) and Crisfield (1990)
# Load                300 N                       450 N                       600 N
#                     x,y,z                       x,y,z                       x,y,z
# Disp Longva         58.56, 40.47, 22.18         51.99, 48.72, 18.45         46.91, 53.64, 15.65 
# Disp Crisfield      58.53, 40.53, 22.16         51.93, 48.79, 18.43         46.84, 53.71, 15.61

# height = [  nodeCoord[end,1]+x_[2][end], 58.56, 58.53, nodeCoord[end,1]+x_[3][end], 51.99, 51.93, nodeCoord[end,1]+x_[4][end], 46.91, 46.84, 
#             nodeCoord[end,2]+y_[2][end], 40.47, 40.53, nodeCoord[end,2]+y_[3][end], 48.72, 48.79, nodeCoord[end,2]+y_[4][end], 53.64, 53.71,
#             nodeCoord[end,3]+z_[2][end], 22.18, 22.16, nodeCoord[end,3]+z_[3][end], 18.45, 18.43, nodeCoord[end,3]+z_[4][end], 15.65, 15.61]
# @show height

# colors = [:red, :blue, :green]
# tbl = (cat = [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9],
#        height,
#        grp = [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3])
# fig = Figure(size = (2000, 400))
# ax = Axis(fig[1,1], xticks = (1:9, ["x @ 300 N", "x @ 450 N", "x @ 600 N", "y @ 300 N", "y @ 450 N", "y @ 600 N", "z @ 300 N", "z @ 450 N", "z @ 600 N"]),
#                     ylabel = "End node coordinate [m]")
# barplot!(ax, tbl.cat, tbl.height,
#         dodge = tbl.grp,
#         color = colors[tbl.grp],
#         bar_labels = :y
#         )
# labels = ["Muscade", "Longva (2015)", "Crisfield (1990)"]
# elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
# title = "Method"
# Legend(fig[1,2], elements, labels, title)
# save(normpath(joinpath(currentDir,"beamTest1-comparison.png")),fig)

# # using Profile,ProfileView,BenchmarkTools
# # mission = :profile
# # if  mission == :time
# #     @btime out = diffed_residual(beam; X,U,A,t,SP)
# # elseif mission == :profile
# #     Profile.clear()
# #     Profile.@profile for i=1:10000
# #         out = diffed_residual(beam; X,U,A,t,SP)
# #     end
# #     ProfileView.view(fontsize=30);
# #     # After clicking on a bar in the flame diagram, you can type warntype_last() and see the result of 
# #     # code_warntype for the call represented by that bar.
# # end
# ;
#end