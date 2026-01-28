module TestBeamElement

using Revise

using Test, Muscade, StaticArrays, LinearAlgebra
using Muscade.Toolbox


L               = 5
model           = Model(:TestModel)
node1           = addnode!(model,𝕣[0,0,0])
node2           = addnode!(model,𝕣[4,3,0])
elnod           = [model.nod[n.inod] for n∈[node1,node2]]
mat             = BarCrossSection(EA=10.,μ=1.)
bar            = Bar3D(elnod;mat)

@testset "constructor" begin
    @test bar.cₘ    ≈ [2.0, 1.5, 0.0]
    @test bar.tgₘ   ≈ [4.0, 3.0, 0.0]
    @test bar.L₀    ≈ 5.0
end

## Testing residual
EA = 10.
L =  2.
μ = 1. 
model           = Model(:TestModel)
node1           = addnode!(model,𝕣[0,0,0])
node2           = addnode!(model,𝕣[L,0,0])
elnod           = [model.nod[n.inod] for n∈[node1,node2]]
mat             = BarCrossSection(EA=EA,μ=μ)
bar            = Bar3D(elnod;mat)
t,SP,dbg  = 0.,(;),(status=:testing,)
U = (SVector{0,𝕣}(),)
A = SVector{0,𝕣}()

t1n2 = 0.1;
x = SVector(0.,0.,0., t1n2,0.0,0.0); X = (x,)
R,FB=Muscade.residual(bar,   X,U,A,t,SP,dbg) 
@testset "residual tension" begin
    @test R        ≈  [ -EA/L*t1n2, 0.0, 0.0, EA/L*t1n2, 0.0, 0.0 ]
    @test FB === nothing
end

displacement =  SVector(0.,0.,0.,  0.,0.,0.); 
velocity     =  SVector(0.,0.,0.,  0.,0.,0.); 
acceleration =  SVector(0.,0.,0.,  0.,0.,0.); 
X            = (displacement,velocity,acceleration)

out = Muscade.diffed_residual(bar; X,U,A,t,SP)
iλ,ix,iu,ia = 1,2,3,4
R = out.R
K = out.∇R[ix][1]
C = out.∇R[ix][2]
M = out.∇R[ix][3]
H = out.∇R[iu][1]

# using Printf
# @printf "\nR\n"
# print_element_array(bar,:X,out.R)    #  R
# @printf "\nK=∂R/∂X₀\n"
# print_element_array(bar,:X,out.∇R[2][1])  # K
# # @printf "\nC=∂R/∂X₁\n"
# # print_element_array(bar,:X,out.∇R[2][2])  # C
# # @printf "\nM=∂R/∂X₂\n"
# # print_element_array(bar,:X,out.∇R[2][3])  # M

@testset "axial stiffness" begin
    # axial force induced by inline displacement of same node
    @test K[1,1]        ≈  EA/L 
    @test K[4,4]        ≈  EA/L 
    # axial force induced by inline displacement of opposite node
    @test K[1,4]        ≈ -EA/L 
    @test K[4,1]        ≈ -EA/L 
end
@testset "spurious stiffness" begin
    # no axial force from anything else than displacements about element axis
    @test norm(K[1, [2,3,5,6]])  ≈ 0.
    @test norm(K[4, [2,3,5,6]])  ≈ 0.
end


# @testset "axial inertia" begin
#     # axial force induced by inline displacement of same node
#     @test M[1,1]        ≈  μ*L/3   
#     @test M[7,7]        ≈  μ*L/3   
#     # axial force induced by inline displacement of opposite node
#     @test M[1,7]        ≈ μ*L/6  
#     @test M[7,1]        ≈ μ*L/6  
# end
# @testset "bending inertia" begin
#     # transverse force induced by translation of same node
#     @test M[8,8]        ≈ 156*μ*L/420    
#     @test M[2,2]        ≈ 156*μ*L/420   
#     @test M[3,3]        ≈ 156*μ*L/420   
#     @test M[9,9]        ≈ 156*μ*L/420   
#     # transverse force induced by translation of the opposite node
#     @test M[3,9]        ≈ 54*μ*L/420   
#     @test M[9,3]        ≈ 54*μ*L/420   
#     @test M[2,8]        ≈ 54*μ*L/420   
#     @test M[8,2]        ≈ 54*μ*L/420     
#     # transverse force induced by rotation of the opposite node
#     @test M[8,6]        ≈  13*μ*L^2/420 
#     @test M[2,12]       ≈ -13*μ*L^2/420 
#     @test M[9,5]        ≈ -13*μ*L^2/420 
#     @test M[3,11]       ≈  13*μ*L^2/420 
#     # # bending moment induced by rotation of same node
#     @test M[5,5]        ≈ 4*μ*L^3/420 
#     @test M[11,11]      ≈ 4*μ*L^3/420 
#     @test M[6,6]        ≈ 4*μ*L^3/420 
#     @test M[12,12]      ≈ 4*μ*L^3/420 
#     # # bending moment induced by rotation of opposite node
#     @test M[5,11]       ≈ -3*μ*L^3/420  
#     @test M[11,5]       ≈ -3*μ*L^3/420  
#     @test M[6,12]       ≈ -3*μ*L^3/420  
#     @test M[12,6]       ≈ -3*μ*L^3/420  
#     # # bending moment induced by translation of opposite node
#     @test M[5,9]        ≈ -13*μ*L^2/420 
#     @test M[11,3]       ≈  13*μ*L^2/420 
#     @test M[6,8]        ≈  13*μ*L^2/420 
#     @test M[12,2]       ≈ -13*μ*L^2/420 
# end
# @testset "torsional inertia" begin
#     # shape function for local roll acceleration not used yet
#     @test M[4,4]        ≈ ι₁*L/4
#     @test M[10,10]      ≈ ι₁*L/4
#     @test M[4,10]       ≈ ι₁*L/4
#     @test M[10,4]       ≈ ι₁*L/4
# end
# @testset "spurious inertia" begin
#     # no axial force from anything else than displacements about element axis
#     @test norm(M[1, [2,3,4,5,6,8,9,10,11,12]])  ≈ 0.
#     @test norm(M[7, [2,3,4,5,6,8,9,10,11,12]])  ≈ 0.
#     # no torsion from anything else than rotations about element axis
#     @test norm(M[4, [1,2,3,5,6,7,8,9,11,12]])   ≈ 0.
#     @test norm(M[10,[1,2,3,5,6,7,8,9,11,12]])   ≈ 0.
#     # no transverse force from anything else than translations in same plane or rotation in orthogonal plane
#     @test norm(M[2, [1,3,4,5,7,9,10,11]])       ≈ 0.
#     @test norm(M[3, [1,2,4,6,7,8,10,12]])       ≈ 0.
#     @test norm(M[8, [1,3,4,5,7,9,10,11]])       ≈ 0.
#     @test norm(M[9, [1,2,4,6,7,8,10,12]])       ≈ 0.
#     # no bending moment force from anything else than translations in same plane or rotation in orthogonal plane
#     @test norm(M[5, [1,2,4,6,7,8,10,12]])       ≈ 0.
#     @test norm(M[6, [1,3,4,5,7,9,10,11]])       ≈ 0.
#     @test norm(M[11,[1,2,4,6,7,8,10,12]])       ≈ 0.
#     @test norm(M[12,[1,3,4,5,7,9,10,11]])       ≈ 0.
# end
# ;

# ## Testing weight
# ## Beam bent upwards with uniform weight (along negative t3)  
# w = 10
# model           = Model(:TestModel)
# node1           = addnode!(model,𝕣[0,0,0])
# node2           = addnode!(model,𝕣[L,0,0])
# elnod           = [model.nod[n.inod] for n∈[node1,node2]]
# mat             = BarCrossSection(EA=EA,EI₂=EI₂,EI₃=EI₃,GJ=GJ,μ=μ,ι₁=ι₁,   w=w)
# bar            = Bar3D(elnod;mat,orient2=SVector(0.,1.,0.))
# x = SVector(0.,     0.,     0.,     0.,     w*L^3/(24*EI₂),     0.,            0.,    0.0,    0.,     0.,     -w*L^3/(24*EI₂),     0.); X = (x,)
# R,FB=Muscade.residual(bar,   X,U,A,t,SP,dbg) 
# @testset "residual weight" begin
#     @test R        ≈  [ 0.0, 0.0, w*L/2, 0.0, 0.0, 0.0, 0.0, 0.0, w*L/2, 0.0, 0.0, 0.0 ]
# end

# ## Testing added mass
# Ca₁ = 1.
# Ca₂ = 2.
# Ca₃ = 3.
# μ = 1. 
# a1,a2,a3 = 4.0,3.0,2.0;
# model           = Model(:TestModel)
# node1           = addnode!(model,𝕣[0,0,0])
# node2           = addnode!(model,𝕣[L,0,0])
# elnod           = [model.nod[n.inod] for n∈[node1,node2]]
# mat             = BarCrossSection(EA=EA ,EI₂=EI₂,EI₃=EI₃,GJ=GJ,μ=μ,ι₁=ι₁,   Ca₁=Ca₁, Ca₂=Ca₂,Ca₃=Ca₃)
# bar            = Bar3D(elnod;mat,orient2=SVector(0.,1.,0.))
# velocity        =  SVector(0.,0.,0.,0.,0.,0.,  0.,0.,0.,0.,0.,0.); 

# @testset "residual addded mass" begin
#     displacement    =  SVector(0.,0.,0.,0.,0.,0.,  0.,0.,0.,0.,0.,0.); 
#     acceleration =  SVector(a1,0.,0.,0.,0.,0.,  a1,0.,0.,0.,0.,0.); 
#     X = (displacement,velocity,acceleration); 
#     R,FB=Muscade.residual(bar,   X,U,A,t,SP,dbg)     
#     @test R        ≈  [ (μ+Ca₁)*a1*L/2, 0., 0.,    0.0, 0.0, 0.0,  (μ+Ca₁)*a1*L/2, 0., 0.,  0.0, 0.0, 0.0 ]

#     displacement =  SVector(0.,     0.,     0.,     0.,     0.,     -(μ+Ca₂)*a2*L^3/(24*EI₃),     0.,    0.,    0.,     0.,      0.,      (μ+Ca₂)*a2*L^3/(24*EI₃))
#     acceleration =  SVector(0.,a2,0.,0.,0.,0.,  0.,a2,0.,0.,0.,0.); 
#     X = (displacement,velocity,acceleration); 
#     R,FB=Muscade.residual(bar,   X,U,A,t,SP,dbg) 
#     @test R        ≈  [ 0., (μ+Ca₂)*a2*L/2, 0.,    0.0, 0.0, 0.0,  0., (μ+Ca₂)*a2*L/2, 0.,  0.0, 0.0, 0.0 ]

#     displacement =  SVector(0.,     0.,     0.,     0.,     (μ+Ca₃)*a3*L^3/(24*EI₂),     0.,     0.,    0.,    0.,     0.,      -(μ+Ca₃)*a3*L^3/(24*EI₂),      0.)
#     acceleration =  SVector(0.,0.,a3,0.,0.,0.,  0.,0.,a3,0.,0.,0.); 
#     X = (displacement,velocity,acceleration); 
#     R,FB=Muscade.residual(bar,   X,U,A,t,SP,dbg) 
#     @test R        ≈  [ 0., 0., (μ+Ca₃)*a3*L/2,    0.0, 0.0, 0.0,  0., 0., (μ+Ca₃)*a3*L/2,  0.0, 0.0, 0.0 ]

# end

# ## Testing damping
# Cl₁ = 1.
# Cl₂ = 2.
# Cl₃ = 3.
# Cq₁ = .1
# Cq₂ = .2
# Cq₃ = .3
# μ = 1. 
# v1,v2,v3 = 1.0,1.1,0.1;
# model           = Model(:TestModel)
# node1           = addnode!(model,𝕣[0,0,0])
# node2           = addnode!(model,𝕣[L,0,0])
# elnod           = [model.nod[n.inod] for n∈[node1,node2]]
# mat             = BarCrossSection(EA=EA ,EI₂=EI₂,EI₃=EI₃,GJ=GJ,μ=μ,ι₁=ι₁,   Cl₁=Cl₁, Cl₂=Cl₂,Cl₃=Cl₃, Cq₁=Cq₁, Cq₂=Cq₂,Cq₃=Cq₃)
# bar            = Bar3D(elnod;mat,orient2=SVector(0.,1.,0.))
# acceleration =  SVector(0,0.,0.,0.,0.,0.,  0.,0.,0.,0.,0.,0.); 
# @testset "residual damping" begin
#     displacement    =  SVector(0.,0.,0.,0.,0.,0.,  0.,0.,0.,0.,0.,0.); 
#     velocity        =  SVector(v1,0.,0.,0.,0.,0.,  v1,0.,0.,0.,0.,0.); 
#     X = (displacement,velocity,acceleration); 
#     R,FB=Muscade.residual(bar,   X,U,A,t,SP,dbg) 
#     @test R        ≈  [ (Cl₁+Cq₁*abs(v1))*v1*L/2, 0., 0.,    0.0, 0.0, 0.0,  (Cl₁+Cq₁*abs(v1))*v1*L/2, 0., 0.,  0.0, 0.0, 0.0 ]

#     displacement =  SVector(0.,     0.,     0.,     0.,     0.,     -(Cl₂+Cq₂*abs(v2))*v2*L^3/(24*EI₃),     0.,    0.,    0.,     0.,      0.,      (Cl₂+Cq₂*abs(v2))*v2*L^3/(24*EI₃))
#     velocity =      SVector(0.,v2,0.,0.,0.,0.,  0.,v2,0.,0.,0.,0.); 
#     X = (displacement,velocity,acceleration); 
#     R,FB=Muscade.residual(bar,   X,U,A,t,SP,dbg) 
#     @test R        ≈  [ 0., (Cl₂+Cq₂*abs(v2))*v2*L/2, 0.,    0.0, 0.0, 0.0,  0., (Cl₂+Cq₂*abs(v2))*v2*L/2, 0.,  0.0, 0.0, 0.0 ]

#     displacement =  SVector(0.,     0.,     0.,     0.,     (Cl₃+Cq₃*abs(v3))*v3*L^3/(24*EI₂),     0.,     0.,    0.,    0.,     0.,      -(Cl₃+Cq₃*abs(v3))*v3*L^3/(24*EI₂),      0.)
#     velocity =      SVector(0.,0.,v3,0.,0.,0.,  0.,0.,v3,0.,0.,0.); 
#     X = (displacement,velocity,acceleration); 
#     R,FB=Muscade.residual(bar,   X,U,A,t,SP,dbg)
#     @test R        ≈  [ 0., 0., (Cl₃+Cq₃*abs(v3))*v3*L/2,    0.0, 0.0, 0.0,  0., 0., (Cl₃+Cq₃*abs(v3))*v3*L/2,  0.0, 0.0, 0.0 ]
# end

end