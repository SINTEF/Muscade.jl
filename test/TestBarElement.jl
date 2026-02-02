module TestBarElement

using Test, Muscade, StaticArrays, LinearAlgebra
using Muscade.Toolbox


EA = 10.
L₀ =  2.
μ = 1. 
model           = Model(:TestModel)
node1           = addnode!(model,𝕣[0,0,0])
node2           = addnode!(model,𝕣[L₀,0,0])
elnod           = [model.nod[n.inod] for n∈[node1,node2]]
mat             = AxisymmetricBarCrossSection(EA=EA,μ=μ)
bar            = Bar3D(elnod;mat)

@testset "constructor" begin
    @test bar.cₘ    ≈ [L₀/2, 0.0, 0.0]
    @test bar.tgₘ   ≈ [L₀, 0.0, 0.0]
    @test bar.L₀    ≈ L₀
    @test bar.wgp   ≈ [0.34785484513745385, 0.6521451548625462, 0.6521451548625462, 0.34785484513745385]
    @test bar.ζgp   ≈ [-0.4305681557970263, -0.16999052179242816, 0.16999052179242816, 0.4305681557970263]
    @test bar.ζnod  ≈ [-0.5, 0.5]
    @test bar.ψ₁    ≈ [0.9305681557970262, 0.6699905217924281, 0.33000947820757187, 0.06943184420297371]
    @test bar.ψ₂    ≈ [0.06943184420297371, 0.33000947820757187, 0.6699905217924281, 0.9305681557970262]
end

Δx = .1
x = SVector(0.,0.,0., Δx,0.0,0.0); X = (x,)
U = (SVector{0,𝕣}(),)
A = SVector{0,𝕣}()
t,SP,dbg  = 0.,(;),(status=:testing,)
R,FB=Muscade.residual(bar,   X,U,A,t,SP,dbg) 

@testset "residual tension" begin
    @test R        ≈  [ -EA/L₀*Δx, 0.0, 0.0, EA/L₀*Δx, 0.0, 0.0 ]
    @test FB === nothing
end

displacement =  SVector(0.,0.,0.,  Δx,0.,0.); 
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
    @test K[1,1]        ≈  EA/L₀ 
    @test K[4,4]        ≈  EA/L₀ 
    @test K[1,4]        ≈ -EA/L₀ 
    @test K[4,1]        ≈ -EA/L₀ 
end
@testset "transverse stiffness" begin
    L = L₀+Δx
    kₜ = (EA/L₀)*Δx/L
    @test K[2,2]        ≈  kₜ
    @test K[3,3]        ≈  kₜ
    @test K[5,5]        ≈  kₜ
    @test K[6,6]        ≈  kₜ
    @test K[2,5]        ≈ -kₜ
    @test K[5,2]        ≈ -kₜ 
    @test K[3,6]        ≈ -kₜ
    @test K[6,3]        ≈ -kₜ
end
@testset "spurious stiffness" begin
    # no axial force from anything else than displacements about element axis
    @test norm(K[[1,4], [2,3,5,6]])  ≈ 0.
    @test norm(K[[2 5], [1,3,4,6]])  ≈ 0.
    @test norm(K[[3 6], [1,2,4,5]])  ≈ 0.
end


@testset "inertia" begin
    @test M[1,1]        ≈  μ*L₀/3   
    @test M[2,2]        ≈  μ*L₀/3   
    @test M[3,3]        ≈  μ*L₀/3   
    @test M[4,4]        ≈  μ*L₀/3   
    @test M[5,5]        ≈  μ*L₀/3   
    @test M[6,6]        ≈  μ*L₀/3   

    @test M[1,4]        ≈  μ*L₀/6   
    @test M[2,5]        ≈  μ*L₀/6   
    @test M[3,6]        ≈  μ*L₀/6   
    @test M[4,1]        ≈  μ*L₀/6   
    @test M[5,2]        ≈  μ*L₀/6   
    @test M[6,3]        ≈  μ*L₀/6   
end
@testset "spurious inertia" begin
    @test norm(M[1, [2,3,5,6]])  ≈ 0. 
    @test norm(M[4, [2,3,5,6]])  ≈ 0. 
    @test norm(M[2, [1,3,4,6]])  ≈ 0. 
    @test norm(M[5, [1,3,4,6]])  ≈ 0. 
    @test norm(M[3, [1,2,4,5]])  ≈ 0. 
    @test norm(M[6, [1,2,4,5]])  ≈ 0. 
end
;

# Testing weight
w = 10
model           = Model(:TestModel)
node1           = addnode!(model,𝕣[0,0,0])
node2           = addnode!(model,𝕣[L₀,0,0])
elnod           = [model.nod[n.inod] for n∈[node1,node2]]
mat             = AxisymmetricBarCrossSection(EA=EA,μ=μ,w=w)
bar            =  Bar3D(elnod;mat)
x = SVector(0.,     0.,     0.,     0.,    0.,    0.); X = (x,)
R,FB=Muscade.residual(bar,   X,U,A,t,SP,dbg) 
@testset "residual weight" begin
    @test R        ≈  [ 0.0, 0.0, w*L₀/2, 0.0, 0.0, w*L₀/2]
end

## Testing inertia and added mass resultants
Caₜ = 2.
Caₙ = 3.
a1,a2,a3 = 4.0,3.0,2.0;
model           = Model(:TestModel)
node1           = addnode!(model,𝕣[0,0,0])
node2           = addnode!(model,𝕣[L₀,0,0])
elnod           = [model.nod[n.inod] for n∈[node1,node2]]
mat             = AxisymmetricBarCrossSection(EA=EA ,μ=μ, Caₜ=Caₜ, Caₙ=Caₙ)
bar            = Bar3D(elnod;mat)
displacement    =  SVector(0.,0.,0.,  0.,0.,0.); 
velocity        =  SVector(0.,0.,0.,  0.,0.,0.); 

@testset "residual addded mass" begin
    acceleration =  SVector(a1,0.,0.,  a1,0.,0.); 
    X = (displacement,velocity,acceleration); 
    R,FB=Muscade.residual(bar,   X,U,A,t,SP,dbg)     
    @test R        ≈  [ (μ+Caₜ)*a1*L₀/2, 0., 0.,    (μ+Caₜ)*a1*L₀/2, 0., 0.  ] 

    acceleration =  SVector(0.,a2,0.,  0.,a2,0.); 
    X = (displacement,velocity,acceleration); 
    R,FB=Muscade.residual(bar,   X,U,A,t,SP,dbg) 
    @test R        ≈  [ 0., (μ+Caₙ)*a2*L₀/2, 0.,    0., (μ+Caₙ)*a2*L₀/2, 0.] 

    acceleration =  SVector(0.,0.,a3,  0.,0.,a3); 
    X = (displacement,velocity,acceleration); 
    R,FB=Muscade.residual(bar,   X,U,A,t,SP,dbg) 
    @test R        ≈  [ 0., 0., (μ+Caₙ)*a3*L₀/2,    0., 0., (μ+Caₙ)*a3*L₀/2] 

end

## Testing damping
Clₜ = 1.
Clₙ = 2.
Cqₜ = 3.
Cqₙ = 4.
v1,v2,v3 = 1.0,1.1,0.1;
model           = Model(:TestModel)
node1           = addnode!(model,𝕣[0,0,0])
node2           = addnode!(model,𝕣[L₀,0,0])
elnod           = [model.nod[n.inod] for n∈[node1,node2]]
mat             = AxisymmetricBarCrossSection(EA=EA, μ=μ,Clₜ=Clₜ, Clₙ=Clₙ, Cqₜ=Cqₜ, Cqₙ=Cqₙ)
bar             = Bar3D(elnod;mat)
acceleration    =  SVector(0,0.,0.,  0.,0.,0.); 
@testset "residual damping" begin
    displacement    =  SVector(0.,0.,0.,  0.,0.,0.); 
    velocity        =  SVector(v1,0.,0.,  v1,0.,0.); 
    X = (displacement,velocity,acceleration); 
    R,FB=Muscade.residual(bar,   X,U,A,t,SP,dbg) 
    @test R        ≈  [ (Clₜ+Cqₜ*abs(v1))*v1*L₀/2, 0., 0.,    (Clₜ+Cqₜ*abs(v1))*v1*L₀/2, 0., 0.]

    velocity =      SVector(0.,v2,0.,  0.,v2,0.); 
    X = (displacement,velocity,acceleration); 
    R,FB=Muscade.residual(bar,   X,U,A,t,SP,dbg) 
    @test R        ≈  [ 0., (Clₙ+Cqₙ*abs(v2))*v2*L₀/2, 0.,    0., (Clₙ+Cqₙ*abs(v2))*v2*L₀/2, 0.]

    velocity =      SVector(0.,0.,v3,  0.,0.,v3); 
    X = (displacement,velocity,acceleration); 
    R,FB=Muscade.residual(bar,   X,U,A,t,SP,dbg)
    @test R        ≈  [ 0., 0., (Clₙ+Cqₙ*abs(v3))*v3*L₀/2,    0., 0., (Clₙ+Cqₙ*abs(v3))*v3*L₀/2]
end


# # Testing drawing
# # Create model
# nel         = 3
# nodeCoord   = hcat( vcat(0:L₀/nel:L₀,zeros(nel)),zeros(Float64,nel+1 + nel,2))
# mat         = AxisymmetricBarCrossSection(EA=EA,μ=μ)
# model       = Model(:TestModel)
# nodid       = addnode!(model,nodeCoord)
# mesh        = hcat(nodid[1:nel],nodid[2:nel+1],nodid[nel+2:2nel+1])
# eleid       = addelement!(model,Bar3D{true},mesh;mat=mat)

# state    = initialize!(model);
# setdof!(state,  [0.0150628,  0.00279812,  0.0323863,   0.00527616 ],nodID=nodid[1:nel+1     ]         ,field=:t1)
# setdof!(state,  [0.0118553, -0.0150178 ,  0.0118042,   0.000708444],nodID=nodid[1:nel+1     ]         ,field=:t2)
# setdof!(state,  [-0.000945, -0.00591755,  0.0133323,  -0.0545434  ],nodID=nodid[1:nel+1     ]         ,field=:t3)
# setdof!(state,  [0.       ,  0.2        , 0.4                     ],nodID=nodid[nel+2:2nel+1],class=:U,field=:t1)
# setdof!(state,  [0.       ,  0.         , 0.                      ],nodID=nodid[nel+2:2nel+1],class=:U,field=:t2)
# setdof!(state,  [1.       ,  1.3        , 1.6                     ],nodID=nodid[nel+2:2nel+1],class=:U,field=:t3)

# state2 = copy(state)
# state2.X[1] .+= 1.

# using GLMakie
# fig      = Figure(size = (1000,1000))
# ax = Axis3(fig[1,1])
# draw!(ax,state)
# display(fig)

end