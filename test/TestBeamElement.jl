module TestEulerBeam3D

using Test, Muscade, StaticArrays, LinearAlgebra

a = SA[1,0,0]
b = SA[0,1,1]
r = Muscade.adjust(a,b)
R = Muscade.Rodrigues(r)
u = R*a

v1      = variate{1,3}(SA[.1,.2,.3])
M1      = Muscade.Rodrigues(v1)
w1,w∂v1 = value_∂{1,3}(Muscade.Rodrigues⁻¹(M1))

v2      = variate{1,3}(SA[1e-7,2e-7,1e-8])
M2      = Muscade.Rodrigues(v2)
w2,w∂v2 = value_∂{1,3}(Muscade.Rodrigues⁻¹(M2))

v3      = variate{1,3}(SA[1e-7,2e-7,1e-8])
M3      = Muscade.Rodrigues(v3)
w3,w∂v3 = value_∂{1,3}(Muscade.Rodrigues⁻¹(M3))


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
mat             = Muscade.BeamCrossSection(EA=10.,EI=3.,GJ=4.)

beam            = Muscade.EulerBeam3D(elnod;mat,orient2=SVector(0.,1.,0.))

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
mat             = Muscade.BeamCrossSection(EA=10.,EI=3.,GJ=4.)
beam            = Muscade.EulerBeam3D(elnod;mat,orient2=SVector(0.,1.,0.))
t,SP,dbg  = 0.,(;),(status=:testing,)
U = (SVector{0,𝕣}(),)
A = SVector{0,𝕣}()

x = SVector(0.,0.,0.,0.,0.,0.,0.1,0.0,0.,0.,0.,0.); X = (x,)
R,FB=residual(beam,   X,U,A,t,SP,dbg) 
@testset "residual tension" begin
    @test R        ≈  [-0.9999999999999998, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9999999999999998, 0.0, 0.0, 0.0, 0.0, 0.0]
    @test FB === nothing
end

x = SVector(0.,0.,0.,0.,0.,0.,0.,0.1,0.,0.,0.,0.); X = (x,)
R,FB=residual(beam,   X,U,A,t,SP,dbg) 
@testset "residual flex" begin
    @test R        ≈  [0.305626505038752, -3.557508839178609, 0.0, 0.0, 0.0, -1.7940357448412423, -0.305626505038752, 3.557508839178609, 0.0, 0.0, 0.0, -1.7940357448412423]
    @test FB === nothing
end

x = SVector(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.1); X = (x,)
R,FB=residual(beam,   X,U,A,t,SP,dbg) 
@testset "residual flex" begin
    @test R        ≈  [0.0, 1.8000000000002365, 0.0, 0.0, 0.0, 0.6000000000000143, 0.0, -1.8000000000002365, 0.0, 0.0, 0.0, 1.2000000000002227]
    @test FB === nothing
end

x = SVector(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.1,0.,0.); X = (x,)
R,FB=residual(beam,   X,U,A,t,SP,dbg) 
@testset "residual torsion" begin
    @test R        ≈ [0.0, 0.0, 0.0, -0.40000000000000124, 0.0, 0.0, 0.0, 0.0, 0.0, 0.40000000000000135, 0.0, 0.0]
    @test FB === nothing
end


end

# using Muscade, StaticArrays, GLMakie
# using Muscade.Elements: EulerBeam3D, BeamCrossSection

# L = 10.0;  # Beam length [m]
# q = 10.0;  # Uniform lateral load [N/m]
# EI = 1e6;  # Bending stiffness [Nm²]
# EA = 1e6;  # Axial stiffness [N]
# GJ = 1e3;  # Torsional stiffness [Nm²]
# Fh = 1e5;  # Axial load
# Fv = 1e2;  # Lateral load

# # Beam clamped at both ends, subjected to uniform distributed load of intensity q : 
# # Analytical solutions from https://mechanics.tamu.edu/wp-content/uploads/2017/03/Appendix-A_Exact-Analytical-Solutions.pdf (contains errors)
# # Checked against https://faculty.arch.tamu.edu/anichols/Courses/Arch%20331/Spring%20last/Notes/Files/24_NS8-2beamdiagrams_RycFHMO.pdf 

# x = (0:L/100:L)
# # Deflection (the two sources agree)
# w = q*L^2*x.^2 .* (1.0 .-x/L).^2 ./ (24.0*EI) 
# # Slope (verified by differentiating the above)
# θ = -q*L^2*x.*(1.0 .- 3.0*x/L + 2.0*x.^2/L^2) / (12.0*EI)
# # Curvature (derived)
# κ = -q*L^2*(1.0 .- 6.0*x/L + 6.0*x.^2/L^2) / (12.0*EI)
# # Bending moment (the two sources do not agree)
# M = -q*L^2*(1.0 .- 6.0*x/L + 6.0*x.^2/L^2) / 12.0
# # Shear force (the two sources agree)
# V = q*L*(1.0 .- 2.0*x/L) / 2.0 

# nel         = 20
# nnod        = nel+1   
# nodeCoord   = hcat((0:L/nel:L),zeros(Float64,nnod,2))
# mat         = Muscade.BeamCrossSection(EA=EA,EI=EI,GJ=GJ)
# model       = Model(:TestModel)
# nodid       = addnode!(model,nodeCoord)
# mesh        = hcat(nodid[1:nnod-1],nodid[2:nnod])
# eleid       = addelement!(model,EulerBeam3D,mesh;mat=mat,orient2=SVector(0.,1.,0.))

# [addelement!(model,Hold,[nodid[1]]  ;field) for field∈[:t1,:t2,:t3,:r1,:r2,:r3]]                                # Clamp end 1
# [addelement!(model,Hold,[nodid[end]];field) for field∈[:t1,:t2,:t3,:r1,:r2,:r3]]                                # Clamp end 2
# [addelement!(model,DofLoad,[nodid[nodeidx]];field=:t2,value=t->-min(1,t)*q*L/nnod) for nodeidx=1:nnod]          # Distributed vertical load q
# # addelement!(model,DofLoad,[nodid[end]];field=:t2,value=t->min(1,t)*Fv)                                        # Vertical force Fv on last node
# # addelement!(model,DofLoad,[nodid[end]];field=:t1,value=t->min(1,t)*Fh)                                          # Horizontal force Fh on last node


# initialstate    = initialize!(model)
# state           = solve(SweepX{0};initialstate,time=[0.,1.])

# w_,dofID = getdof(state[2],field=:t2,nodID=nodid[1:nnod])
# θ_,dofID = getdof(state[2],field=:r3,nodID=nodid[1:nnod])
# req = @request gp(resultants(m))
# out = getresult(state[2],req,eleid)
# Mgp1_ = [ out[idxEl].gp[1][:resultants][:m][2] for idxEl ∈ 1:nel]
# Mgp2_ = [ out[idxEl].gp[2][:resultants][:m][2] for idxEl ∈ 1:nel]
# xgp1 = (L/nel)*((0.5-1.0/(2*sqrt(3))):1:nel)
# xgp2 = (L/nel)*((0.5+1.0/(2*sqrt(3))):1:nel)
# req = @request gp(κ)
# out = getresult(state[2],req,eleid)
# κgp1_ = [ out[idxEl].gp[1][:κ][2] for idxEl ∈ 1:nel]
# κgp2_ = [ out[idxEl].gp[2][:κ][2] for idxEl ∈ 1:nel]

# fig      = Figure(size = (1000,1000))
# display(fig) # open interactive window (gets closed down by "save")
# ax=Axis(fig[1,1], ylabel="Deflection w [m]",        yminorgridvisible = true,xminorgridvisible = true,xticks = (0:L/nel:L))
# lines!(fig[1,1],x,            -w,                       label="Analytical solution")
# scatter!(fig[1,1],(0:L/nel:L),  w_[:],                  label="Muscade/beam");
# axislegend()
# ax=Axis(fig[2,1], ylabel="Rotation θ [deg]",        yminorgridvisible = true,xminorgridvisible = true,xticks = (0:L/nel:L))
# lines!(fig[2,1],x,            θ*180/pi,                 label="Analytical solution")
# scatter!(fig[2,1],(0:L/nel:L),  θ_[:]*180/pi,           label="Muscade/beam");
# ax=Axis(fig[3,1], ylabel="Curvature κ [m⁻¹]",       yminorgridvisible = true,xminorgridvisible = true,xticks = (0:L/nel:L))
# lines!(fig[3,1],x,            κ,                        label="Analytical solution")
# scatter!(fig[3,1],[xgp1;xgp2],  [κgp1_;κgp2_],          label="Muscade/beam");
# ax=Axis(fig[4,1], ylabel="Bending moment M [Nm]",   yminorgridvisible = true,xminorgridvisible = true,xticks = (0:L/nel:L), xlabel="Position on beam [m]")
# lines!(fig[4,1],x,            M,                        label="Analytical solution")
# scatter!(fig[4,1],[xgp1;xgp2],  [Mgp1_;Mgp2_],          label="Muscade/beam");
# # ax=Axis(fig[4,1])
# display(fig) # open interactive window (gets closed down by "save")
