# # Static analysis of a beam
#
using Muscade, StaticArrays#, GLMakie
using Printf

include("BeamElements.jl")

# # Analysis

# Defining the beam properties and loading. 
# Beam simply supported at both ends  
L = 20;  # Beam length [m]
EI = 1.;  # Bending stiffness [Nm²]
EA = 1.#e6;  # Axial stiffness [N]
GJ = 1.#e3;  # Torsional stiffness [Nm²]
#m  = 1;  # Mass per unit length [kg/m]
nModes = 5; 

# Improved beam element requires
# 1) Add inertia loads, drag loads, Use `Muscade.motion` to create Adiffs that will facilitate the dynamic computation.
# 2) split EIx, EIy, and add mass per unit length as beam properties
# 2) Add U-dofs, using a "isoparametric" formulation (?)
            # Alt 1: if U is to be global, u has to be rotated in local, and then interpolated (differently for axial and lateral)
            # Alt 2 (prefered): separate force element that corresponds to a force at each node, which will be interpolated and used by beam element in R. Enables defining costs on U. 
# 4) performance.  Liberal use of nested Adiff makes code simple, but not fast... Replace the Jacobian in some way (Vegard)

# Then create a model 
nel         = 20
Nnod        = nel+1   
nodeCoord   = hcat((0:L/nel:L),zeros(Float64,Nnod,2))
mat         = BeamCrossSection(EA=EA,EI=EI,GJ=GJ)
# beam        = EulerBeam3D(elnod;mat,orient2=SVector(0.,1.,0.))
model       = Model(:TestModel)
nodid       = addnode!(model,nodeCoord)
mesh        = hcat(nodid[1:Nnod-1],nodid[2:Nnod])
eleid       = addelement!(model,EulerBeam3D,mesh;mat=mat,orient2=SVector(0.,1.,0.))
[addelement!(model,Hold,[nodid[1]]  ;field) for field∈[:t1,:t2,:t3]]                                # Simply supported end 1
[addelement!(model,Hold,[nodid[end]];field) for field∈[:t1,:t2,:t3]]                                # Simply supported end 2
# [addelement!(model,DofLoad,[nodid[nodeidx]];field=:t2,value=t->-min(1,t)*q*L/Nnod) for nodeidx=1:Nnod];          # Distributed vertical load q


t= 0.
SP = (;)
dbg  = (status=:testing,)
x = SVector(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.)
X = (x,x,x)
U = (SVector{0,Float64}(),)
A = SVector{0,Float64}()
R,FB = Muscade.residual(model.eleobj[1][1],   X,U,A,t,SP,dbg)

out2 = diffed_residual(model.eleobj[1][1]; X,U,A)
@printf "\nR\n"
print_element_array(model.eleobj[1][1],:X,out2.R)
@printf "\n∂R/∂X₀\n"
print_element_array(model.eleobj[1][1],:X,out2.∇R[2][3])
# @test out2.R   ≈ [-1.]
# @test out2.∇R[2][1] ≈ [0.;;]


# # Solve the static analysis problem 
# initialstate    = initialize!(model);
# state           = solve(SweepX{0};initialstate,time=[0.,1.])


# # Fetch results
# w_ = getdof(state[2];field=:t2,nodID=nodid[1:Nnod])
# θ_ = getdof(state[2];field=:r3,nodID=nodid[1:Nnod])
# req = @request gp(resultants(m))
# out = getresult(state[2],req,eleid)
# Mgp1_ = [ out[idxEl].gp[1][:resultants][:m][2] for idxEl ∈ 1:nel]
# Mgp2_ = [ out[idxEl].gp[2][:resultants][:m][2] for idxEl ∈ 1:nel]
# xgp1 = (L/nel)*((0.5-1.0/(2*sqrt(3))):1:nel)
# xgp2 = (L/nel)*((0.5+1.0/(2*sqrt(3))):1:nel)
# req = @request gp(κ)
# out = getresult(state[2],req,eleid)
# κgp1_ = [ out[idxEl].gp[1][:κ][2] for idxEl ∈ 1:nel]
# κgp2_ = [ out[idxEl].gp[2][:κ][2] for idxEl ∈ 1:nel];

# # Display results
# fig      = Figure(size = (1000,1000))
# ax = Axis(fig[1,1], ylabel="Deflection w [m]",        yminorgridvisible = true,xminorgridvisible = true,xticks = (0:L/nel:L))
# lines!(fig[1,1], x,            -w,                       label="Analytical solution")
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

# Extract mass and stiffness matrix (espy)

# Solve eigenvalue problem
 
# # Compare with analytical solutions for the natural frequency of a simply supported beam 
# # See e.g. https://roymech.org/Useful_Tables/Vibrations/Natural_Vibrations_derivation.html
# fₙ(k) = sqrt(EI/m)*(k^2*π)/(2*L^2)


# fig      = Figure(size = (1000,1000))
# ax = Axis(fig[1,1], ylabel="Natural frequency [Hz]", yminorgridvisible = true,xminorgridvisible = true,xticks = (1:nModes))
# scatter!(fig[1,1], collect(1:nModes), [fₙ(k) for k∈1:nModes], label="Analytical solution")
# axislegend()
# currentDir = @__DIR__
# if occursin("build", currentDir)
#     save(normpath(joinpath(currentDir,"..","src","assets","beamModes.png")),fig)
# elseif occursin("examples", currentDir)
#     save(normpath(joinpath(currentDir,"beamModes.png")),fig)
# end
