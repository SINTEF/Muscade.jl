# # Static analysis of a beam
#
using Muscade, StaticArrays, GLMakie

include("BeamElements.jl")

# # Analysis

# Defining the beam properties and loading. 
# Beam clamped at both ends, subjected to uniform distributed load of intensity q  

L = 10.0;  # Beam length [m]
q = 10.0;  # Uniform lateral load [N/m]
EI = 1e6;  # Bending stiffness [Nm²]
EA = 1e6;  # Axial stiffness [N]
GJ = 1e3;  # Torsional stiffness [Nm²]
#src # Fh = 1e5;  # Axial load
#src # Fv = 1e2;  # Lateral load

# Analytical solutions from [here](https://mechanics.tamu.edu/wp-content/uploads/2017/03/Appendix-A_Exact-Analytical-Solutions.pdf) (which contains errors), and checked against [this source](https://faculty.arch.tamu.edu/anichols/Courses/Arch%20331/Spring%20last/Notes/Files/24_NS8-2beamdiagrams_RycFHMO.pdf)
x = (0:L/100:L);
# Deflection (the two sources agree)
w = q*L^2*x.^2 .* (1.0 .-x/L).^2 ./ (24.0*EI); 
# Slope (verified by differentiating the above)
θ = -q*L^2*x.*(1.0 .- 3.0*x/L + 2.0*x.^2/L^2) / (12.0*EI); 
# Curvature (derived)
κ = -q*L^2*(1.0 .- 6.0*x/L + 6.0*x.^2/L^2) / (12.0*EI); 
# Bending moment (the two sources do not agree)
M = -q*L^2*(1.0 .- 6.0*x/L + 6.0*x.^2/L^2) / 12.0; 
# Shear force (the two sources agree)
V = q*L*(1.0 .- 2.0*x/L) / 2.0; 

# Create the model 
nel         = 20
nnodes        = nel+1   
nodeCoord   = hcat((0:L/nel:L),zeros(Float64,nnodes,2))
mat         = BeamCrossSection(EA=EA,EI=EI,GJ=GJ)
model       = Model(:TestModel)
nodid       = addnode!(model,nodeCoord)
mesh        = hcat(nodid[1:nnodes-1],nodid[2:nnodes])
eleid       = addelement!(model,EulerBeam3D,mesh;mat=mat,orient2=SVector(0.,1.,0.))

[addelement!(model,Hold,[nodid[1]]  ;field) for field∈[:t1,:t2,:t3,:r1,:r2,:r3]]                                # Clamp end 1
[addelement!(model,Hold,[nodid[end]];field) for field∈[:t1,:t2,:t3,:r1,:r2,:r3]]                                # Clamp end 2
[addelement!(model,DofLoad,[nodid[nodeidx]];field=:t2,value=t->-min(1,t)*q*L/nnodes) for nodeidx=1:nnodes];          # Distributed vertical load q
#src # addelement!(model,DofLoad,[nodid[end]];field=:t2,value=t->min(1,t)*Fv)                                        # Vertical force Fv on last node
#src # addelement!(model,DofLoad,[nodid[end]];field=:t1,value=t->min(1,t)*Fh)                                          # Horizontal force Fh on last node

# Solve the problem 
initialstate    = initialize!(model);
state           = solve(SweepX{0};initialstate,time=[0.,1.])

# Fetch results
w_ = getdof(state[2];field=:t2,nodID=nodid[1:nnodes])
θ_ = getdof(state[2];field=:r3,nodID=nodid[1:nnodes])
req = @request gp(resultants(m))
out = getresult(state[2],req,eleid)
Mgp1_ = [ out[idxEl].gp[1][:resultants][:m][2] for idxEl ∈ 1:nel]
Mgp2_ = [ out[idxEl].gp[2][:resultants][:m][2] for idxEl ∈ 1:nel]
xgp1 = (L/nel)*((0.5-1.0/(2*sqrt(3))):1:nel)
xgp2 = (L/nel)*((0.5+1.0/(2*sqrt(3))):1:nel)
req = @request gp(κ)
out = getresult(state[2],req,eleid)
κgp1_ = [ out[idxEl].gp[1][:κ][2] for idxEl ∈ 1:nel]
κgp2_ = [ out[idxEl].gp[2][:κ][2] for idxEl ∈ 1:nel];

# Display results
fig      = Figure(size = (1000,1000))
ax = Axis(fig[1,1], ylabel="Deflection w [m]",        yminorgridvisible = true,xminorgridvisible = true,xticks = (0:L/nel:L))
lines!(fig[1,1], x,            -w,                       label="Analytical solution")
scatter!(fig[1,1],(0:L/nel:L),  w_[:],                  label="Muscade/beam");
axislegend()
ax=Axis(fig[2,1], ylabel="Rotation θ [deg]",        yminorgridvisible = true,xminorgridvisible = true,xticks = (0:L/nel:L))
lines!(fig[2,1],x,            θ*180/pi,                 label="Analytical solution")
scatter!(fig[2,1],(0:L/nel:L),  θ_[:]*180/pi,           label="Muscade/beam");
ax=Axis(fig[3,1], ylabel="Curvature κ [m⁻¹]",       yminorgridvisible = true,xminorgridvisible = true,xticks = (0:L/nel:L))
lines!(fig[3,1],x,            κ,                        label="Analytical solution")
scatter!(fig[3,1],[xgp1;xgp2],  [κgp1_;κgp2_],          label="Muscade/beam");
ax=Axis(fig[4,1], ylabel="Bending moment M [Nm]",   yminorgridvisible = true,xminorgridvisible = true,xticks = (0:L/nel:L), xlabel="Position on beam [m]")
lines!(fig[4,1],x,            M,                        label="Analytical solution")
scatter!(fig[4,1],[xgp1;xgp2],  [Mgp1_;Mgp2_],          label="Muscade/beam");

currentDir = @__DIR__
if occursin("build", currentDir)
    save(normpath(joinpath(currentDir,"..","src","assets","beam.png")),fig)
elseif occursin("examples", currentDir)
    save(normpath(joinpath(currentDir,"beam.png")),fig)
end
# ![Result](assets/beam.png)