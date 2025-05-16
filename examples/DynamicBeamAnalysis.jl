using Revise
using Muscade, StaticArrays, GLMakie
using Printf

include("BeamElements.jl")

# Beam simply supported at both ends  
L = 1;  # Beam length [m]
q = 0.0;  # Uniform lateral load [N/m]
EI₂ = 1;  # Bending stiffness [Nm²]
EI₃ = 1e3;  # Bending stiffness [Nm²]
EA = 1e3;  # Axial stiffness [N]
GJ = 1e3;  # Torsional stiffness [Nm²]
μ = 1;
ι₁= 1;

nel         = 20
Nnod        = nel+1   
nodeCoord   = hcat((0:L/nel:L),zeros(Float64,Nnod,2))
mat         = BeamCrossSection(EA=EA,EI₂=EI₂,EI₃=EI₃,GJ=GJ,μ=μ,ι₁=ι₁)
model       = Model(:TestModel)
nodid       = addnode!(model,nodeCoord)
mesh        = hcat(nodid[1:Nnod-1],nodid[2:Nnod])
eleid       = addelement!(model,EulerBeam3D,mesh;mat=mat,orient2=SVector(0.,1.,0.))
[addelement!(model,Hold,[nodid[1]]  ;field) for field∈[:t1,:t2,:t3,:r1]]                                # Simply supported end 1
[addelement!(model,Hold,[nodid[end]];field) for field∈[:t1,:t2,:t3,:r1]]                                # Simply supported end 2
[addelement!(model,DofLoad,[nodid[nodeidx]];field=:t2,value=t->-min(1,t)*q*L/Nnod) for nodeidx=1:Nnod];          # Distributed vertical load q


# Solve the static analysis problem 
initialstate    = initialize!(model);
state           = solve(SweepX{0};initialstate,time=[0., 1.]);

x_s = getdof(state[2];field=:t1,nodID=nodid[1:Nnod])
y_s = getdof(state[2];field=:t2,nodID=nodid[1:Nnod])
z_s = getdof(state[2];field=:t3,nodID=nodid[1:Nnod])

# Solve eigenvalue problem
nmod            = 20
res             = solve(EigX{ℝ};state=state[2],nmod)

# Compare with analytical solutions for the natural frequency of a simply supported beam 
# See e.g. https://roymech.org/Useful_Tables/Vibrations/Natural_Vibrations_derivation.html
# fₙ(k) = sqrt(EI/m)*(k^2*π)/(2*L^2)

fₙ(k) = sqrt(EI₂/μ)*(k^2*π)/(2*L^2)
@show [fₙ(i)            for i∈1:nmod];
@show [res.ω[i]/(2π)    for i∈1:nmod]; 

# Display a chosen eigenvector
imod            = [2,    2]
A               = [1,    0] 
eigres           = increment(state[2],res,imod,A);
x_eig = getdof(eigres;field=:t1,nodID=nodid[1:Nnod])
y_eig = getdof(eigres;field=:t2,nodID=nodid[1:Nnod])
z_eig = getdof(eigres;field=:t3,nodID=nodid[1:Nnod])

fig      = Figure(size = (2000,1000))
ax = Axis(fig[1,1], ylabel="Modeshape x [m]",        yminorgridvisible = true,xminorgridvisible = true,xticks = (0:L/nel:L))
scatter!(fig[1,1],(0:L/nel:L),  x_[:],                  label="Modeshape");
ax=Axis(fig[2,1], ylabel="Modeshape y [m]",        yminorgridvisible = true,xminorgridvisible = true,xticks = (0:L/nel:L))
scatter!(fig[2,1],(0:L/nel:L),  y_[:],           label="Modeshape");
ax=Axis(fig[3,1], ylabel="Modeshape z [m]",       yminorgridvisible = true,xminorgridvisible = true,xticks = (0:L/nel:L))
scatter!(fig[3,1],(0:L/nel:L),  z_[:],          label="Modeshape");


ax = Axis3(fig[1:3,2],aspect=:equal)
scatter!(ax,nodeCoord[:,1],                 nodeCoord[:,2],                 nodeCoord[:,3],                 label="As meshed");
scatter!(ax,nodeCoord[:,1]+x_s[:],          nodeCoord[:,2]+y_s[:],          nodeCoord[:,3]+z_s[:],          label="Static equilibrium");
scatter!(ax,nodeCoord[:,1]+x_s[:]+x_eig[:], nodeCoord[:,2]+y_s[:]+y_eig[:], nodeCoord[:,3]+z_s[:]+ z_eig[:],label="Modeshape");
# xlims!(ax, 0,1); ylims!(ax, 0,1); zlims!(ax, 0,1); 
axislegend()
display(fig)

