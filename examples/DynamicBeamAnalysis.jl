# # Modal analysis of a beam
#

# Beam simply supported at both ends  
using Revise
using Muscade, StaticArrays, GLMakie
using Printf
include("BeamElements.jl");

L = 1;  # Beam length [m]
q = 0.0;  # Uniform lateral load [N/m]
EI₂ = 1;  # Bending stiffness [Nm²]
EI₃ = 1;  # Bending stiffness [Nm²]
EA = 1e6;  # Axial stiffness [N]
GJ = 1e6;  # Torsional stiffness [Nm²]
μ = 1;
ι₁= 1;

# Create model
nel         = 50
Nnod        = nel+1   
nodeCoord   = hcat((0:L/nel:L),zeros(Float64,Nnod,2))
mat         = BeamCrossSection(EA=EA,EI₂=EI₂,EI₃=EI₃,GJ=GJ,μ=μ,ι₁=ι₁)
model       = Model(:TestModel)
nodid       = addnode!(model,nodeCoord)
mesh        = hcat(nodid[1:Nnod-1],nodid[2:Nnod])
eleid       = addelement!(model,EulerBeam3D,mesh;mat=mat,orient2=SVector(0.,1.,0.)); 

# Set boundary conditions and constraints 
[addelement!(model,Hold,[nodid[1]]  ;field) for field∈[:t1,:t2,:t3,:r1]]                                # Simply supported end 1
[addelement!(model,Hold,[nodid[end]];field) for field∈[:t1,:t2,:t3,:r1]]                                # Simply supported end 2
[addelement!(model,Hold,[nodid[nodeidx]];field=:t3) for nodeidx∈2:Nnod-1]                               # Enforce beam motions in one dimension to obtain planar modeshapes
[addelement!(model,DofLoad,[nodid[nodeidx]];field=:t2,value=t->sin(t)*q*L/Nnod) for nodeidx=1:Nnod];    # Distributed vertical load q

# Static analysis
initialstate    = initialize!(model);
state           = solve(SweepX{0};initialstate,time=[0.]);

# Solve eigenvalue problem
nmod            = 15
res             = solve(EigX{ℝ};state=state[1],nmod);

# Analytical solutions for the natural frequency of a simply supported beam 
# See e.g. https://roymech.org/Useful_Tables/Vibrations/Natural_Vibrations_derivation.html
fₙ(k) = √(EI₂/μ)*(k^2*π)/(2*L^2)
Φₙ(k,x) = sin.(k*π/L.*x);

# Display solution and comparison against analytical solution
fig      = Figure(size = (2000,1000))
axes = [Axis(fig[idxLine,1], yminorgridvisible = false,xminorgridvisible = false ) for idxLine=1:3]
axes[1].title = "Modeshapes of a simply supported beam. Muscade (" *string(nel)*" elements): markers. Analytical solution: lines. "
for idxMod=1:nmod
    eigres  = increment(state[1],res,[idxMod],[1]);
    t2_eig  = getdof(eigres;field=:t2,nodID=nodid[1:Nnod])
    δ       = sign(Φₙ(idxMod,0:L/nel:L)'*t2_eig) * maximum(t2_eig)
    selectAxis = axes[mod(idxMod-1,3)+1]
    labelStr= "Mode "*string(idxMod)*", Muscade: "*string(round(res.ω[idxMod]/(2π),digits=3))*" Hz, Analytical: " *string(round(fₙ(idxMod),digits=3))* " Hz"
    scatter!(selectAxis,(0:L/nel:L),  t2_eig[:]/δ,          label=labelStr  );
    lines!(  selectAxis,(0:L/nel:L),  Φₙ(idxMod,0:L/nel:L)                  );
end
for ax∈axes; 
    xlims!(ax,0,1); ylims!(ax, -2,2); axislegend(ax)
end

currentDir = @__DIR__
if occursin("build", currentDir)
    save(normpath(joinpath(currentDir,"..","src","assets","beamModes.png")),fig)
elseif occursin("examples", currentDir)
    save(normpath(joinpath(currentDir,"beamModes.png")),fig)
end
# ![Result](assets/beamModes.png)

#src Dynamic analysis
#src T               = 0.01 *(1:1000)
#src dynAnalysis     = solve(SweepX{2};initialstate=state[1],time=T) 
#src ERROR: MethodError: no method matching motion{2}(::Tuple{SVector{12, Float64}, SVector{12, ∂ℝ{1, 1, Float64}}, SVector{12, ∂ℝ{1, 1, Float64}}})
#src Closest candidates are:
#src   (::Type{motion{P}} where P)()
#src    @ Muscade C:\Users\thsa\code\Muscade.jl\src\Taylor.jl:1
#src   motion{P}(::Tuple{Vararg{SVector{N, R}, ND}}) where {ND, P, N, R}
#src    @ Muscade C:\Users\thsa\code\Muscade.jl\src\Taylor.jl:24