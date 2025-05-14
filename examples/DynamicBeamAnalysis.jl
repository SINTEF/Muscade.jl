using Revise
using Muscade, StaticArrays, GLMakie
using Printf

include("BeamElements.jl")

# Beam simply supported at both ends  
L = 1;  # Beam length [m]
q = 0.0;  # Uniform lateral load [N/m]
EI = 1;  # Bending stiffness [Nm²]
EA = 1;  # Axial stiffness [N]
GJ = 1;  # Torsional stiffness [Nm²]
μ = 1


nel         = 5
Nnod        = nel+1   
nodeCoord   = hcat((0:L/nel:L),zeros(Float64,Nnod,2))
mat         = BeamCrossSection(EA=EA,EI=EI,GJ=GJ)
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

# Solve eigenvalue problem
nmod            = 2
res             = solve(EigX;state=state[2],nmod)

# Display eigenvectors
imod            = [1,    2]
A               = [1,    0] 
eigres           = increment(state[2],res,imod,A);

# Analytical solution
fₙ(k) = sqrt(EI/μ)*(k^2*π)/(2*L^2)
@show [fₙ(i)            for i∈1:nmod];
@show [res.ω[i]/(2π)    for i∈1:nmod]; 


x_ = getdof(eigres;field=:t1,nodID=nodid[1:Nnod])
y_ = getdof(eigres;field=:t2,nodID=nodid[1:Nnod])
z_ = getdof(eigres;field=:t3,nodID=nodid[1:Nnod])




# t= 0.
# SP = (;)
# dbg  = (status=:testing,)
# x = SVector(0.,0.,0.,0.,0.,0.,
#             0.,0.,0.,0.,0.,0.)
# X = (x,x,x)
# U = (SVector{0,Float64}(),)
# A = SVector{0,Float64}()
# # R,FB = Muscade.residual(model.eleobj[1][1],   X,U,A,t,SP,dbg)
# out2 = diffed_residual(model.eleobj[1][1]; X,U,A)
# @printf "\nR\n"
# print_element_array(model.eleobj[1][1],:X,out2.R)
# @printf "\n∂R/∂X₀\n"
# print_element_array(model.eleobj[1][1],:X,out2.∇R[2][1])
# @printf "\n∂R/∂X₂\n"
# print_element_array(model.eleobj[1][1],:X,out2.∇R[2][3])


# # Fetch and display results
# nnodes = Nnod
# w_ = getdof(state[2];field=:t2,nodID=nodid[1:nnodes])
# θ_ = getdof(state[2];field=:r3,nodID=nodid[1:nnodes])
# req = @request gp(resultants(mᵢ))
# out = getresult(state[2],req,eleid)
# Mgp1_ = [ out[idxEl].gp[1][:resultants][:mᵢ][2] for idxEl ∈ 1:nel]
# Mgp2_ = [ out[idxEl].gp[2][:resultants][:mᵢ][2] for idxEl ∈ 1:nel]
# Mgp3_ = [ out[idxEl].gp[3][:resultants][:mᵢ][2] for idxEl ∈ 1:nel]
# Mgp4_ = [ out[idxEl].gp[4][:resultants][:mᵢ][2] for idxEl ∈ 1:nel]
# xgp1 = (L/nel)*( (0.5-1/2*sqrt(3/7+2/7*sqrt(6/5))) :1:nel)
# xgp2 = (L/nel)*( (0.5-1/2*sqrt(3/7-2/7*sqrt(6/5))) :1:nel)
# xgp3 = (L/nel)*( (0.5+1/2*sqrt(3/7-2/7*sqrt(6/5))) :1:nel)
# xgp4 = (L/nel)*( (0.5+1/2*sqrt(3/7+2/7*sqrt(6/5))) :1:nel)
# req = @request gp(κ)
# out = getresult(state[2],req,eleid)
# κgp1_ = [ out[idxEl].gp[1].κ[1][2] for idxEl ∈ 1:nel]
# κgp2_ = [ out[idxEl].gp[2].κ[1][2] for idxEl ∈ 1:nel];
# κgp3_ = [ out[idxEl].gp[3].κ[1][2] for idxEl ∈ 1:nel];
# κgp4_ = [ out[idxEl].gp[4].κ[1][2] for idxEl ∈ 1:nel];

# fig      = Figure(size = (1000,1000))
# ax = Axis(fig[1,1], ylabel="Deflection w [m]",        yminorgridvisible = true,xminorgridvisible = true,xticks = (0:L/nel:L))
# scatter!(fig[1,1],(0:L/nel:L),  w_[:],                  label="Muscade/beam");
# axislegend()
# ax=Axis(fig[2,1], ylabel="Rotation θ [deg]",        yminorgridvisible = true,xminorgridvisible = true,xticks = (0:L/nel:L))
# scatter!(fig[2,1],(0:L/nel:L),  θ_[:]*180/pi,           label="Muscade/beam");
# ax=Axis(fig[3,1], ylabel="Curvature κ [m⁻¹]",       yminorgridvisible = true,xminorgridvisible = true,xticks = (0:L/nel:L))
# scatter!(fig[3,1],[xgp1;xgp2;xgp3;xgp4],  [κgp1_;κgp2_;κgp3_;κgp4_],          label="Muscade/beam");
# ax=Axis(fig[4,1], ylabel="Bending moment M [Nm]",   yminorgridvisible = true,xminorgridvisible = true,xticks = (0:L/nel:L), xlabel="Position on beam [m]")
# scatter!(fig[4,1],[xgp1;xgp2;xgp3;xgp4],  [Mgp1_;Mgp2_;Mgp3_;Mgp4_],          label="Muscade/beam");
# display(fig)

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
