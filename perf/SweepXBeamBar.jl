# Benchmarking performance of Muscade analyses on a real problems
# Perform the dynamic analysis of a steel catenary riser (SCR), using both beam and bar elements, contact elements, and constraints for top motions.
# Store computational time for future reference, and generate flame graphs


using Muscade, StaticArrays, Muscade.Toolbox, Interpolations, CSV, DataFrames
using BenchmarkTools, Dates
using Profile,ProfileSVG

include("CommonProperties.jl")
include("BenchmarkHistory.jl")

# Create Muscade beam element model
BeamModel   = Model(:CatenaryRiserBeam);
topNode     = addnode!(BeamModel,[0,0,0])
x1_mat      = BeamCrossSection(EA=x1_EA, EI₂=x1_EI₂, EI₃=x1_EI₃, GJ=x1_GJ, μ=x1_μ, ι₁=x1_ι₁, w=x1_w, g̃=ramping_gravityField, Ca₂=x1_Ca₂, Cq₂=x1_Cq₂, Ca₃=x1_Ca₃, Cq₃=x1_Cq₃)
x2_mat      = BeamCrossSection(EA=x2_EA, EI₂=x2_EI₂, EI₃=x2_EI₃, GJ=x2_GJ, μ=x2_μ, ι₁=x2_ι₁, w=x2_w, g̃=ramping_gravityField, Ca₂=x2_Ca₂, Cq₂=x2_Cq₂, Ca₃=x2_Ca₃, Cq₃=x2_Cq₃)
xSection    = [x1_mat,  x2_mat, x1_mat]; 
nodeList, elementList, nodeCoord = MeshLine!(BeamModel, topNode, π, EulerBeam3D, xSection, segLength, nel; orient2=SVector(0.,1,0));
[addelement!(BeamModel,Hold,[nodeList[1][1]]  ;field)  for field∈[:t1,:t2,:t3,:r1]]; 
addelement!(BeamModel,DofConstraint,[nodeList[3][end]]; λclass=:X, xinod=(1,),xfield=(:t1,), λinod=(1,),λfield=(:λt1,), gap=horizMove, mode=Muscade.equal)
addelement!(BeamModel,DofConstraint,[nodeList[3][end]]; λclass=:X, xinod=(1,),xfield=(:t3,), λinod=(1,),λfield=(:λt3,), gap=vertMove , mode=Muscade.equal);
[addelement!(BeamModel,SoilContact,[nodeList[1][idxNod]],z₀=0.,Kh=1.0e3,Kv=1.0e4,Ch=0.,Cv=0.) for idxNod = 1:length(nodeList[1])];


# Create Muscade bar element model
BarModel    = Model(:CatenaryRiserBar);
topNode     = addnode!(BarModel,[0,0,0])
x1_mat      = AxisymmetricBarCrossSection(EA=x1_EA, μ=x1_μ, w=x1_w, g̃=ramping_gravityField, Caₙ=x1_Ca₂, Cqₙ=x1_Cq₂)
x2_mat      = AxisymmetricBarCrossSection(EA=x2_EA, μ=x2_μ, w=x2_w, g̃=ramping_gravityField, Caₙ=x2_Ca₂, Cqₙ=x2_Cq₂)
xSection    = [x1_mat,  x2_mat, x1_mat]; 
nodeList, elementList, nodeCoord = MeshLine!(BarModel, topNode, π, Bar3D, xSection, segLength, nel);
[addelement!(BarModel,Hold,[nodeList[1][1]]  ;field)  for field∈[:t1,:t2,:t3]]; 
addelement!(BarModel,DofConstraint,[nodeList[3][end]]; λclass=:X, xinod=(1,),xfield=(:t1,), λinod=(1,),λfield=(:λt1,), gap=horizMove, mode=Muscade.equal)
addelement!(BarModel,DofConstraint,[nodeList[3][end]]; λclass=:X, xinod=(1,),xfield=(:t3,), λinod=(1,),λfield=(:λt3,), gap=vertMove , mode=Muscade.equal);
[addelement!(BarModel,SoilContact,[nodeList[1][idxNod]],z₀=0.,Kh=1.0e3,Kv=1.0e4,Ch=0.,Cv=0.) for idxNod = 1:length(nodeList[1])];


staticLoadSteps = (-10:.2:0)*1.;    
nStaticLoadSteps= length(staticLoadSteps)
dynamicLoadSteps= (0.1:0.3:100)*1.0

for model ∈ [BarModel BeamModel]
local initialstate = initialize!(model);

println(model.ID)
println("   Static analysis")
local staticStates        = solve(SweepX{0};initialstate,time=staticLoadSteps,verbose=false);

println("   First run of the dynamic analysis to trigger compilation")
solve(SweepX{2}; initialstate=staticStates[nStaticLoadSteps], time=dynamicLoadSteps, verbose=false);

println("   Second run to assess computational time")
local benchmark_seconds   = @belapsed solve(SweepX{2}; initialstate=$staticStates[nStaticLoadSteps], time=$dynamicLoadSteps, verbose=false);
newEntryInBenchmarkHistory(model.ID,benchmark_seconds)

println("   Third run for detailed profiling")
Profile.clear()
@profile solve(SweepX{2}; initialstate=staticStates[nStaticLoadSteps], time=dynamicLoadSteps, verbose=false);
ProfileSVG.save("Profile_$(model.ID).svg",width=1800, maxdepth=100, maxframes=100000)
end