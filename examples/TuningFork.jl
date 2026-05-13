using Muscade, StaticArrays, GLMakie, Muscade.Toolbox, Interpolations

# Dimensions of the tuning fork and number of elements
l₀=3e-2; n₀=6       # base - length and number of elements 
l₁=1e-2; n₁=2       # fork - length and number of elements per side
l₂=1e-1; n₂=20      # prong - length and number of elements per side
h = 5e-3;           # width of the tuning tuning fork
b = 1e-3;           # thickness 
E = 210e9;          # Young modulus
G = 79.3e9;         # shear modulus
μ = 7850. *h*b;     # mass per unit length, assumes steel
beamCrossSection = BeamCrossSection(EA=E*h*b, EI₂=E*h^3*b/12, EI₃=E*b*h^3/12, GJ=G*h*b^3/3, μ=μ, ι₁=μ/12*(h^2+b^2))

# Define parametrized point mass element (inducing an inertia force)
struct LumpedMass <: AbstractElement; m :: 𝕣; end
LumpedMass(nod::Vector{Node};m::𝕣) = LumpedMass(m)
@espy function Muscade.residual(o::LumpedMass, X,U,A, t,SP,dbg); return SVector((o.m+A[1]).*∂2(X)),noFB; end
Muscade.doflist( ::Type{LumpedMass})  = (inod =(1,1,1,2), class=(:X,:X,:X,:A), field=(:t1,:t2,:t3,:mass))

# External excitation load
pling = linear_interpolation([-10, 0, 1e-2, 2e-2, 10],[0.,0,20.,0,0]); 
@functor with(pling) pling_(t) = pling(t)

# Time steps
Δt = 1e-4;
nDynamicLoadSteps = 1000;
timeVec = Δt:Δt:nDynamicLoadSteps*Δt

# Meshing: listing nodes
nodeCoord = vcat(
    hcat((0:n₀) * l₀/n₀     , zeros(n₀+1,1) , zeros(n₀+1,1) ), # base nodes
    hcat(ones(n₁,1)*l₀      , (1:n₁)*l₁/n₁  , zeros(n₁,1)   ), # upper fork
    hcat(ones(n₁,1)*l₀      , -(1:n₁)*l₁/n₁ , zeros(n₁,1)   ), # lower fork
    hcat(l₀ .+ (1:n₂)*l₂/n₂ , ones(n₂,1)*l₁ , zeros(n₂,1)   ), # upper prong 
    hcat(l₀ .+ (1:n₂)*l₂/n₂ , -ones(n₂,1)*l₁, zeros(n₂,1)   )  # lower prong 
)

model       = Model(:TuningFork);

# Add beam elements
nodid = addnode!(model,nodeCoord)
mesh = vcat(
    hcat(nodid[1:n₀],nodid[2:n₀+1]),                                                                                # base 
    hcat(nodid[n₀+1:n₀+n₁],nodid[n₀+2:n₀+n₁+1]),                                                                    # upper fork 
    [nodid[n₀+1] nodid[n₀+n₁+2]],         hcat(nodid[n₀+n₁+2:n₀+2n₁],nodid[n₀+n₁+3:n₀+2n₁+1])           ,           # lower fork
    [nodid[n₀+n₁+1] nodid[n₀+2n₁+2]],     hcat(nodid[n₀+2n₁+2:n₀+2n₁+n₂],nodid[n₀+2n₁+3:n₀+2n₁+n₂+1])   ,           # upper prong
    [nodid[n₀+2n₁+1] nodid[n₀+2n₁+n₂+2]], hcat(nodid[n₀+2n₁+n₂+2:n₀+2n₁+2n₂],nodid[n₀+2n₁+n₂+3:n₀+2n₁+2n₂+1])       # lower prong
)
eleid = addelement!(model,EulerBeam3D, mesh;    mat=beamCrossSection,orient2=SVector(0.,0,1))
[addelement!(model,Hold,[nodid[1]]  ;field)  for field∈[:t1,:t2,:t3,:r1,:r2,:r3]]; 

# Add the node that will contain the parameter to optimize
Anod = addnode!(model,[-1., 0, 0])

# Create variations of the model
XAmodel =   deepcopy(model)
XUAmodel =  deepcopy(model)

# Add parasitic lump mass on both XA and XUA models
spuriousMass = 3e-3;
spuriousMassLocation = n₀+2n₁+n₂
addelement!(XAmodel, LumpedMass,  [nodid[spuriousMassLocation] Anod]; m=spuriousMass)
addelement!(XUAmodel, LumpedMass, [nodid[spuriousMassLocation] Anod]; m=spuriousMass)

# Add external load to X and XA models
addelement!(model,  DofLoad,[nodid[n₀+2n₁+1]];field=:t2,value=pling_)
addelement!(XAmodel,DofLoad,[nodid[n₀+2n₁+1]];field=:t2,value=pling_)

# Establish target response
initialState    = initialize!(model;time=0.);
staticState     = solve(SweepX{0};initialstate=initialState, time=[0.])
dynamicStates   = solve(SweepX{2};initialstate=staticState[end], time=timeVec)
vibTarget       = getdof(dynamicStates;field=:t2,nodID=[nodid[n₀+2n₁+n₂+1]])
target          = Muscade.FunctionFromVector(timeVec,vibTarget)

# Cost function deviation to target 
@functor with(σᵥ=1e-4,target) Xcost(x,t)=((x-target(t))/σᵥ)^2

# Add costs on the coorecting mass, XA and XUA model
@functor with(σₘ=1e-3)          Acost(a) = 0.5*(a/σₘ)^2
addelement!(XAmodel,  SingleAcost,[Anod]; field=:mass, cost=Acost)
addelement!(XUAmodel, SingleAcost,[Anod]; field=:mass, cost=Acost)

# Add costs on the deviation to target measurements, XA and XUA model
addelement!(XAmodel, SingleDofCost,[nodid[n₀+2n₁+n₂+1]]; class=:X, field=:t2, cost=Xcost)
addelement!(XUAmodel,SingleDofCost,[nodid[n₀+2n₁+n₂+1]]; class=:X, field=:t2, cost=Xcost)

# Add cost to unknown loads XUA model
@functor with(σᵤ=1e-1) Ucost(u,t) = 0.5*(u/σᵤ)^2*exp(t)
addelement!(XUAmodel,SingleDofCost,[nodid[n₀+2n₁+n₂+1]]; class=:U,field=:t2           ,    cost=Ucost)

# Run analyses

# XA model (before optimization)
XAinitialState    = initialize!(XAmodel;time=0.);
XAdynamicStates   = solve(SweepX{2};initialstate=XAinitialState, time=timeVec)
# XA model (after optimization)
# optimXAstate  = solve(SweepXA{2,2,2};  initialstate=XAinitialState,time=timeVec, verbose=false,catcherror=true,maxAiter=50,maxΔa=1e-10)

# Contents of the three models
Muscade.describe(model,:eletyp)
Muscade.describe(XAmodel,:eletyp)
Muscade.describe(XUAmodel,:eletyp)

# Gather and compare results (add tuned XA results when ready)
analysisResults = (dynamicStates, XAdynamicStates); 
t1 = Matrix{Float64}(undef,length(analysisResults),length(timeVec))
t2 = Matrix{Float64}(undef,length(analysisResults),length(timeVec))
t3 = Matrix{Float64}(undef,length(analysisResults),length(timeVec))
transmittedForce = Matrix{Float64}(undef,length(analysisResults),length(timeVec))
req = @request gp(resultants(fᵢ))
for idx ∈ 1:length(analysisResults)
t1[idx,:] = getdof(analysisResults[idx];field=:t1,nodID=[nodid[n₀+2n₁+n₂+1]])[:]
t2[idx,:] = getdof(analysisResults[idx];field=:t2,nodID=[nodid[n₀+2n₁+n₂+1]])[:]
t3[idx,:] = getdof(analysisResults[idx];field=:t3,nodID=[nodid[n₀+2n₁+n₂+1]])[:]
out = getresult(analysisResults[idx],req,[eleid[1]]);
transmittedForce[idx,:] = [ out[idxEl].gp[1][:resultants][:fᵢ] for idxEl ∈ 1:size(out,2)];
 end

# Plot
fig2      = Figure(size = (1000,1000))
ax5 = Axis(fig2[1,1],ylabel="Excitation [N]")
lines!(ax5,timeVec,pling(timeVec)   )
ax1 = Axis(fig2[2,1],ylabel="xₑ [mm]")
lines!(ax1,timeVec,t1[1,:]*1e3, label="target")
lines!(ax1,timeVec,t1[2,:]*1e3, label="detuned XA")
axislegend(ax1)
ax2 = Axis(fig2[3,1],ylabel="yₑ [mm]")
lines!(ax2,timeVec,t2[1,:]*1e3, label="target")
lines!(ax2,timeVec,t2[2,:]*1e3, label="detuned XA")
# ax3 = Axis(fig2[3,1],ylabel="zₑ [mm]")
# lines!(ax3,timeVec,t3[1,:]*1e3, label="target")
# lines!(ax3,timeVec,t3[2,:]*1e3, label="detuned XA")
# ax4 = Axis(fig2[4,1],ylabel=L"F_{trans} [N]")
# lines!(ax4,timeVec,transmittedForce[1,:], label="target")
# lines!(ax4,timeVec,transmittedForce[2,:], label="detuned XA")
display(fig2)


# DirectXUA
# XUAinitialState    = initialize!(XUAmodel;time=0.);
# XUAdynamicStates   = solve(DirectXUA{2,0,1};initialstate=[XUAinitialState], time=[timeVec])


# More plotting for the future
# α      = 2π*(0:19)/20; circle = .005*[cos.(α) sin.(α)]'
# fig      = Figure(size = (1000,1000))
# ax = Axis3(fig[1,1],xgridvisible=false,ygridvisible=false,zgridvisible=false,aspect=:data)
# draw!(ax,initialState,EulerBeam3D=(;style=:shape))
# draw!(ax,staticState[end],EulerBeam3D=(;style=:shape))
# display(fig);
# Produce an animation
# # α      = 2π*(0:19)/20; circle = .1*[cos.(α) sin.(α)]'
# fig2      = Figure(size = (1000,1000))
# ax2 = Axis3(fig2[1,1],xgridvisible=false,ygridvisible=false,zgridvisible=false)
# xlims!(ax2,0,l₀+l₂); ylims!(ax2,-l₁,l₁);# zlims!(ax2,-,10)
# graphic = draw!(ax2,dynamicStates[1],EulerBeam3D=(;style=:shape))
# # ax2.azimuth[]=-π/2+π/180*10;
# # ax2.elevation[]=0+π/180*10;
# framerate = 20
# loadStepsIterator = 1:3:nDynamicLoadSteps
# record(fig2, "tuningFork.mp4", loadStepsIterator;
#         framerate = framerate) do stateIdx
#            draw!(graphic,dynamicStates[stateIdx])
# end