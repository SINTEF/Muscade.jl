using Muscade, StaticArrays, GLMakie, Muscade.Toolbox

l₀=3e-2; n₀=6
l₁=1e-2; n₁=2
l₂=1e-1; n₂=20

h = 5e-3; 
b = 1e-3; 
E = 210e9;
G = 79.3e9; 
μ = 7850. *h*b; 
g = 9.81;
beamCrossSection = BeamCrossSection(EA=E*h*b, EI₂=E*h^3*b/12, EI₃=E*b*h^3/12, GJ=G*h*b^3/3, μ=μ, ι₁=μ/12*(h^2+b^2); w=μ*g, Cl₁=0*1e3, Cl₂=0*1e3, Cl₃=0*1e3)

nodeCoord = vcat(
    hcat((0:n₀) * l₀/n₀     , zeros(n₀+1,1) , zeros(n₀+1,1) ),
    hcat(ones(n₁,1)*l₀      , (1:n₁)*l₁/n₁  , zeros(n₁,1)   ),
    hcat(ones(n₁,1)*l₀      , -(1:n₁)*l₁/n₁ , zeros(n₁,1)   ),
    hcat(l₀ .+ (1:n₂)*l₂/n₂ , ones(n₂,1)*l₁ , zeros(n₂,1)   ),
    hcat(l₀ .+ (1:n₂)*l₂/n₂ , -ones(n₂,1)*l₁, zeros(n₂,1)   )
)



model       = Model(:TuningFork);

nodid = addnode!(model,nodeCoord)

mesh = vcat(
    hcat(nodid[1:n₀],nodid[2:n₀+1]),
    hcat(nodid[n₀+1:n₀+n₁],nodid[n₀+2:n₀+n₁+1]), 
    [nodid[n₀+1] nodid[n₀+n₁+2]],         hcat(nodid[n₀+n₁+2:n₀+2n₁],nodid[n₀+n₁+3:n₀+2n₁+1])           ,
    [nodid[n₀+n₁+1] nodid[n₀+2n₁+2]],     hcat(nodid[n₀+2n₁+2:n₀+2n₁+n₂],nodid[n₀+2n₁+3:n₀+2n₁+n₂+1])   ,
    [nodid[n₀+2n₁+1] nodid[n₀+2n₁+n₂+2]], hcat(nodid[n₀+2n₁+n₂+2:n₀+2n₁+2n₂],nodid[n₀+2n₁+n₂+3:n₀+2n₁+2n₂+1])
)

eleid = addelement!(model,EulerBeam3D, mesh;    mat=beamCrossSection,orient2=SVector(0.,0,1))
[addelement!(model,Hold,[nodid[1]]  ;field)  for field∈[:t1,:t2,:t3,:r1,:r2,:r3]]; 

pling = Muscade.FunctionFromVector([-10, 0, 1e-2, 2e-2, 10],[0,0,10,0,0]); @functor with() pling_(t) = pling(t)
addelement!(model,DofLoad,[nodid[n₀+2n₁+n₂+1]];field=:t2,value=pling_)

XAmodel = deepcopy(model)

Δt = 5e-4;
nDynamicLoadSteps = 100;
timeVec = Δt:Δt:nDynamicLoadSteps*Δt

# X model
initialState    = initialize!(model;time=0.);
dynamicStates   = solve(SweepX{2};initialstate=initialState, time=timeVec)
vibTarget       = getdof(dynamicStates;field=:t2,nodID=[nodid[n₀+2n₁+n₂+1]])
target          = Muscade.FunctionFromVector(timeVec,vibTarget)


# XA model
Anod = addnode!(XAmodel,[-1., 0, 0])
struct LumpedMass <: AbstractElement
    m :: 𝕣
end
LumpedMass(nod::Vector{Node};m::𝕣) = LumpedMass(m)
@espy function Muscade.residual(o::LumpedMass, X,U,A, t,SP,dbg) 
    return SVector((o.m+A[1])*∂2(X)),noFB
end
Muscade.doflist( ::Type{LumpedMass})  = (inod =(1,1,1,2), class=(:X,:X,:X,:A), field=(:t1,:t2,:t3,:mass))
addelement!(XAmodel, LumpedMass, [nodid[n₀+2n₁+n₂] Anod]; m=.01)

@functor with(σₘ=1e-4,target) Acost(a) = 0.5*(a / σₘ)^2
addelement!(XAmodel,SingleAcost  ,[Anod];   field=:mass,           cost=Acost)
@functor with(σᵥ=1e-4,target) xcost(x,t)=((x-target(t))/σᵥ)^2
addelement!(XAmodel,SingleDofCost,[nodid[n₀+2n₁+n₂+1]];class=:X ,field=:t2    ,cost=xcost)


# XA model (initial)
XAinitialState    = initialize!(XAmodel;time=0.);
XAdynamicStates   = solve(SweepX{2};initialstate=XAinitialState, time=timeVec)
# XA model (optimized)
optimXAstate  = solve(SweepXA{2,2,2};  initialstate=XAinitialState,time=timeVec, verbose=true,catcherror=true,maxAiter=50,maxΔa=1e-10)



t1Vib = getdof(dynamicStates;field=:t1,nodID=[nodid[n₀+2n₁+n₂+1],nodid[n₀+2n₁+2n₂+1]])
t2Vib = getdof(dynamicStates;field=:t2,nodID=[nodid[n₀+2n₁+n₂+1],nodid[n₀+2n₁+2n₂+1]])
t3Vib = getdof(dynamicStates;field=:t3,nodID=[nodid[n₀+2n₁+n₂+1],nodid[n₀+2n₁+2n₂+1]])
XAt1Vib = getdof(XAdynamicStates;field=:t1,nodID=[nodid[n₀+2n₁+n₂+1],nodid[n₀+2n₁+2n₂+1]])
XAt2Vib = getdof(XAdynamicStates;field=:t2,nodID=[nodid[n₀+2n₁+n₂+1],nodid[n₀+2n₁+2n₂+1]])
XAt3Vib = getdof(XAdynamicStates;field=:t3,nodID=[nodid[n₀+2n₁+n₂+1],nodid[n₀+2n₁+2n₂+1]])
# optimXAt1Vib = getdof(optimXAdynamicStates;field=:t1,nodID=[nodid[n₀+2n₁+n₂+1],nodid[n₀+2n₁+2n₂+1]])
# optimXAt2Vib = getdof(optimXAdynamicStates;field=:t2,nodID=[nodid[n₀+2n₁+n₂+1],nodid[n₀+2n₁+2n₂+1]])
# optimXAt3Vib = getdof(optimXAdynamicStates;field=:t3,nodID=[nodid[n₀+2n₁+n₂+1],nodid[n₀+2n₁+2n₂+1]])


# req = @request gp(resultants(fᵢ))
# out = getresult(dynamicStates,req,[eleid[1]])

fig2      = Figure(size = (1000,1000))
ax1 = Axis(fig2[1,1])
lines!(ax1,timeVec,t1Vib[1,:]   ,label="target")
lines!(ax1,timeVec,XAt1Vib[1,:] ,label="detuned XA")
ax2 = Axis(fig2[2,1])
lines!(ax2,timeVec,t2Vib[1,:]   ,label="target")
lines!(ax2,timeVec,XAt1Vib[1,:] ,label="detuned XA")
ax3 = Axis(fig2[3,1])
lines!(ax3,timeVec,t3Vib[1,:]   ,label="target")
lines!(ax3,timeVec,XAt1Vib[1,:], label="detuned XA")
display(fig2)



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