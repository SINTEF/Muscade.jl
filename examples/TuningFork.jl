using Muscade, StaticArrays, GLMakie, Muscade.Toolbox

# Dimensions of the tuning fork and number of elements
l₀=2.5e-2; n₀=1       # base - length and number of elements 
l₁=1e-2; n₁=1       # fork - length and number of elements per side
l₂=9.5e-2; n₂=5      # prong - length and number of elements per side
h = 5e-3;           # width of the tuning tuning fork
b = 5e-3;           # thickness 
E = 210e9;          # Young modulus
G = 79.3e9;         # shear modulus
μ = 7850. *h*b;     # mass per unit length, assumes steel
beamCrossSection = BeamCrossSection(EA=E*h*b, EI₂=E*h^3*b/12, EI₃=E*b*h^3/12, GJ=G*h*b^3/3, μ=μ, ι₁=μ/12*(h^2+b^2),Cl₂=1.,Cl₃=1.)

# Approximation of the fork fundamental frequency (cantilever beam bending)
f₁ = 1.875^2/(2*π*l₂^2)*√(E*h^3*b/12/μ)

# Define parametrized point mass element (inducing an inertia force)
struct LumpedMass <: AbstractElement; m :: 𝕣; end
LumpedMass(nod::Vector{Node};m::𝕣) = LumpedMass(m)
@espy function Muscade.residual(o::LumpedMass, X,U,A, t,SP,dbg); return SVector((o.m+A[1]).*∂2(X)),noFB; end
Muscade.doflist( ::Type{LumpedMass})  = (inod =(1,1,1,2), class=(:X,:X,:X,:A), field=(:t1,:t2,:t3,:mass))

# Time steps
Δt = 1/f₁ / 8;
nDynamicLoadSteps = 25*8;
timeVec = Δt:Δt:nDynamicLoadSteps*Δt;

# Meshing: listing nodes
XnodeCoord = vcat(
    hcat((0:n₀) * l₀/n₀     , zeros(n₀+1,1) , zeros(n₀+1,1) ), # base nodes
    hcat(ones(n₁,1)*l₀      , (1:n₁)*l₁/n₁  , zeros(n₁,1)   ), # upper fork
    hcat(ones(n₁,1)*l₀      , -(1:n₁)*l₁/n₁ , zeros(n₁,1)   ), # lower fork
    hcat(l₀ .+ (1:n₂)*l₂/n₂ , ones(n₂,1)*l₁ , zeros(n₂,1)   ), # upper prong 
    hcat(l₀ .+ (1:n₂)*l₂/n₂ , -ones(n₂,1)*l₁, zeros(n₂,1)   )  # lower prong 
)

# External excitation load
pling = Muscade.FunctionFromVector([-10, 0, 5/f₁, (5+1/2)/f₁, 10],[0,0,10,0,0]); @functor with() pling_(t) = pling(t)
plingNode = n₀+2n₁+n₂+1; # where do we hit the tuning fork

model       = Model(:TuningFork);

Xnodid = addnode!(model,XnodeCoord)
# Add the node that will contain the parameter to optimize
Anod = addnode!(model,[-1., 0, 0])
# Add the node that will contain the load to estimate
Unodid = addnode!(model,[-2., 0, 0])

# Add beam elements
mesh = vcat(
    hcat(Xnodid[1:n₀],Xnodid[2:n₀+1]),                                                                                # base 
    hcat(Xnodid[n₀+1:n₀+n₁],Xnodid[n₀+2:n₀+n₁+1]),                                                                    # upper fork 
    [Xnodid[n₀+1] Xnodid[n₀+n₁+2]],         hcat(Xnodid[n₀+n₁+2:n₀+2n₁],Xnodid[n₀+n₁+3:n₀+2n₁+1])           ,           # lower fork
    [Xnodid[n₀+n₁+1] Xnodid[n₀+2n₁+2]],     hcat(Xnodid[n₀+2n₁+2:n₀+2n₁+n₂],Xnodid[n₀+2n₁+3:n₀+2n₁+n₂+1])   ,           # upper prong
    [Xnodid[n₀+2n₁+1] Xnodid[n₀+2n₁+n₂+2]], hcat(Xnodid[n₀+2n₁+n₂+2:n₀+2n₁+2n₂],Xnodid[n₀+2n₁+n₂+3:n₀+2n₁+2n₂+1])       # lower prong
)
nel = n₀+2n₁+2n₂
plingElement = plingNode-1

eleid = Vector{Muscade.EleID}(undef,nel)
for idx=1:nel
    if idx == plingElement
       eleid[idx] = addelement!(model,EulerBeam3D{true}, vcat(mesh[idx,:],Unodid); mat=beamCrossSection,orient2=SVector(0.,0,1))
    else
        eleid[idx] = addelement!(model,EulerBeam3D{false}, mesh[idx,:];             mat=beamCrossSection,orient2=SVector(0.,0,1))
    end
end
[addelement!(model,Hold,[Xnodid[1]]  ;field)  for field∈[:t1,:t2,:t3,:r1,:r2,:r3]]; 

# Create variations of the model
XAmodel =   deepcopy(model)
XUAmodel =  deepcopy(model)

# Add parasitic lump mass on both XA and XUA models
spuriousMass = 2e-3;
spuriousMassLocation = plingNode-1
addelement!(XAmodel, LumpedMass,  [Xnodid[spuriousMassLocation] Anod]; m=spuriousMass)
addelement!(XUAmodel, LumpedMass, [Xnodid[spuriousMassLocation] Anod]; m=spuriousMass)

# Add external load to X and XA models
addelement!(model,   DofLoad,[Xnodid[plingNode]];field=:t2,value=pling_)
addelement!(XAmodel, DofLoad,[Xnodid[plingNode]];field=:t2,value=pling_)
addelement!(XUAmodel,DofLoad,[Xnodid[plingNode]];field=:t2,value=pling_)

# Establish target response
initialState    = initialize!(model;time=0.);
staticState     = solve(SweepX{0};initialstate=initialState, time=[0.])
dynamicStates   = solve(SweepX{2};initialstate=staticState[end], time=timeVec)
vibTarget       = getdof(dynamicStates;field=:t2,nodID=[Xnodid[plingNode]])
target          = Muscade.FunctionFromVector(timeVec,vibTarget)

# Add costs on the deviation to target measurements, XA and XUA model
@functor with(σᵥ=1e-6,target) Xcost(x,t)=((x-target(t))/σᵥ)^2
addelement!(XAmodel, SingleDofCost,[Xnodid[plingNode]]; class=:X, field=:t2, cost=Xcost)
addelement!(XUAmodel,SingleDofCost,[Xnodid[plingNode]]; class=:X, field=:t2, cost=Xcost)

# Add costs on the coorecting mass, XA and XUA model
@functor with(σₘ=1e-3)          Acost(a) = 0.5*(a/σₘ)^2
addelement!(XAmodel,  SingleAcost,[Anod]; field=:mass, cost=Acost)

struct SingleDecayAcost{Field,Tcost,Tcostargs} <: AbstractElement
    cost     :: Tcost     
    costargs :: Tcostargs
    fac      :: 𝕣1
end
SingleDecayAcost(nod::Vector{Node};field::Symbol,fac,cost::Functor ,costargs=()) = SingleDecayAcost{field,typeof(cost),typeof(costargs)}(cost,costargs,fac)
Muscade.doflist(::Type{<:SingleDecayAcost{Field,Tcost,Tcostargs}}) where{Field,Tcost,Tcostargs} = (inod=(1,),class=(:A,),field=(Field,))
@espy function Muscade.lagrangian(o::SingleDecayAcost,Λ,X,U,A,t,SP,dbg) 
    iter  = min(length(o.fac),default{:iter}(SP,length(o.fac)))
    ☼cost = o.cost(    A[1]  ,o.costargs...)
    return cost*o.fac[iter],noFB
end
n=15; fac = [2^(n-i) for i∈1:n]
addelement!(XUAmodel,SingleDecayAcost  ,[Anod]; field=:mass,fac, cost=Acost)
# addelement!(XUAmodel,SingleAcost,       [Anod]; field=:mass, cost=Acost)

# Add cost to unknown loads XUA model
@functor with(σᵤₚ=1e-1,f₁) UcostPling(u,t) = 0.5*(u/σᵤₚ)^2*exp(t*f₁/10)
addelement!(XUAmodel,SingleDofCost,[Unodid]; class=:U,field=:t2, cost=UcostPling)

@functor with(σᵤₙ=1e-6) UcostNoLoad(u,t) = 0.5*(u/σᵤₙ)^2
[addelement!(XUAmodel,SingleDofCost,[Unodid]; class=:U,field=field, cost=UcostNoLoad) for field∈(:t1, :t3)]

# Run analyses

# XA model (before optimization)
XAinitialState    = initialize!(XAmodel;time=0.);
XAdynamicStates   = solve(SweepX{2};initialstate=XAinitialState, time=timeVec)
# XA model (after optimization)
optimXAstate  = solve(SweepXA{2,2,2}; initialstate=XAinitialState,verbose=false,time=timeVec, maxAiter=20,maxΔa=1e-10);

# # DirectXUA
XUAinitialState    = initialize!(XUAmodel;time=0.);
XUAstaticState     = solve(SweepX{0};initialstate=XUAinitialState, time=[timeVec[1]])
# # XUAdynamicStates   = solve(SweepX{2};initialstate=XUAinitialState, time=timeVec)
# optimXUAstate   = solve(DirectXUA{2,0,1};initialstate=[XUAstaticState[end]], 
#         saveiter=true, maxiter=10, time=[timeVec],
#         maxΔx=1e-6,maxΔλ=Inf,maxΔu=1e-5,maxΔa=1e-10);


# figXUA      = Figure(size = (1000,1000))
# ax = Axis3(figXUA[1,1],aspect=:data)
# for (idx,t) ∈ enumerate(timeVec)
# draw!(ax,optimXUAstate[5][1][idx])
# end
# zlims!(ax,-.001,.001)
# display(figXUA)

# Contents of the three models
Muscade.describe(model,:eletyp)
Muscade.describe(XAmodel,:eletyp)
Muscade.describe(XUAmodel,:eletyp)

# Gather and compare results (add tuned XA results when ready)
analysisResults = (dynamicStates, XAdynamicStates,optimXAstate)#,optimXUAstate[10][1]); 
t1 = Matrix{Float64}(undef,length(analysisResults),length(timeVec))
t2 = Matrix{Float64}(undef,length(analysisResults),length(timeVec))
t3 = Matrix{Float64}(undef,length(analysisResults),length(timeVec))
transmittedForce = Matrix{Float64}(undef,length(analysisResults),length(timeVec))
req = @request gp(resultants(fᵢ))
for idx ∈ 1:length(analysisResults)
    t1[idx,:] = getdof(analysisResults[idx];field=:t1,nodID=[Xnodid[plingNode]])[:]
    t2[idx,:] = getdof(analysisResults[idx];field=:t2,nodID=[Xnodid[plingNode]])[:]
    t3[idx,:] = getdof(analysisResults[idx];field=:t3,nodID=[Xnodid[plingNode]])[:]
 end

println("Estimated mass to remove [g]:")
println(-getdof(optimXAstate[end];class=:A,field=:mass,nodID=[Anod])[1]*1e3)
println("Expected [g]:")
println(spuriousMass*1e3)

# figExc      = Figure(size = (1000,1000))
# ax = Axis(figExc[1,1])
# for idx ∈ 1:5
#     excEst1 = getdof(optimXUAstate[idx][1];class=:U,field=:t2,nodID=[Unodid])[:]
#     lines!(ax,timeVec,excEst1,color=:black)
# end
# display(figExc)

# excEst = getdof(optimXUAstate[10][1];class=:U,field=:t2,nodID=[Unodid])[:]

# Plot
fig2      = Figure(size = (1000,1000))
ax5 = Axis(fig2[1,1],ylabel="Excitation [N]")
lines!(ax5,timeVec,pling.(timeVec), label="Actual",color=:black)
# lines!(ax5,timeVec,excEst, label="Estimated")
ax1 = Axis(fig2[2,1],ylabel="xₑ (axial) [mm]")
scatter!(ax1,timeVec,t1[1,:]*1e3, label="target")
lines!(ax1,timeVec,t1[2,:]*1e3, label="detuned XA")
lines!(ax1,timeVec,t1[3,:]*1e3, label="tuned XA")
# lines!(ax1,timeVec,t1[4,:]*1e3, label="detuned XUA")
axislegend(ax1)
ax2 = Axis(fig2[3,1],ylabel="yₑ (lateral) [mm]")
scatter!(ax2,timeVec,t2[1,:]*1e3, label="target")
lines!(ax2,timeVec,t2[2,:]*1e3, label="detuned XA")
lines!(ax2,timeVec,t2[3,:]*1e3, label="tuned XA")
ax2 = Axis(fig2[4,1],ylabel="zₑ (out-of-plane) [mm]")
scatter!(ax2,timeVec,t3[1,:]*1e3, label="target")
lines!(ax2,timeVec,t3[2,:]*1e3, label="detuned XA")
lines!(ax2,timeVec,t3[3,:]*1e3, label="tuned XA")
display(fig2)



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