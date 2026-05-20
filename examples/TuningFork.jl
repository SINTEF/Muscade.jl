# # Tuning vibrations of a tuning fork
# 
# Tuning forks should resonate at a desired frequency when subjected to an impulsive load. Achieving the right frequency cab be done by adding material or filing material off the prongs.
# We start from a target solution establihed using SweepX. A sprious mass is then introduced, parametrized by an A-dof. We use the SweepXA and DirectXUA solvers to estimate how much mass should be removed. In SweepXA, the excitation is assume to be known, while it is estimated by DirectXUA. 

using Muscade, StaticArrays, GLMakie, Muscade.Toolbox

# Dimensions of the tuning fork and number of elements
l₀=2.5e-2;  n₀=1    # base - length and number of elements 
l₁=1e-2;    n₁=1    # fork - length and number of elements per side
l₂=9.5e-2;  n₂=5    # prong - length and number of elements per side
h = 5e-3;           # width of the tuning fork cross-section
b = 5e-3;           # thickness of the cross-section
E = 210e9;          # Young modulus
G = 79.3e9;         # shear modulus
μ = 7850. *h*b;     # mass per unit length, assumes steel
beamCrossSection = BeamCrossSection(EA=E*h*b, EI₂=E*h^3*b/12, EI₃=E*b*h^3/12, GJ=G*h*b^3/3, μ=μ, ι₁=μ/12*(h^2+b^2),Cl₂=1.,Cl₃=1.); 

# Approximation of the fork fundamental frequency (cantilever beam bending)
f₁ = 1.875^2/(2*π*l₂^2)*√(E*h^3*b/12/μ) 

# Define an A-parametrized point mass element (inducing an inertia force)
struct LumpedMass <: AbstractElement; m :: 𝕣; end
LumpedMass(nod::Vector{Node};m::𝕣) = LumpedMass(m)
@espy function Muscade.residual(o::LumpedMass, X,U,A, t,SP,dbg); return SVector((o.m+A[1]).*∂2(X)),noFB; end
Muscade.doflist( ::Type{LumpedMass})  = (inod =(1,1,1,2), class=(:X,:X,:X,:A), field=(:t1,:t2,:t3,:mass));

# Time step and simulation length settings
Δt = 1/f₁ / 8;
nDynamicLoadSteps = 25*8;
timeVec = Δt:Δt:nDynamicLoadSteps*Δt;

# Coordinates of the nodes.
XnodeCoord = vcat(
    hcat((0:n₀) * l₀/n₀     , zeros(n₀+1,1) , zeros(n₀+1,1) ), # base nodes
    hcat(ones(n₁,1)*l₀      , (1:n₁)*l₁/n₁  , zeros(n₁,1)   ), # upper fork
    hcat(ones(n₁,1)*l₀      , -(1:n₁)*l₁/n₁ , zeros(n₁,1)   ), # lower fork
    hcat(l₀ .+ (1:n₂)*l₂/n₂ , ones(n₂,1)*l₁ , zeros(n₂,1)   ), # upper prong 
    hcat(l₀ .+ (1:n₂)*l₂/n₂ , -ones(n₂,1)*l₁, zeros(n₂,1)   )  # lower prong 
)

# External excitation load at the end of the upper prong 
pling = Muscade.FunctionFromVector([-10, 0, 5/f₁, (5+1/2)/f₁, 10],[0,0,10,0,0]); @functor with() pling_(t) = pling(t)
plingNode = n₀+2n₁+n₂+1; # where do we hit the tuning fork

# # Building the model
model       = Model(:TuningFork);

# Add the nodes describing the beam model, the node that will contain the parameter to optimize, and  the node that will contain the load to estimate in DirectXUA
Xnod = addnode!(model,XnodeCoord)
Anod = addnode!(model,Vector{𝕣}(undef,3))
Unod = addnode!(model,Vector{𝕣}(undef,3)); 

# List of nodes that will be used to create the beam elements
nel = n₀+2n₁+2n₂
mesh = vcat(
    hcat(Xnod[1:n₀],Xnod[2:n₀+1]),                                                                              # base 
    hcat(Xnod[n₀+1:n₀+n₁],Xnod[n₀+2:n₀+n₁+1]),                                                                  # upper fork 
    [Xnod[n₀+1] Xnod[n₀+n₁+2]],         hcat(Xnod[n₀+n₁+2:n₀+2n₁],Xnod[n₀+n₁+3:n₀+2n₁+1])           ,           # lower fork
    [Xnod[n₀+n₁+1] Xnod[n₀+2n₁+2]],     hcat(Xnod[n₀+2n₁+2:n₀+2n₁+n₂],Xnod[n₀+2n₁+3:n₀+2n₁+n₂+1])   ,           # upper prong
    [Xnod[n₀+2n₁+1] Xnod[n₀+2n₁+n₂+2]], hcat(Xnod[n₀+2n₁+n₂+2:n₀+2n₁+2n₂],Xnod[n₀+2n₁+n₂+3:n₀+2n₁+2n₂+1])       # lower prong
)

# In the XUA analysis we estimate the distributed load on the element close to the node where the load was applied
plingElement = n₀+2n₁+n₂;

# Add beam elements to the model, enabling Udof only for the element that is hit. 
eleid = Vector{Muscade.EleID}(undef,nel)
for idx=1:nel
    if idx == plingElement
       eleid[idx] = addelement!(model,EulerBeam3D{true}, vcat(mesh[idx,:],Unod); mat=beamCrossSection,orient2=SVector(0.,0,1))
    else
        eleid[idx] = addelement!(model,EulerBeam3D{false}, mesh[idx,:];             mat=beamCrossSection,orient2=SVector(0.,0,1))
    end
end
[addelement!(model,Hold,[Xnod[1]]  ;field)  for field∈[:t1,:t2,:t3,:r1,:r2,:r3]]; 

# # Create variations of the model
XAmodel =   deepcopy(model)
XUAmodel =  deepcopy(model);

# Add parasitic lump mass on both XA and XUA models
spuriousMass = 1e-3;
spuriousMassLocation = plingNode-1
addelement!(XAmodel, LumpedMass,  [Xnod[spuriousMassLocation] Anod]; m=spuriousMass)
addelement!(XUAmodel, LumpedMass, [Xnod[spuriousMassLocation] Anod]; m=spuriousMass);

# Add known external load to X and XA models
addelement!(model,   DofLoad,[Xnod[plingNode]];field=:t2,value=pling_)
addelement!(XAmodel, DofLoad,[Xnod[plingNode]];field=:t2,value=pling_);

# Establish target response
initialState    = initialize!(model;time=0.);
dynamicStates   = solve(SweepX{2};initialstate=initialState, time=timeVec, verbose=false)
vibTarget       = getdof(dynamicStates;field=:t2,nodID=[Xnod[plingNode]])
target          = Muscade.FunctionFromVector(timeVec,vibTarget);

# Add costs on the deviation to target measurements, in the XA and XUA models
@functor with(σᵥ=1e-6,target) Xcost(x,t)=((x-target(t))/σᵥ)^2
addelement!(XAmodel, SingleDofCost,[Xnod[plingNode]]; class=:X, field=:t2, cost=Xcost)
addelement!(XUAmodel,SingleDofCost,[Xnod[plingNode]]; class=:X, field=:t2, cost=Xcost);

# Add costs on the coorecting mass, in the XA and XUA models
@functor with(σₘ=5e-3,timeVec)          Acost_(a) = 0.5*(a/σₘ)^2/length(timeVec)
addelement!(XAmodel,  SingleAcost,[Anod]; field=:mass, cost=Acost_);

# Prioritize estimation of U-loads during the first DirectXUA iterations, rather than updates of the A-parameter. 
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
n=5; fac = [2^(n-i) for i∈1:n]
addelement!(XUAmodel,SingleDecayAcost  ,[Anod]; field=:mass,fac, cost=Acost_);
#src addelement!(XUAmodel,SingleAcost,       [Anod]; field=:mass, cost=Acost_)

# Add cost to unknown loads XUA model, force per unit legnth in the direction of the excitation 
@functor with(σᵤₚ=1e-1,f₁,l₂,n₂) UcostPling(u,t) = 0.5*(u*(l₂/n₂)/σᵤₚ)^2*exp(t*f₁/5)
addelement!(XUAmodel,SingleDofCost,[Unod]; class=:U,field=:t2, cost=UcostPling);

# and along the orthogonal dofs
@functor with(σᵤₙ=1e-6) UcostNoLoad(u,t) = 0.5*(u/σᵤₙ)^2
[addelement!(XUAmodel,SingleDofCost,[Unod]; class=:U,field=field, cost=UcostNoLoad) for field∈(:t1, :t3)];

# # Run analyses

# XA model (before estimating the mass)
XAinitialState    = initialize!(XAmodel;time=0.);
XAdynamicStates   = solve(SweepX{2};initialstate=XAinitialState, time=timeVec, verbose=false); 

# XA model (after having estimated the mass)
optimXAstate  = solve(SweepXA{2}; initialstate=XAinitialState, time=timeVec, 
                maxAiter=20,maxΔa=1e-10, 
                verbose=false);

# DirectXUA (slack convergence criteria and number of iterations)
XUAinitialState    = initialize!(XUAmodel;time=0.);
optimXUAstate   = solve(DirectXUA{2,0,1};initialstate=[XUAinitialState], time=[timeVec],
                maxiter=15, maxΔx=1e-1,maxΔλ=Inf,maxΔu=1e1,maxΔa=1e-1, 
                verbose=false);


#  Gather results for comparison
analysisResults = (dynamicStates, XAdynamicStates,optimXAstate,optimXUAstate[end]); 
t1 = Matrix{Float64}(undef,length(analysisResults),length(timeVec))
t2 = Matrix{Float64}(undef,length(analysisResults),length(timeVec))
t3 = Matrix{Float64}(undef,length(analysisResults),length(timeVec))
transmittedForce = Matrix{Float64}(undef,length(analysisResults),length(timeVec))
req = @request gp(resultants(fᵢ))
for idx ∈ 1:length(analysisResults)
    t1[idx,:] = getdof(analysisResults[idx];field=:t1,nodID=[Xnod[plingNode]])[:]
    t2[idx,:] = getdof(analysisResults[idx];field=:t2,nodID=[Xnod[plingNode]])[:]
    t3[idx,:] = getdof(analysisResults[idx];field=:t3,nodID=[Xnod[plingNode]])[:]
 end
excEstt2 = getdof(optimXUAstate[end];class=:U,field=:t2,nodID=[Unod])[:]*l₂/n₂;

# Text output
println("Estimated mass to remove from SweepXA analysis [g]:")
println(-getdof(optimXAstate[end];class=:A,field=:mass,nodID=[Anod])[1]*1e3)
println("Estimated mass to remove from DirectXUA analysis [g]:")
println(-getdof(optimXUAstate[end];class=:A,field=:mass,nodID=[Anod])[1]*1e3)
println("Expected [g]:")
println(spuriousMass*1e3)

# Plot the axial and lateral displacement of the control node, and the estimated force.
fig      = Figure(size = (1000,1000))
ax0 = Axis(fig[1,1],ylabel="Excitation [N]")
lines!(ax0,timeVec,pling.(timeVec), label="Actual",color=:black)
lines!(ax0,timeVec,excEstt2, label="Estimated (XUA)")
axislegend(ax0)
ax1 = Axis(fig[2,1],ylabel="xₑ (axial) [mm]")
scatter!(ax1,timeVec,t1[1,:]*1e3, label="Target")
lines!(ax1,timeVec,t1[2,:]*1e3, label="Detuned config.")
lines!(ax1,timeVec,t1[3,:]*1e3, label="Tuned (XA)")
lines!(ax1,timeVec,t1[4,:]*1e3, label="Tuned (XUA)")
axislegend(ax1)
ax2 = Axis(fig[3,1],ylabel="yₑ (lateral) [mm]")
scatter!(ax2,timeVec,t2[1,:]*1e3, label="Target")
lines!(ax2,timeVec,t2[2,:]*1e3, label="Detuned config.")
lines!(ax2,timeVec,t2[3,:]*1e3, label="Tuned (XA)")
lines!(ax2,timeVec,t2[4,:]*1e3, label="Tuned (XUA)")

currentDir = @__DIR__
if occursin("build", currentDir)
    save(normpath(joinpath(currentDir,"..","src","assets","tuningFork.png")),fig)
elseif occursin("examples", currentDir)
    save(normpath(joinpath(currentDir,"tuningFork.png")),fig)
end

# ![Result](assets/tuningFork.png)