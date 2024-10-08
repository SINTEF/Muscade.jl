module TestFloaterMotions
using Test, Muscade, Muscade.FloaterMotions,StaticArrays, LinearAlgebra,Interpolations,GLMakie

D = length(floatermotion)
fold(x::SVector{6}) = SMatrix{3,3}( x[1],x[2],x[3],
                                    x[2],x[4],x[5],
                                    x[3],x[5],x[6])

# Define stiffness, damping and mass matrix for the true system
# K         = fold(SVector{6}([1.0,    0.0,     -0.3,     1.0,     0.0,     1.0]))
# C         = fold(SVector{6}([0.25,    0.0,     0.05,     0.15,     0.0,     0.03]))
# M         = fold(SVector{6}([1.0,    0.0,     -0.1,     0.5,     0.0,     0.3]))
K         = fold(SVector{6}([1.0,    0.0,    0.0,     1.0,     0.0,     5.0]))
C         = fold(SVector{6}([0.25,   0.0,    0.0,     0.5,     0.0,     0.5]))
M         = fold(SVector{6}([1.0,    0.0,    0.0,     5.00,     0.0,    1.0]))

#Solve direct problem
model     = Model(:MooredFloater)
n1        = addnode!(model,𝕣[0,0,0])  
e1        = addelement!(model,FloaterOnCalmWater,[n1]; K,C,M)
initialstate    = initialize!(model;time=0.)
initialstate    = setdof!(initialstate,[-.5];   field=:surge,   nodID=[n1], order=0)                                          
initialstate    = setdof!(initialstate,[.8];    field=:sway,    nodID=[n1], order=0)                                          
initialstate    = setdof!(initialstate,[-5.3];  field=:yaw,     nodID=[n1], order=0)                                          
T               = 0.1 *(1:500)
state           = solve(SweepX{2};  initialstate,time= T,verbose=false)
surge   = [s.X[1][1] for s∈state]
sway    = [s.X[1][2] for s∈state]
yaw     = [s.X[1][3] for s∈state]

# Create fake measurements
surgeMeas = surge   + .0 * randn(length(T))
swayMeas = sway     + .0 * randn(length(T))
yawMeas = yaw       + .0 * randn(length(T))

# Create intial guesses for M and C
maxDevToModel = 0.0; 
devToModelM = 1. .+ maxDevToModel*((rand(6).-.5)*2); 
devToModelC = 1. .+ maxDevToModel*((rand(6).-.5)*2);
Mguess = M .* fold(devToModelM[@SVector [i for i∈1:6 ]])
Cguess = C .* fold(devToModelC[@SVector [i for i∈1:6 ]])

# Create XUA model
modelXUA     = Model(:MooredFloater)
n1        = addnode!(modelXUA,𝕣[0,0,0])  
e1        = addelement!(modelXUA,FloaterOnCalmWater,[n1]; K,C=Cguess,M=Mguess)

# Assign costs to unknown forces
Quu       = @SVector [0.1 ^-2 for i=1:3 ]  
e2        = [addelement!(modelXUA,SingleDofCost     ,[n1]; class=:U,field=f           ,    cost=(u,t)-> 0.5*Quu[i]*u^2)  for (i,f)∈enumerate(floatermotion)]
# Assign costs to variations of model parameters (wrt guess).
fac       = [100,90,80,70,60,50,40,30,20,10,1] 
QCaa      = @SVector [.1 ^-2 for i=1:6 ]  
e3        = [addelement!(modelXUA,SingleDecayAcost  ,[n1];          field=f,fac,           cost=(a  )-> 0.5*QCaa[i]/length(T)*a^2) for (i,f)∈enumerate((:C11,:C12,:C16,:C22,:C26,:C66))] 
QMaa      = @SVector [.1 ^-2 for i=1:6 ]  
e4        = [addelement!(modelXUA,SingleDecayAcost  ,[n1];          field=f,fac,           cost=(a  )-> 0.5*QMaa[i]/length(T)*a^2) for (i,f)∈enumerate((:M11,:M12,:M16,:M22,:M26,:M66))] 
# Assign costs to measurement errors
surgeInt    = linear_interpolation(T, surgeMeas)
swayInt     = linear_interpolation(T, swayMeas)
yawInt      = linear_interpolation(T, yawMeas)
@once devSurge(surge,t)     = 5e-2 ^-2 * (surge-surgeInt(t))^2
@once devSway(sway,t)       = 5e-2 ^-2 * (sway-swayInt(t))^2
@once devYaw(yaw,t)         = 5e-2 ^-2 * (yaw-yawInt(t))^2
e5             = addelement!(modelXUA,SingleDofCost,[n1];class=:X,field=:surge,    cost=devSurge)
e6             = addelement!(modelXUA,SingleDofCost,[n1];class=:X,field=:sway,     cost=devSway)
e7             = addelement!(modelXUA,SingleDofCost,[n1];class=:X,field=:yaw,      cost=devYaw)

#Setting scale to improve convergence
myScaling = (   X=(surge=1.,sway=1.,yaw=1.),
                A=( C11=1e-1,C12=1e-1,C16=1e-1,C22=1e-1,C26=1e-1,C66=1e-1,
                    M11=1e-1,M12=1e-1,M16=1e-1,M22=1e-1,M26=1e-1,M66=1e-1))
setscale!(modelXUA;scale=myScaling,Λscale=1e4)

#Solve inverse problem
initialstateXUA    = initialize!(modelXUA;time=0.)
stateXUA         = solve(DirectXUA{2,0,1};initialstate=initialstateXUA,time=T,maxiter=100,saveiter=true,
                        maxΔx=1e0,maxΔλ=Inf,maxΔu=1e2,maxΔa=1e0)
niter=findlastassigned(stateXUA)

# Fetch and display estimated model parameters
Mest      = Mguess .* fold(exp10.(SVector{6}(stateXUA[niter][1].A[1:6 ]))) 
Cest      = Cguess .* fold(exp10.(SVector{6}(stateXUA[niter][1].A[7:12]))) 
@show Mguess ./ M
@show Mest ./ M
@show Cguess ./ C
@show Cest ./ C;

# Fetch response and loads 
surgeRec   = [s.X[1][1] for s∈stateXUA[niter]]
swayRec    = [s.X[1][2] for s∈stateXUA[niter]]
yawRec     = [s.X[1][3] for s∈stateXUA[niter]]
surgeExtF   = [s.U[1][1] for s∈stateXUA[niter]]
swayExtF    = [s.U[1][2] for s∈stateXUA[niter]]
yawExtF     = [s.U[1][3] for s∈stateXUA[niter]]
req = @request r₂,r₁,r₀  
loads = getresult(stateXUA[niter],req,[e1])
inertiaLoads = [loads[i][:r₂] for i∈1:length(T)]
dampingLoads = [loads[i][:r₁] for i∈1:length(T)]
stiffnessLoads = [loads[i][:r₀] for i∈1:length(T)]

# Display response
fig      = Figure(size = (2000,1000))
display(fig) # open interactive window (gets closed down by "save")
ax=Axis(fig[1,1], ylabel="Surge [m]",        yminorgridvisible = true,xminorgridvisible = true)
lines!(fig[1,1],T,surge,                          label="Exact direct solution")
scatter!(fig[1,1],T,surgeMeas,                          label="Simulated measurements")
lines!(fig[1,1],T,surgeRec,                       label="Inverse solution")
axislegend()

ax=Axis(fig[2,1], ylabel="Sway [m]",        yminorgridvisible = true,xminorgridvisible = true)
lines!(fig[2,1],T,sway,                         label="Measurement")
scatter!(fig[2,1],T,swayMeas,                          label="Simulated measurements")
lines!(fig[2,1],T,swayRec,                       label="Inverse solution")

ax=Axis(fig[3,1], ylabel="Yaw [deg]",       yminorgridvisible = true,xminorgridvisible = true,xlabel="Time [s]")
lines!(fig[3,1],T,yaw,                        label="Measurement")
scatter!(fig[3,1],T,yawMeas,                          label="Simulated measurements")
lines!(fig[3,1],T,yawRec,                       label="Inverse solution")

# Display loads 
ax=Axis(fig[1,2], ylabel="Surge force [N]",        yminorgridvisible = true,xminorgridvisible = true)
lines!(fig[1,2],T,-[inertiaLoads[i][1] for i∈1:length(T)], label="Inertia")
lines!(fig[1,2],T,-[dampingLoads[i][1] for i∈1:length(T)], label="Damping")
lines!(fig[1,2],T,-[stiffnessLoads[i][1] for i∈1:length(T)], label="Stiffness")
lines!(fig[1,2],T,surgeExtF,                       label="Unknown")
axislegend()

ax=Axis(fig[2,2], ylabel="Sway force [N]",        yminorgridvisible = true,xminorgridvisible = true)
lines!(fig[2,2],T,-[inertiaLoads[i][2] for i∈1:length(T)], label="Inertia")
lines!(fig[2,2],T,-[dampingLoads[i][2] for i∈1:length(T)], label="Damping")
lines!(fig[2,2],T,-[stiffnessLoads[i][2] for i∈1:length(T)], label="Stiffness")
lines!(fig[2,2],T,swayExtF,                       label="Unknown")

ax=Axis(fig[3,2], ylabel="Yaw moment [Nm]",       yminorgridvisible = true,xminorgridvisible = true,xlabel="Time [s]")
lines!(fig[3,2],T,-[inertiaLoads[i][3] for i∈1:length(T)], label="Inertia")
lines!(fig[3,2],T,-[dampingLoads[i][3] for i∈1:length(T)], label="Damping")
lines!(fig[3,2],T,-[stiffnessLoads[i][3] for i∈1:length(T)], label="Stiffness")
lines!(fig[3,2],T,yawExtF,                       label="Unknown")
display(fig) # open interactive window (gets closed down by "save")

end
;