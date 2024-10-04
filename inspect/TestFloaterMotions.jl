# module TestFloaterMotions
# using Profile,ProfileView,BenchmarkTools
using Test, Muscade, Muscade.FloaterMotions,StaticArrays, LinearAlgebra,Interpolations,GLMakie

D = length(floatermotion)
fold(x::SVector{6}) = SMatrix{3,3}( x[1],x[2],x[3],
                                    x[2],x[4],x[5],
                                    x[3],x[5],x[6])
K         = fold(SVector{6}([1.0,    0.0,     -0.3,     1.0,     0.0,     1.0]))
C         = fold(SVector{6}([0.3,    0.0,     0.1,     0.2,     0.0,     0.05]))
M         = fold(SVector{6}([1.0,    0.0,     -0.1,     0.5,     0.0,     0.3]))
model     = Model(:MooredFloater)
n1        = addnode!(model,𝕣[0,0,0])  
e1        = addelement!(model,FloaterOnCalmWater,[n1]; K,C,M)
initialstate    = Muscade.State{1,3,1,Nothing}(initialize!(model;time=0.))  # recast to force the state to have 2nd derivatives, 
initialstate.X[2][1] = 0.                                           # so we can set initial velocity
initialstate.X[1][1] = 1.                                           # so we can set initial position
initialstate.X[2][2] = 0.                                           # so we can set initial velocity
initialstate.X[1][2] = 2.                                           # so we can set initial position
initialstate.X[2][3] = 0.                                           # so we can set initial velocity
initialstate.X[1][3] = -1.3                                           # so we can set initial position
T               = 0.15 *(1:200)
state           = solve(SweepX{2};  initialstate,time= T,verbose=false)
surge   = [s.X[1][1] for s∈state]
sway    = [s.X[1][2] for s∈state]
yaw     = [s.X[1][3] for s∈state]

# Create fake measurements
surgeMeas = surge   + .05 * randn(length(T))
swayMeas = sway     + .05 * randn(length(T))
yawMeas = yaw       + .05 * randn(length(T))


# Create intial guesses for M and C
maxDevToModel = 0.5; 
devToModelM = 1. .+ maxDevToModel*((rand(6).-.5)*2); 
devToModelC = 1. .+ maxDevToModel*((rand(6).-.5)*2);
Cguess = C .* fold(devToModelC[@SVector [i for i∈1:6 ]])
Mguess = M .* fold(devToModelM[@SVector [i for i∈1:6 ]])


fac       = [100,90,80,70,60,50,40,30,20,10,1]
modelXUA     = Model(:MooredFloater)
n1        = addnode!(modelXUA,𝕣[0,0,0])  
e1        = addelement!(modelXUA,FloaterOnCalmWater,[n1]; K,C=Cguess,M=Mguess)

Quu       = @SVector [0.1 ^-2 for i=1:3 ]  
e2        = [addelement!(modelXUA,SingleDofCost     ,[n1]; class=:U,field=f           ,    cost=(u,t)-> 0.5*Quu[i]*u^2)  for (i,f)∈enumerate(floatermotion)]
QCaa      = @SVector [.1 ^-2 for i=1:6 ]  #NB: uncertainty on ΞA
e3        = [addelement!(modelXUA,SingleDecayAcost  ,[n1];          field=f,fac,           cost=(a  )-> 0.5*QCaa[i]*a^2) for (i,f)∈enumerate((:C11,:C12,:C16,:C22,:C26,:C66))] 
QMaa      = @SVector [.1 ^-2 for i=1:6 ]  
e4        = [addelement!(modelXUA,SingleDecayAcost  ,[n1];          field=f,fac,           cost=(a  )-> 0.5*QMaa[i]*a^2) for (i,f)∈enumerate((:M11,:M12,:M16,:M22,:M26,:M66))] 


surgeInt    = linear_interpolation(T, surgeMeas)
swayInt     = linear_interpolation(T, swayMeas)
yawInt      = linear_interpolation(T, yawMeas)
@once devSurge(surge,t)     = 1e-2 ^-2 * (surge-surgeInt(t))^2
@once devSway(sway,t)       = 1e-2 ^-2 * (sway-swayInt(t))^2
@once devYaw(yaw,t)         = 1e-2 ^-2 * (yaw-yawInt(t))^2
e5             = addelement!(modelXUA,SingleDofCost,[n1];class=:X,field=:surge,    cost=devSurge)
e6             = addelement!(modelXUA,SingleDofCost,[n1];class=:X,field=:sway,     cost=devSway)
e7             = addelement!(modelXUA,SingleDofCost,[n1];class=:X,field=:yaw,      cost=devYaw)
initialstateXUA    = initialize!(modelXUA;time=0.)

NDX       = 3
NDU       = 1
NA        = 1
stateXUA         = solve(DirectXUA{3,1,1};initialstate=initialstateXUA,time=T,maxiter=200,saveiter=true,
                        maxΔx=1e-3,maxΔλ=1e-2,maxΔu=1e-3,maxΔa=1e-3)


# As = [s.A[1][1] for s∈stateXUA]

# iter=10
# step=1
# idof=1
# ider=1

# # if saveiter
# stateXUA[iter][step].A[idof]
# getdof(stateXUA[iter][step],field=:surge,nodID=[n1]) 
# stateXUA[iter][step].X[ider][idof] #X[1] value, then derivative
# stateXUA[iter][step].U[ider][idof]

# if not saveiter
# stateXUA[step].A[idof]
# getdof(stateXUA[iter][step],field=:surge,nodID=[n1]) 
# stateXUA[step].X[ider][idof] #X[1] value, then derivative
# stateXUA[step].U[ider][idof]

niter = findlastassigned(stateXUA)
surgeRec   = [s.X[1][1] for s∈stateXUA[niter]]
swayRec    = [s.X[1][2] for s∈stateXUA[niter]]
yawRec     = [s.X[1][3] for s∈stateXUA[niter]]

surgeExtF   = [s.U[1][1] for s∈stateXUA[niter]]
swayExtF    = [s.U[1][2] for s∈stateXUA[niter]]
yawExtF     = [s.U[1][3] for s∈stateXUA[niter]]

# C11     = [exp10(stateXUA[niter].A[1]) for s∈stateXUA[niter]]



# Display solution 
fig      = Figure(size = (1000,1000))
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
display(fig) # open interactive window (gets closed down by "save")


# Display solution 
# fig      = Figure(size = (1000,1000))
# display(fig) # open interactive window (gets closed down by "save")
ax=Axis(fig[1,2], ylabel="Surge [N]",        yminorgridvisible = true,xminorgridvisible = true)
lines!(fig[1,2],T,surgeExtF,                       label="Inverse solution")
axislegend()
ax=Axis(fig[2,2], ylabel="Sway [N]",        yminorgridvisible = true,xminorgridvisible = true)
lines!(fig[2,2],T,swayExtF,                       label="Inverse solution")
ax=Axis(fig[3,2], ylabel="Yaw [Nm]",       yminorgridvisible = true,xminorgridvisible = true,xlabel="Time [s]")
lines!(fig[3,2],T,yawExtF,                       label="Inverse solution")
display(fig) # open interactive window (gets closed down by "save")



# # Display solution 
# fig      = Figure(size = (1000,1000))
# display(fig) # open interactive window (gets closed down by "save")
# ax=Axis(fig[1,1], ylabel="M11",        yminorgridvisible = true,xminorgridvisible = true)
# lines!(fig[1,1],T,C11,                       label="Inverse solution")
# axislegend()
# # ax=Axis(fig[2,1], ylabel="Sway [m]",        yminorgridvisible = true,xminorgridvisible = true)
# # lines!(fig[2,1],T,swayExtF,                       label="Inverse solution")
# # ax=Axis(fig[3,1], ylabel="Yaw [deg]",       yminorgridvisible = true,xminorgridvisible = true,xlabel="Time [s]")
# # lines!(fig[3,1],T,yawExtF,                       label="Inverse solution")
# display(fig) # open interactive window (gets closed down by "save")







# Profile.clear()
# Profile.@profile for i=1:1
#     local stateXUA         = solve(DirectXUA{3,1,1};initialstate=state,time=0.:3.,maxiter=5)
# end
# ProfileView.view(fontsize=30);
# end
# After clicking on a bar in the flame diagram, you can type warntype_last() and see the result of 
# code_warntype for the call represented by that bar.



# end
;