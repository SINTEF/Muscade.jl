# module TestFloaterMotions
# using Profile,ProfileView,BenchmarkTools
using Test, Muscade, Muscade.FloaterMotions,StaticArrays, LinearAlgebra

D = length(floatermotion)

K         =      SMatrix{D,D}(I(D))
C         = 0.05*SMatrix{D,D}(I(D))
M         =      SMatrix{D,D}(I(D))
Quu       = @SVector [1 for i=1:D ]  # DecayUcost only allows priors with uncorrelated U's.  Should be changed to take std dev as input, not inverse variance...
QCaa      = @SVector [1 for i=1:D ]  # same with A's
fac       = [100,10,1]

model     = Model(:SeaSick)
n1        = addnode!(model,𝕣[0,0,0])  
e1        = addelement!(model,FloaterOnCalmWater,[n1]; K,C,M)
e2        = [addelement!(model,SingleDofCost     ,[n1]; class=:U,field=f           ,    cost=(u,t)->u^2) for f∈floatermotion]
e3        = [addelement!(model,SingleDecayAcost  ,[n1];          field=f,fac,cost=(a  )->a^2) for f∈(:C11,:C12,:C16,:C22,:C26,:C66)] 
e4        = [addelement!(model,SingleDecayAcost  ,[n1];          field=f,fac,cost=(a  )->a^2) for f∈(:M11,:M12,:M16,:M22,:M26,:M66)] 
# need to add `SingleDofCost` elements to include measurements of motion
state     = initialize!(model)   

NDX       = 3
NDU       = 1
NA        = 1

stateXUA         = solve(DirectXUA{3,1,1};initialstate=state,time=0.:6.,maxiter=5)

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