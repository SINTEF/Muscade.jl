using Profile,ProfileView,BenchmarkTools
using Muscade

model           = Model(:TestModel)
n1              = addnode!(model,𝕣[ 0, 0])  # moving node
n2              = addnode!(model,𝕣[10, 0])  # anchor 1 
n3              = addnode!(model,𝕣[ 0,10])  # anchor 2 
n4              = addnode!(model,𝕣[     ])  # A-nod for springs
e1              = addelement!(model,Spring{2},[n1,n2,n4], EI=1)
e2              = addelement!(model,Spring{2},[n1,n3,n4], EI=1)
@once f1 f1(t)  = t
e3              = addelement!(model,DofLoad  ,[n1], field=:tx1      ,value=f1)
e4              = addelement!(model,DofHold  ,[n2], field=:tx1)
e5              = addelement!(model,DofHold  ,[n2], field=:tx2)
e6              = addelement!(model,DofHold  ,[n3], field=:tx1)
e7              = addelement!(model,DofHold  ,[n3], field=:tx2)
@once f2 f2(x)  = 1x^2
@once f3 f3(a)  = 0.11a^2
e8              = addelement!(model,XdofCost ,[n1], field=:tx1      ,cost=f2)
e9              = addelement!(model,XdofCost ,[n1], field=:tx2      ,cost=f2)
e10             = addelement!(model,AdofCost ,[n4], field=:ΞL₀      ,cost=f3)
e10             = addelement!(model,AdofCost ,[n4], field=:ΞEI      ,cost=f3)
stateX          = solve(staticX;  model,time=[.5,1.],verbose=false)
state           = solve(staticXUA;model,initial=stateX,maxYiter= 0,verbose=false)
#@btime state    = solve(staticXUA;model,initial=stateX,maxYiter= 0,verbose=false) 

Profile.clear()
Profile.@profile for i=1:25000
    local stateX          = solve(staticX;  model,time=[.5,1.],verbose=false)
#    local state  = solve(staticXUA;model,initial=stateX,maxYiter= 50,verbose=false);
end
ProfileView.view(fontsize=30);

# After clicking on a bar, you can type warntype_last() and see the result of 
# code_warntype for the call represented by that bar.
;