using Profile,ProfileView,BenchmarkTools
using Muscade

model           = Model(:TestModel)
n1              = addnode!(model,ð£[ 0, 0])  # moving node
n2              = addnode!(model,ð£[10, 0])  # anchor 1 
n3              = addnode!(model,ð£[ 0,10])  # anchor 2 
n4              = addnode!(model,ð£[     ])  # A-nod for springs
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
e10             = addelement!(model,AdofCost ,[n4], field=:ÎLâ      ,cost=f3)
e10             = addelement!(model,AdofCost ,[n4], field=:ÎEI      ,cost=f3)
state     = solve(staticXUA;model,time=[.5,1.],verbose=false)
@btime state     = solve(staticXUA;model,time=[.5,1.],verbose=false)

Profile.clear()
Profile.@profile for i=1:1000
    local state  = solve(staticXUA;model,time=[.5,1.],verbose=false);
end
ProfileView.view(fontsize=30);
# After clicking on a bar, you can type warntype_last() and see the result of 
# code_warntype for the call represented by that bar.
;