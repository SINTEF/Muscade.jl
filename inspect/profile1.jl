using Profile,ProfileView,BenchmarkTools
using Muscade

include("../Test/SomeElements.jl")

model           = Model(:TestModel)
n1              = addnode!(model,𝕣[ 0, 0])  # moving node
n2              = addnode!(model,𝕣[10, 0])  # anchor 1 
n3              = addnode!(model,𝕣[ 0,10])  # anchor 2 
n4              = addnode!(model,𝕣[     ])  # A-nod for springs
e1              = addelement!(model,Spring{2},[n1,n2,n4], EI=1)
e2              = addelement!(model,Spring{2},[n1,n3,n4], EI=1)
@once f1(t)     = t
e3              = addelement!(model,DofLoad  ,[n1], field=:tx1      ,value=f1)
e4              = addelement!(model,Hold  ,[n2], field=:tx1)
e5              = addelement!(model,Hold  ,[n2], field=:tx2)
e6              = addelement!(model,Hold  ,[n3], field=:tx1)
e7              = addelement!(model,Hold  ,[n3], field=:tx2)
@once f2(x,t)   = 1x^2
@once f3(a)     = 0.11a^2
e8              = addelement!(model,SingleDofCost ,[n1], class=:X ,field=:tx1      ,cost=f2)
e9              = addelement!(model,SingleDofCost ,[n1], class=:X ,field=:tx2      ,cost=f2)
e10             = addelement!(model,SingleDofCost ,[n4], class=:A ,field=:ΞL₀      ,cost=f3)
e10             = addelement!(model,SingleDofCost ,[n4], class=:A ,field=:ΞEI      ,cost=f3)
initialstate    = initialize!(model)
@show typeof(initialstate)

stateX          = solve(StaticX;  initialstate,time=[.5,1.],verbose=false)
#state           = solve(StaticXUA;initialstate=stateX,maxYiter= 50,verbose=false)
#@btime state    = solve(StaticXUA;initialstate=stateX,maxYiter= 0,verbose=false) 

Profile.clear()
Profile.@profile for i=1:25000#1000#25000
    local stateX          = solve(StaticX;  initialstate,time=[.5,1.],verbose=false)
#    local state  = solve(StaticXUA;initialstate=stateX,maxYiter= 50,verbose=false);
end
ProfileView.view(fontsize=30);

# After clicking on a bar in the flame diagram, you can type warntype_last() and see the result of 
# code_warntype for the call represented by that bar.
;


