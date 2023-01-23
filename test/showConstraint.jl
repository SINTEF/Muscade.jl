using Muscade,GLMakie

f1(x)           = .15*x^2-.2x-1 +2
f2(x)           = .05*x^2+x +2 +2

g1(x,t)         = x[2]-f1(x[1]) # start outside feasible domain 😃
g2(x,t)         = x[2]-f2(x[1])

gravity(t)      = -2.
gravity2(t)     = 0
model           = Model(:TestModel)
n1              = addnode!(model,𝕣[0,0]) 
e1              = addelement!(model,Constraint,[n1],xinod=(1,1),xfield=(:t1,:t2),λinod=1, λclass=:X, λfield=:λ1,g=g1,mode=inequal)
e2              = addelement!(model,Constraint,[n1],xinod=(1,1),xfield=(:t1,:t2),λinod=1, λclass=:X, λfield=:λ2,g=g2,mode=inequal)
e3              = addelement!(model,DofLoad   ,[n1],field=:t1,value=gravity2)
e4              = addelement!(model,DofLoad   ,[n1],field=:t2,value=gravity)
e5              = addelement!(model,QuickFix  ,[n1],inod=(1,1),field=(:t1,:t2),res=(x,u,a,t)->0.1x)

state           = solve(staticX;model,time=[0.],verbose=false,γ0=.5,γfac=.25,saveiter=true,maxiter=1000) # because there is zero physical stiffness in this model, setting γ0=0 gives singularity if one or more constraint is inactive
last            = findlastassigned(state)
@show last
X               = state[last].X[1][1:2] 
x1              = [0,[state[i].X[1][1] for i ∈1:last]...]
x2              = [0,[state[i].X[1][2] for i ∈1:last]...]

fig      = Figure(resolution = (800,800))
display(fig) # open interactive window (gets closed down by "save")
axe = Axis(fig[1,1],title="Modified interior point")
a1 = -4:.01:0
a2 = -4:.01:0
b = f1.(a1)
c = f2.(a2)
scatter!(axe,x1,x2)
lines!(  axe,x1,x2,color = :black, linewidth = 1)
lines!(  axe,a1,b ,color = :red, linewidth = 1)
lines!(  axe,a2,c ,color = :blue, linewidth = 1)

;