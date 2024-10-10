#module TestOneDofViscous
using StaticArrays,Test, Muscade, LinearAlgebra,Interpolations,GLMakie

struct SingleDecayAcost{Field,Tcost,Tcostargs} <: AbstractElement
    cost     :: Tcost     
    costargs :: Tcostargs
    fac      :: 𝕣1
end
SingleDecayAcost(nod::Vector{Node};field::Symbol,fac,cost::Function ,costargs=()) = SingleDecayAcost{field,typeof(cost),typeof(costargs)}(cost,costargs,fac)
Muscade.doflist(::Type{<:SingleDecayAcost{Field,Tcost,Tcostargs}}) where{Field,Tcost,Tcostargs} = (inod=(1,),class=(:A,),field=(Field,))
@espy function Muscade.lagrangian(o::SingleDecayAcost,Λ,X,U,A,t,SP,dbg) 
    iter  = min(length(o.fac),default{:iter}(SP,length(o.fac)))
    ☼cost = o.cost(    A[1]  ,o.costargs...)
    return cost*o.fac[iter],noFB
end
struct OneDofViscous <: AbstractElement
    K   :: SMatrix{1,1,𝕣}
    C   :: SMatrix{1,1,𝕣}
end
OneDofViscous(nod::Vector{Node};K,C)  = OneDofViscous(K,C)

Muscade.doflist(::Type{<:OneDofViscous}) = (  inod = (1,1,1),                                             
                                                 class= (:X,:U,:A), 
                                                 field= (:surge,:surge,:C))

@espy function Muscade.residual(o::OneDofViscous,   X,U,A,t,SP,dbg) 
    x,x′    = ∂0(X),∂1(X)   
    ☼u         = ∂0(U)
    a          = exp10.(A)
    ☼r₁        = (o.C.*a[@SVector [i for i∈1:1]])∘x′
    ☼r₀        = o.K∘x
    return r₀+r₁-u,  noFB
end

K         = SMatrix{1,1}(1.0)
C         = SMatrix{1,1}(1.0)

#Solve direct problem
model     = Model(:OneDofViscous)
n1        = addnode!(model,𝕣[0,0,0])  
e1        = addelement!(model,OneDofViscous,[n1]; K,C)
initialstate    = initialize!(model;time=0.)
initialstate    = setdof!(initialstate,[1.0];   field=:surge,   nodID=[n1], order=0)                                          
T               = 0.2 *(1:6)
#T               = 0.01 *(1:200)
state           = solve(SweepX{2};  initialstate,time= T,verbose=false)
surge   = [s.X[1][1] for s∈state]

# Create fake measurements
surgeMeas = surge   + .0 * randn(length(T))

# Create intial guesses for M and C
maxDevToModel = 0.0; 
devToModelC = 1. .+ maxDevToModel*((rand(1).-.5)*2);
Cguess = C .* devToModelC[@SVector [i for i∈1:1 ]]

# Create XUA model
modelXUA     = Model(:OneDofViscous)
n1        = addnode!(modelXUA,𝕣[0,0,0])  
e1        = addelement!(modelXUA,OneDofViscous,[n1]; K,C=Cguess)

# Assign costs to unknown forces
Quu       = @SVector [0.1 ^-2 for i=1:1]  
e2        = addelement!(modelXUA,SingleDofCost     ,[n1]; class=:U,field=:surge,           cost=(u,t)-> 0.5*Quu[1]*u^2)  
# Assign costs to variations of model parameters (wrt guess).
fac       = [100,90,80,70,60,50,40,30,20,10,1] 
QCaa      = @SVector [.1 ^-2 for i=1:1]  
e3        = addelement!(modelXUA,SingleDecayAcost  ,[n1];          field=:C,fac,           cost=(a  )-> 0.5*QCaa[1]/length(T)*a^2) 
# Assign costs to measurement errors
surgeInt    = linear_interpolation(T, surgeMeas)
@once devSurge(surge,t)     = 5e-2 ^-2 * (surge-surgeInt(t))^2
#what if we do not use Interplation? --> no difference
# @once devSurge(surge,t) = 5e-2 ^-2 *(surge-exp(-t))^2
e5             = addelement!(modelXUA,SingleDofCost,[n1];class=:X,field=:surge,    cost=devSurge)

# #Setting scale to improve convergence
# myScaling = (   X=(surge=1.,    sway=1.,    yaw=1.),
#                 U=(surge=1e-1,  sway=1e-1,  yaw=1e-1),
#                 A=( C11=1e-1,C12=1e-1,C16=1e-1,C22=1e-1,C26=1e-1,C66=1e-1,
#                     M11=1e-1,M12=1e-1,M16=1e-1,M22=1e-1,M26=1e-1,M66=1e-1))
# setscale!(modelXUA;scale=myScaling,Λscale=1e2)

#Solve inverse problem
initialstateXUA    = initialize!(modelXUA;time=0.)
stateXUA         = solve(DirectXUA{1,0,0};initialstate=initialstateXUA,time=T,
                        maxiter=1,saveiter=true,fastresidual=true,
                        maxΔx=Inf,maxΔλ=Inf,maxΔu=Inf,maxΔa=Inf)
niter=findlastassigned(stateXUA)

# Fetch and display estimated model parameters
Cest      = Cguess .* exp10.(SVector{1}(stateXUA[niter][1].A[1])) 
@show Cguess ./ C
@show Cest ./ C;

# Fetch response and loads 
surgeRec   = [s.X[1][1] for s∈stateXUA[niter]]  # cf. getdof
surgeExtF   = [s.U[1][1] for s∈stateXUA[niter]]
req = @request r₂,r₁,r₀  
loads = getresult(stateXUA[niter],req,[e1])
dampingLoads = [loads[i].r₁ for i∈1:length(T)]
stiffnessLoads = [loads[i].r₀ for i∈1:length(T)]

# Display response
fig      = Figure(size = (2000,1000))
display(fig) # open interactive window (gets closed down by "save")
ax=Axis(fig[1,1], ylabel="Surge [m]",        yminorgridvisible = true,xminorgridvisible = true)
lines!(fig[1,1],T,surge,                          label="Exact direct solution")
scatter!(fig[1,1],T,surgeMeas,                          label="Simulated measurements")
lines!(fig[1,1],T,surgeRec,                       label="Inverse solution")
axislegend()


# Display loads 
ax=Axis(fig[1,2], ylabel="Surge force [N]",        yminorgridvisible = true,xminorgridvisible = true)
lines!(fig[1,2],T,-[dampingLoads[i][1] for i∈1:length(T)], label="Damping")
lines!(fig[1,2],T,-[stiffnessLoads[i][1] for i∈1:length(T)], label="Stiffness")
lines!(fig[1,2],T,surgeExtF,                       label="Unknown")
axislegend()


display(fig) # open interactive window (gets closed down by "save")

#end
;