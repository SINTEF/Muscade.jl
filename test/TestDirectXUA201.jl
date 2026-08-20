module TestDirectXUA201

using Test
using Muscade
include("SomeElements.jl")

# Δt              = 0.01
# t               = Δt:Δt:1000*Δt
Δt              = 0.1
t               = Δt:Δt:100*Δt
# Δt              = 1.
# t               = Δt:Δt:10*Δt

### Target

K,C,M           = 1.,.4,.3

model           = Model(:Target)
node            = addnode!(model,𝕣[])
osc             = addelement!(model,AdjustableSdofOscillator,[node]; K,C,M)

initialstate    = initialize!(model;time=0.)
x,x′,x″         = 0.,1.,0.   
initialstate    = setdof!(initialstate,[x′];field=:tx1,nodID=[node],order=1)  # initial speed

state0          = solve(SweepX{2};  initialstate,time= t,verbose=false,catcherror=true)
xtarget         = getdof(state0;field=:tx1,nodID=[node])
target          = Muscade.FunctionFromVector(t,xtarget)

### Optimize

K,C,M           = 1.,.6,.3

model           = Model(:Optim)
node            = addnode!(model,𝕣[])
osc             = addelement!(model,AdjustableSdofOscillator,[node]; K,C,M)

@functor with(σ=1e4        ) acost(a  )=( a           /σ)^2
@functor with(σ=1e-3       ) ucost(u,t  )=( u           /σ)^2
@functor with(σ=1e-4,target) xcost(x,t)=((x-target(t))/σ)^2

cX              = addelement!(model,SingleDofCost,[node];class=:X ,field=:tx1    ,cost=xcost)
cU              = addelement!(model,SingleDofCost,[node];class=:U ,field=:tu1    ,cost=ucost)
cC              = addelement!(model,SingleAcost  ,[node];field=:ΞC,               cost=acost)

initialstate    = initialize!(model;time=0.)
x,x′,x″         = 0.,1.,0.   
initialstate2    = setdof!(initialstate,[x′];field=:tx1,nodID=[node],order=1)  # initial speed

statewrong      = solve(SweepX{2};  initialstate=initialstate2,time= t,verbose=false,catcherror=true)
state           = solve(DirectXUA{2,0,1};  primerstate=[initialstate],time= [t],verbose=false,catcherror=true,maxiter=50,maxΔa=1e-5,maxΔx=1e-5,maxΔu=1e-5,maxΔλ=Muscade.∞)

A              = getdof(state[1][1];class=:A,field=:ΞC,nodID=[node])
xfitted        = getdof(state[1];field=:tx1,nodID=[node])
xwrong         = getdof(statewrong;field=:tx1,nodID=[node])

# # Adjust the damping in the green model to match the red measurements
# using GLMakie
# fig       = Figure(size = (1000,800))
# axeX      = Axis(fig[1,1])
# lines!(  axeX,t,xfitted[1,:],color=:black)
# lines!(  axeX,t,xwrong[1,:],color=:lightgreen)
# scatter!(axeX,t[1:2:end],xtarget[1,1:2:end],color=:red)
# display(fig)

@testset "Retrieve Adof" begin
    @test  A ≈ [log10(0.4)-log10(0.6)] rtol= 1e-2
    @test  xfitted≈xtarget rtol= 1e-3
end

### Trajectory

traj = Vector{Muscade.State}(undef,length(t))
intermediateState = nothing
let primer = initialstate
    for (it,ti) in enumerate(t)
        # # Solve X
        # itermax = 10
        # local stateX = solve(DirectXUA{0,0,0};
        #     initialstate=[primer],
        #     time=[ti:eps():ti+eps()],
        #     verbose=false,
        #     maxiter=itermax,
        #     maxΔx=1e-5,
        #     maxΔu=1e-5,
        #     maxΔa=Inf,
        #     maxΔλ=Inf,
        #     saveiter=true,
        # )
        
        # local laststepX = findlastassigned(stateX)
        # traj[it] = stateX[laststepX][1][end]
        # global intermediateState = stateX[laststepX][1][end]
        # primer = intermediateState
        
        traj[it] = initialstate
    end
end

statetraj = solve(DirectXUA{2,0,1};  primerstate=[traj],time= [t],verbose=false,catcherror=true,maxiter=20,maxΔa=1e-5,maxΔx=1e-5,maxΔu=1e-5,maxΔλ=Muscade.∞)
A              = getdof(statetraj[1][1];class=:A,field=:ΞC,nodID=[node])
xfittedtraj        = getdof(statetraj[1];field=:tx1,nodID=[node])
xtraj        = getdof(traj;field=:tx1,nodID=[node])


@testset "Initial full trajectory" begin
    @test  A ≈ [log10(0.4)-log10(0.6)] rtol= 1e-2
    @test  xfittedtraj≈xtarget rtol= 1e-3
    @test  xfittedtraj≈xfitted rtol= 1e-10
    @test  statetraj[1][1].A === statetraj[1][2].A 
end

end