#module TestSweepXA2
using Test
using Muscade


include("SomeElements.jl")
include("../examples/DryFriction.jl")


K,C,M           = .3,1.,1.

model           = Model(:TestModel)
node            = addnode!(model,𝕣[])
osc             = addelement!(model,AdjustableSdofOscillator,[node]; K,C,M)

@functor with(     ) acost(a,σ)=(a/σ)^2
@functor with(σ=.1) xcost(x,t)=(x/σ)^2
cC              = addelement!(model,SingleAcost  ,[node];field=:ΞK,costargs=(.2,),cost=acost)
cM              = addelement!(model,SingleAcost  ,[node];field=:ΞC,costargs=(.1,),cost=acost)
cX              = addelement!(model,SingleDofCost,[node];class=:X ,field=:tx1    ,cost=xcost)

initialstate    = initialize!(model;time=0.)
x,x′,x″         = 0.,1.,0.   
initialstate    = setdof!(initialstate,[x′];field=:tx1,nodID=[node],order=1)  # initial speed




# # Test assembly process
# Δt              = 0.4                          
# model,dis       = initialstate.model,initialstate.dis
# out,asm,Xdofgr  = prepare(Muscade.AssemblySweepXA{2},model,dis)  
# β,γ             = 1/4,1/2
# out.c= (a₁=γ/(β*Δt), a₂=γ/β, a₃=(γ/2β-1)*Δt, b₁=1/(β*Δt^2), b₂=1/(β*Δt), b₃=1/2β) 
# state           = Muscade.State{3,3,1}(copy(initialstate,SP=(γ=0.,))) 
# Muscade.assemble!{:newmark}(out,asm,dis,model,state,(;))
# a                 = out.c.a₂*x′ + out.c.a₃*x″
# b                 = out.c.b₂*x′ + out.c.b₃*x″
# R                 = K*x+C*x′+M*x″ 
# @testset "assemble!" begin
#     @test out.Lλ    ≈ [R-C*a-M*b]
#     @test out.Lx    ≈ [0.]   
#     @test out.Lr[]  ≈ 0.0
#     @test out.La    ≈ [0.,0.]
#     @test out.Lλx   ≈ [K + out.c.a₁*C + out.c.b₁*M ;] 
#     @test out.Lλa   ≈ [2.302585092994046 0.0;]
#     @test out.Lxx   ≈ [200.;] 
#     @test out.Lxr   ≈ [0.]
#     @test out.Lrr[] ≈ 0.0
#     @test out.Lax   ≈ [0.;0.]
#     @test out.Laa   ≈ [50. 0.;0. 200.]
# end

t               = 2.:1:21
state           = solve(SweepXA{2};  initialstate,time= t,verbose=false,catcherror=true)


;
# X  = getdof(state;field=:tx1,nodID=[node],order=0)
# X′ = getdof(state;field=:tx1,nodID=[node],order=1)
# X″ = getdof(state;field=:tx1,nodID=[node],order=2)

# @testset "SDof oscillator output" begin
#     @test X[ 1:10:end] ≈ [0.3653456491624315, 0.039495394592936224, -0.9800856952523974, 0.0204208015740059, 0.10202734361080085, -0.08876358375431506, 0.020279004568508258, 0.009410960025011333, -0.01164572216979256, 0.0044937076208295505]
#     @test X′[1:10:end] ≈ [0.8267282458121575, -0.5543584424772612, 0.1761431693326722, 0.17142199935793515, -0.07990836783415763, 0.009615311853025504, 0.01797222270368847, -0.012951778973228531, 0.0032459856194525815, 0.0015257230558188115]
#     @test X″[1:10:end] ≈ [-0.8663587709392137, -0.033410494488486674, 0.1512397675677002, -0.08357963580025408, -0.012670847567464316, 0.025533223926927903, -0.0130068667401185, 0.0010595839803349055, 0.0027793256172318247, -0.002010048120257385]
# end

# # using GLMakie
# # fig      = Figure(size = (2000,1500))
# # axe      = Axis(fig[1,1],title="Test",xlabel="time",ylabel="x")
# # oedge    = lines!(  axe,T,X , linewidth = 1)
# # oedge    = lines!(  axe,T,X′, linewidth = 1)
# # oedge    = lines!(  axe,T,X″, linewidth = 1)
# #save("C:\\Users\\philippem\\.julia\\dev\\Muscade\\test\\testDynamic.jpg",fig)

# #######

# model           = Model(:TestModel)
# node            = addnode!(model,𝕣[])
# ele             = addelement!(model,SdofOscillator,[node], K₁=1.,K₂=.3,C₁=1.,C₂=2.,M₁=3.)
# drag            = addelement!(model,DryFriction,[node], fieldx=:tx1,friction=0.1)
# initialstate    = Muscade.State{1,3,1}(copy(initialize!(model),time=0.,SP=nothing))  # recast to force the state to have 2nd derivatives, 
# initialstate.X[2][1] = 1.                                                   # so we can set initial velocity
# T               = 0.15 *(1:100)
# state           = solve(SweepX{2};  initialstate,time= T,verbose=false,catcherror=true)

# X  = getdof(state;field=:tx1,nodID=[node],order=0)
# X′ = getdof(state;field=:tx1,nodID=[node],order=1)
# X″ = getdof(state;field=:tx1,nodID=[node],order=2)

# @testset "DryFriction oscillator output, Δx=0" begin
#     @test X[ 1:10:end] ≈ [0.14456443480393982,  0.7975704237766761,  0.5807448605308462,  0.025592798296295478, -0.4191499109805254, -0.47980495035837123, -0.327180958763213, -0.14250971878268207, -0.02302664416796214,  0.0049513031770612665]
#     @test X′[1:10:end] ≈ [  0.9275257973858644,  0.06890381030669016,   -0.3051199796365325,   -0.38287482503031733,   -0.17497655979555798,    0.05871470543387,    0.12727994128778658,    0.10817353690164826,    0.04828950336575319,    0.0]
#     @test X″[1:10:end] ≈ [ -0.9663227015218087,    -0.38893509283263966,    -0.15433355292283862,     0.05463308919551511,     0.19339562928015905,     0.08171052489330786,     0.011795478883263008,    -0.0317198507523637,    -0.043361893119507146,    -0.029111361795204654]
# end

# # ########

# model           = Model(:TestModel)
# node            = addnode!(model,𝕣[])
# ele             = addelement!(model,SdofOscillator,[node], K₁=1.,K₂=.3,C₁=1.,C₂=2.,M₁=3.)
# drag            = addelement!(model,DryFriction,[node], fieldx=:tx1,friction=0.1,Δx=0.3)
# initialstate    = Muscade.State{1,3,1}(copy(initialize!(model),time=0.,SP=nothing))  # recast to force the state to have 2nd derivatives, 
# initialstate.X[2][1] = 1.                                           # so we can set initial velocity
# T               = 0.15 *(1:100)
# state           = solve(SweepX{2};  initialstate,time= T,verbose=false,catcherror=true)

# X  = [s.X[1][1] for s∈state]
# X′ = [s.X[2][1] for s∈state]
# X″ = [s.X[3][1] for s∈state]

# @testset "DryFriction oscillator output, Δx=0.3" begin
#     @test X[ 1:10:end] ≈ [0.14469296858597538,    0.8008980273267399,    0.5415497021192464,   -0.09635105317563947,   -0.5596719823602733,   -0.5665780815035283,   -0.28989109303335375,   -0.028124724913396883,    0.07731217827508936,    0.04100891968258263]
#     @test X′[1:10:end] ≈ [0.9292395811463382,    0.0693102861285251,   -0.362977355522562,   -0.42183677927675134,   -0.16079100719211595,    0.12680704937085094,    0.20617807686338632,    0.1271272660824127,    0.016424801988853523,   -0.05507494559856239]
#     @test X″[1:10:end] ≈ [-0.9434722513821608,    -0.3907491466781534,    -0.17587190387027593,     0.08650341365728759,     0.22492856089769642,     0.1294376784731921,    -0.013913519287633938,    -0.07718750825920974,    -0.06535656010969325,    -0.026317647770282403]
# end

# # using GLMakie
# # fig      = Figure(size = (2000,1500))
# # axe      = Axis(fig[1,1],title="Test",xlabel="time",ylabel="x")
# # oedge    = lines!(  axe,T,X , linewidth = 1)
# # oedge    = lines!(  axe,T,X′, linewidth = 1)
# # oedge    = lines!(  axe,T,X″, linewidth = 1)
# # save("C:\\Users\\philippem\\.julia\\dev\\Muscade\\test\\testDynamicWithFriction.jpg",fig)

#end