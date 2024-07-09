module TestNewmarkX
using Test
using Muscade

include("SomeElements.jl")

model           = Model(:TestModel)
node            = addnode!(model,𝕣[])
ele             = addelement!(model,SdofOscillator,[node], K₁=1.,K₂=.3,C₁=1.,C₂=2.,M₁=3.)
initialstate    = Muscade.State{1,3,1}(initialize!(model;time=0.))  # recast to force the state to have 2nd derivatives, 
initialstate.X[2][1] = 1.                                           # so we can set initial velocity
T               = 0.4 *(1:100)
state           = solve(NewmarkX;  initialstate,time= T,verbose=false,catcherror=true)

X  = [s.X[1][1] for s∈state]
X′ = [s.X[2][1] for s∈state]
X″ = [s.X[3][1] for s∈state]

@testset "SDof oscillator output" begin
    @test X[ 1:10:end] ≈ [0.3653456491624315, 0.039495394592936224, -0.9800856952523974, 0.0204208015740059, 0.10202734361080085, -0.08876358375431506, 0.020279004568508258, 0.009410960025011333, -0.01164572216979256, 0.0044937076208295505]
    @test X′[1:10:end] ≈ [0.8267282458121575, -0.5543584424772612, 0.1761431693326722, 0.17142199935793515, -0.07990836783415763, 0.009615311853025504, 0.01797222270368847, -0.012951778973228531, 0.0032459856194525815, 0.0015257230558188115]
    @test X″[1:10:end] ≈ [-0.8663587709392137, -0.033410494488486674, 0.1512397675677002, -0.08357963580025408, -0.012670847567464316, 0.025533223926927903, -0.0130068667401185, 0.0010595839803349055, 0.0027793256172318247, -0.002010048120257385]
end

# using GLMakie
# fig      = Figure(size = (2000,1500))
# axe      = Axis(fig[1,1],title="Test",xlabel="time",ylabel="x")
# oedge    = lines!(  axe,T,X , linewidth = 1)
# oedge    = lines!(  axe,T,X′, linewidth = 1)
# oedge    = lines!(  axe,T,X″, linewidth = 1)
#save("C:\\Users\\philippem\\.julia\\dev\\Muscade\\test\\testDynamic.jpg",fig)

#######

model           = Model(:TestModel)
node            = addnode!(model,𝕣[])
ele             = addelement!(model,SdofOscillator,[node], K₁=1.,K₂=.3,C₁=1.,C₂=2.,M₁=3.)
drag            = addelement!(model,DryFriction,[node], drag=0.1,field=:x)
initialstate    = Muscade.State{1,3,1}(initialize!(model;time=0.))  # recast to force the state to have 2nd derivatives, 
initialstate.X[2][1] = 1.                                           # so we can set initial velocity
T               = 0.15 *(1:100)
state           = solve(NewmarkX;  initialstate,time= T,verbose=false,catcherror=true)

X  = [s.X[1][1] for s∈state]
X′ = [s.X[2][1] for s∈state]
X″ = [s.X[3][1] for s∈state]

@testset "DryFriction oscillator output" begin
    @test X[ 1:10:end] ≈ [0.14456443480393982,  0.7975704237766761,  0.5807448605308462,  0.025592798296295478, -0.4191499109805254, -0.47980495035837123, -0.327180958763213, -0.14250971878268207, -0.02302664416796214,  0.0049513031770612665]
    @test X′[1:10:end] ≈ [  0.9275257973858644,  0.06890381030669016,   -0.3051199796365325,   -0.38287482503031733,   -0.17497655979555798,    0.05871470543387,    0.12727994128778658,    0.10817353690164826,    0.04828950336575319,    0.0]
    @test X″[1:10:end] ≈ [ -0.9663227015218087,    -0.38893509283263966,    -0.15433355292283862,     0.05463308919551511,     0.19339562928015905,     0.08171052489330786,     0.011795478883263008,    -0.0317198507523637,    -0.043361893119507146,    -0.029111361795204654]
end

# using GLMakie
# fig      = Figure(size = (2000,1500))
# axe      = Axis(fig[1,1],title="Test",xlabel="time",ylabel="x")
# oedge    = lines!(  axe,T,X , linewidth = 1)
# oedge    = lines!(  axe,T,X′, linewidth = 1)
# oedge    = lines!(  axe,T,X″, linewidth = 1)
# save("C:\\Users\\philippem\\.julia\\dev\\Muscade\\test\\testDynamicWithFriction.jpg",fig)

end