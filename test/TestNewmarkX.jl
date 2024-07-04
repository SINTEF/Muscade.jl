module TestNewmarkX
using Test
using Muscade

include("SomeElements.jl")

model           = Model(:TestModel)
node            = addnode!(model,𝕣[])
ele             = addelement!(model,SdofOscillator,[node], K₁=1.,K₂=.3,C₁=1.,C₂=2.,M₁=3.)

initialstate    = Muscade.State{1,3,1}(initialize!(model;time=0.))
initialstate.X[2][1] = 1.
state           = solve(NewmarkX;  initialstate,time= 0.4 *(1:100),verbose=false,catcherror=true)#,maxiter=2)#,maxΔx=Muscade.∞)

X  = [s.X[1][1] for s∈state]
X′ = [s.X[2][1] for s∈state]
X″ = [s.X[3][1] for s∈state]

@testset "SDof oscillator output" begin
    @test X[1:10:end] ≈ [0.3653456491624315, 0.039495394592936224,   -0.9800856952523974, 0.0204208015740059, 0.10202734361080085,   -0.08876358375431506, 0.020279004568508258, 0.009410960025011333,   -0.01164572216979256, 0.0044937076208295505]
    @test X′[1:10:end] ≈ [0.8267282458121575, -0.5543584424772612,  0.1761431693326722,  0.17142199935793515, -0.07990836783415763,  0.009615311853025504,  0.01797222270368847, -0.012951778973228531,  0.0032459856194525815,  0.0015257230558188115]
    @test X″[1:10:end] ≈ [-0.8663587709392137, -0.033410494488486674,  0.1512397675677002, -0.08357963580025408, -0.012670847567464316,  0.025533223926927903, -0.0130068667401185,  0.0010595839803349055,  0.0027793256172318247, -0.002010048120257385]
end

# using GLMakie
# fig      = Figure(size = (2000,1500))
# axe      = Axis(fig[1,1],title="Test",xlabel="time",ylabel="x")
# oedge    = lines!(  axe,T,X, linewidth = 1)
# oedge    = lines!(  axe,T,X′, linewidth = 1)
# oedge    = lines!(  axe,T,X″, linewidth = 1)
# save("C:\\Users\\philippem\\.julia\\dev\\Muscade\\test\\testDynamic.jpg",fig)

end