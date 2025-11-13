module TestSweepX1
using Test
using Muscade

include("SomeElements.jl")


K,C,M           = 1.,.3,0.

model           = Model(:TestModel)
node            = addnode!(model,𝕣[])
osc             = addelement!(model,AdjustableSdofOscillator,[node]; K,C,M)

initialstate    = initialize!(model;time=0.)
initialstate    = setdof!(initialstate,[1.];field=:tx1,nodID=[node],order=0)  # initial speed

t               = .1:.1:2
state           = solve(SweepX{1};  initialstate,time=t,maxiter= 5,verbose=false,catcherror=true)

x = [state[i].X[1][1] for i = 1:length(t)]
x′= [state[i].X[2][1] for i = 1:length(t)]

@testset "Exponential decay" begin
    @test x' ≈ [0.857143  0.612245  0.437318  0.31237  0.223121  0.159372  0.113837  0.0813124  0.0580803  0.0414859  0.0296328  0.0211663  0.0151188  0.0107991  0.00771366  0.00550976  0.00393554  0.0028111  0.00200793  0.00143424] rtol=1e-5
    @test x′' ≈ [ -2.85714  -2.04082  -1.45773  -1.04123  -0.743738  -0.531241  -0.379458  -0.271041  -0.193601  -0.138286  -0.098776  -0.0705543  -0.0503959  -0.0359971  -0.0257122  -0.0183659  -0.0131185  -0.00937034  -0.0066931  -0.00478079] rtol = 1e-5
end

# using GLMakie
# fig      = Figure(size = (600,450))
# display(fig) # open interactive window (gets closed down by "save")
# axe      = Axis(fig[1,1],title="Test",xlabel="t")
# ox    = lines!(  axe,t,x    ,color = :black, linewidth = 1.)
# ox′    = lines!(  axe,t,x′    ,color = :red, linewidth = 1.)
end