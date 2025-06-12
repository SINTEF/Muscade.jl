module TestEigX
using Muscade, StaticArrays,Test,Random

include("SomeElements.jl")
# Create the model 
nel         = 10
nnod        = nel+1  
L           = 10 
EI          = 1
M           = 10/nnod
C           = 3/nnod
nodeCoord   = reshape(range(0,L,length=nnod),(nnod,1))
model       = Model()
xnod        = addnode!(model,nodeCoord)
anod        = addnode!(model,𝕣[])
mesh_spring = hcat(xnod[1:nnod-1],xnod[2:nnod],[anod for i=1:nel])
mesh_mass   = reshape(xnod,(nnod,1))
spring      = addelement!(model,Spring{1},mesh_spring;EI)
mass        = addelement!(model,SdofOscillator,mesh_mass;M₁=M,C₁=C)
hold        = addelement!(model,SdofOscillator,[xnod[1]];K₁=1.)
#hold        = addelement!(model,Hold,[xnod[1]];field=:tx1)

# Solve the problem 
initialstate    = initialize!(model);
nmod            = 5
Random.seed!(1234) # eigensolvers use a random start.  Freeze it for reproducible tests

eiginc          = solve(EigX{ℝ};state=initialstate,nmod,verbose=false)
imod            = [1,    2]
A               = [1,4+5im] 
state           = increment(initialstate,eiginc,imod,A)

@testset "EigXR" begin
    @test eiginc.ω[1:nmod]    ≈ [ 0.143146493914704,    0.42677293340604844,    0.7024494006334382,    0.9650405646123343,    1.209654848799371]
    @test state.X[1][1:2:end] ≈ [1.7338129265935287, 4.177182755457795, 4.192402275285353, 1.9005876862087123, -1.1387295193620908, -2.9254875682906434]
    @test state.X[2][1:2:end] ≈ [-0.8521204984886861, -2.015355547777684, -1.8990540291889357, -0.5770550457522317, 1.1113125586121952, 2.0941148057104684]
    @test state.X[3][1:2:end] ≈ [-0.2937262712953869, -0.6962620392508903, -0.6613336756803813, -0.2137967604399907, 0.36006958839215825, 0.694478296097603]
end

eiginc          = solve(EigX{ℂ};state=initialstate,nmod,verbose=false)
imod            = [1,    2]
A               = [1,4+5im] 
state           = increment(initialstate,eiginc,imod,A)

@testset "EigXC" begin
    @test abs.(eiginc.p[1:nmod])    ≈ abs.([0.10517722364770647 - 4.762644452370236e-17im, 0.1948227763522931 + 9.598886497048079e-17im, 0.14999999999999958 - 0.3995436605528897im, 0.1499999999999998 + 0.39954366055289037im, 0.14999999999999933 + 0.6862471569706331im])
    @test state.X[1][1:2:end] ≈  [0.6824241336467091, 1.9966601248302318, 3.162812968667382, 4.094394443525039, 4.722313413114684, 5.0] atol=1e-6
    @test state.X[2][1:2:end] ≈  [0.12071650663891864, 0.35319652886061903, 0.5594815803033334, 0.7242724423928726, 0.8353473307560867, 0.8844683290568793] atol=1e-6
    @test state.X[3][1:2:end] ≈  [0.02223145453650511, 0.06504555832867998, 0.10303553062324605, 0.13338382896767445, 0.15383965890798573, 0.1628859051167022] atol=1e-6
end

end