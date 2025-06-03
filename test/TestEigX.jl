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

@testset "EigX" begin
    @test eiginc.ω[1:nmod]    ≈ [ 0.143146493914704,    0.42677293340604844,    0.7024494006334382,    0.9650405646123343,    1.209654848799371]
    @test state.X[1][1:2:end] ≈ [-0.7213643765727504, -1.737944612133856, -1.7442767943801336, -0.7907521223050666, 0.4737759750318221, 1.2171685211838659]
    @test state.X[2][1:2:end] ≈ [0.3545303894837429, 0.8385020529010675, 0.7901140341220919, 0.24008758208134098, -0.4623689665273218, -0.8712703649413108]
    @test state.X[3][1:2:end] ≈ [0.12220676482804212, 0.289684442982145, 0.2751522654759765, 0.08895156129160814, -0.14980934227064358, -0.28894230480333694]
    end

eiginc          = solve(EigX{ℂ};state=initialstate,nmod,verbose=false)
imod            = [1,    2]
A               = [1,4+5im] 
state           = increment(initialstate,eiginc,imod,A)

@testset "EigX" begin
    @test abs.(eiginc.p[1:nmod])    ≈ abs.([0.10517722364770647 - 4.762644452370236e-17im, 0.1948227763522931 + 9.598886497048079e-17im, 0.14999999999999958 - 0.3995436605528897im, 0.1499999999999998 + 0.39954366055289037im, 0.14999999999999933 + 0.6862471569706331im])
    @test state.X[1][1:2:end] ≈ [-0.22624258166643885, -0.6619483677080589, -1.0485604715290562, -1.357405642021717, -1.565578245763244, -1.6576390730604207] atol=1e-6
    @test state.X[2][1:2:end] ≈ [-0.04892234955705188, -0.1431387018093662, -0.22673911136446892, -0.2935233182920837, -0.33853824348180706, -0.3584453358618972] atol=1e-6
    @test state.X[3][1:2:end] ≈  [-0.010040786515169252, -0.029377680343219508, -0.04653576601415812, -0.06024250680679059, -0.06948133646069415, -0.07356705324527212] atol=1e-6
end


end