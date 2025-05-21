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
    @test real(eiginc.p[1:nmod])    ≈ [0.10517722364770647 , 0.1948227763522931 , 0.14999999999999958 , 0.1499999999999998 , 0.14999999999999933 ] 
    @test imag(eiginc.p[1:nmod])   ≈  [- 4.762644452370236e-17, 9.598886497048079e-17, - 0.3995436605528897, 0.39954366055289037, 0.6862471569706331] 
    @test state.X[1][1:2:end] ≈ [-0.295349913218206, -0.8641449877271711, -1.3688503816966837, -1.7720344049295693, -2.043794742866822, -2.163976174464458] atol=1e-6
    @test state.X[2][1:2:end] ≈ [-0.058639602696423314, -0.17156977701560847, -0.27177540585304977, -0.35182469612812317, -0.40578075818226267, -0.4296419177254722] atol=1e-6
    @test state.X[3][1:2:end] ≈ [-0.011539889743192326, -0.033763808398807195, -0.05348361984462822, -0.06923679587785343, -0.07985499549808367, -0.08455071541451746] atol=1e-6
    end


end