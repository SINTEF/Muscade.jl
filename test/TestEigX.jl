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
res             = solve(EigX{ℝ};state=initialstate,nmod,verbose=false)
imod            = [1,    2]
A               = [1,4+5im] 

state           = increment(initialstate,res,imod,A)

@testset "EigX" begin
    @test res.ω[1:nmod]       ≈ [ 0.143146493914704,    0.42677293340604844,    0.7024494006334382,    0.9650405646123343,    1.209654848799371]
    @test state.X[1][1:2:end] ≈ [-0.6077935403710729, -1.405655132897898, -1.2179130626616075, -0.10935212813168926, 1.259675911289766, 2.0492817965071035]
    @test state.X[2][1:2:end] ≈ [0.3545303894837429, 0.8385020529010675, 0.7901140341220919, 0.24008758208134098, -0.4623689665273218, -0.8712703649413108]
    @test state.X[3][1:2:end] ≈ [0.11987959405446318, 0.28287552627158646, 0.2643665890321386, 0.07498904939514335, -0.1659131539866164, -0.30599307029388034]
end

res             = solve(EigX{ℂ};state=initialstate,nmod,verbose=false)
imod            = [1,    2]
A               = [1,4+5im] 

state           = increment(initialstate,res,imod,A)

@testset "EigX" begin
    @test res.p[1:nmod]       ≈ [0.10517722364770678 + 1.7075171462914056e-18im, 0.19482277635229214 + 1.20559891185176e-16im, 0.14999999999999986 + 0.39954366055289015im, 0.14999999999999986 - 0.39954366055289053im, 0.14999999999999922 + 0.686247156970632im]
    @test state.X[1][1:2:end] ≈ [-0.3049538176434667, -0.8922444233466716, -1.413361341916891, -1.829655715455044, -2.1102528946958445, -2.234342270501683]
    @test state.X[2][1:2:end] ≈ [-0.06416546640074926, -0.18773754009678098, -0.29738597928593286, -0.38497866084177756, -0.4440192362831729, -0.4701289362222911]
    @test state.X[3][1:2:end] ≈ [-0.013000856029516719, -0.038038354071799764, -0.06025472140646839, -0.07800227170142804, -0.08996475034094181, -0.09525495500901532]
    end


end