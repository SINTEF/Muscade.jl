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
    @test state.X[1][1:2:end] ≈ [0.7213643765727511,      1.7379446121338569,      1.7442767943801332,      0.7907521223050662,     -0.473775975031821,     -1.2171685211838656]
    @test state.X[2][1:2:end] ≈ [-0.35453038948374327,    -0.8385020529010682,    -0.7901140341220925,    -0.24008758208134254,     0.4623689665273183,     0.8712703649413073]
    @test state.X[3][1:2:end] ≈ [ -0.12220676482804241,    -0.2896844429821456,    -0.27515226547597676,    -0.08895156129160818,     0.14980934227064358,     0.28894230480333727]
end

res             = solve(EigX{ℂ};state=initialstate,nmod,verbose=false)
imod            = [1,    2]
A               = [1,4+5im] 

state           = increment(initialstate,res,imod,A)

@testset "EigX" begin
    @test res.p[1:nmod]       ≈ [ 0.10589751672082914 - 3.3615962156811485e-17im,    0.19410248327916924 - 7.132917772279973e-17im,    0.14999999999999786 - 0.40440236912020056im,     0.1499999999999976 + 0.4044023691202014im,    0.15000000000000172 - 0.6993716619361788im]
    @test state.X[1][1:2:end] ≈ [ -0.47361538562180927,    -0.701490253292311,    -1.1139019374721226,    -1.4434435839966486,    -1.6655985565824047,    -1.7308909726007045]
    @test state.X[2][1:2:end] ≈ [ -0.08652580617119972,    -0.12815675235648696,    -0.20350112361509923,    -0.2637057907313993,    -0.30429175706922545,    -0.3162201679789166]
    @test state.X[3][1:2:end] ≈ [ -0.016222591349616317,    -0.024027913915759348,    -0.0381541150979112,    -0.04944179625554147,    -0.057051197144856625,    -0.05928763670203735]
end


end