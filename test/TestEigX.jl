#module TestEigX
using Muscade, StaticArrays

include("SomeElements.jl")
using KrylovKit: KrylovKit
# Create the model 
nel         = 10
nnod        = nel+1  
L           = 10 
EI          = 1
M           = 10/nnod
nodeCoord   = reshape(range(0,L,length=nnod),(nnod,1))
model       = Model()
xnod        = addnode!(model,nodeCoord)
anod        = addnode!(model,𝕣[])
mesh_spring = hcat(xnod[1:nnod-1],xnod[2:nnod],[anod for i=1:nel])
mesh_mass   = reshape(xnod,(nnod,1))
spring      = addelement!(model,Spring{1},mesh_spring;EI)
mass        = addelement!(model,SdofOscillator,mesh_mass;M₁=M)
hold        = addelement!(model,SdofOscillator,[xnod[1]];K₁=1.)
#hold        = addelement!(model,Hold,[xnod[1]];field=:tx1)
# Solve the problem 
initialstate    = Muscade.initialize!(model);

state           = solve(EigX     ;initialstate,
                                maxiter     = 200,    # number of 
                                krylovdim   = 100,#nnod÷5,
                                tol         = 1e-3,
                                verbosity   = 3,
                                orth        = KrylovKit.ModifiedGramSchmidtIR(),
                                issymmetric = true,
                                isposdef    = true)

#end