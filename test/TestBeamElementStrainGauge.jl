#module TestBeamElements

using Test, Muscade, StaticArrays, LinearAlgebra
include("../examples/BeamElements.jl")


###
L               = 5
model           = Model(:TestModel)
node1           = addnode!(model,𝕣[0,0,0])
node2           = addnode!(model,𝕣[4,0,0])
elnod           = [model.nod[n.inod] for n∈[node1,node2]]
mat             = BeamCrossSection(EA=10.,EI₂=3.,EI₃=3.,GJ=4.,μ=1.,ι₁=1.0)
P               = SMatrix{3,5}(0.,.5,0.,  0.,0,.5,   0.,-.5,0.,  0.,0,-.5,  0.,.5,0.   )
D               = SMatrix{3,5}(1.,0.,0.,  1.,0.,0.,  1.,0.,0.,   1.,0.,0.,  1/√2,0,1/√2)
L               = 0.1
σ               = 1e-5
εₘ(t)           = SVector(1e-4*cos(t),1e-4*sin(t),-1e-4*cos(t),-1e-4*sin(t),5e-5*cos(t))
instrumentedbeam = StrainGaugeOnEulerBeam3D(elnod;P,D,L,σ,εₘ,elementkwargs=(mat=mat,orient2=SVector(0.,1.,0.)))
@testset "constructor" begin
    @test instrumentedbeam.eleobj.cₘ    ≈ [2.0, 0.0, 0.0]
    @test instrumentedbeam.eleobj.rₘ    ≈ I
    @test instrumentedbeam.eleobj.ζnod  ≈ [-0.5, 0.5]
    @test instrumentedbeam.eleobj.tgₘ   ≈ [4.0, 0.0, 0.0]
    @test instrumentedbeam.eleobj.tgₑ   ≈ [4.0, 0.0, 0.0]
    @test instrumentedbeam.L            ≈ L 
    @test instrumentedbeam.q            ≈ .5*σ^-2
    @test instrumentedbeam.P            ≈ P
    @test instrumentedbeam.D            ≈ D
    @test instrumentedbeam.E            ≈ SVector(1.,1.,1.,1.,.5) 
    @test instrumentedbeam.K1           ≈ SVector(0.,0.,0.,0.,.25) 
    @test instrumentedbeam.K2           ≈ SVector(-0.5, 0, 0.5, 0,-.25)
    @test instrumentedbeam.K3           ≈ SVector(0, -0.5, 0, 0.5,0) 
    @test typeof(instrumentedbeam)      ==  StrainGaugeOnEulerBeam3D{5, EulerBeam3D{BeamCrossSection, false}, typeof(εₘ), @NamedTuple{strain::@NamedTuple{ε::Nothing, κ::Nothing}}}
end

Λ   = SVector(0.,0,0,0,0,0,0,0,0,0,0,0)
X   = (SVector(0,0,0,0,0.1,0.0, 0,0,0,0,-0.1,-0.0),)
U   = (SVector{0,𝕣}(),)
A   = SVector{0,𝕣}()

X   = (SVector(0,0,0, 0,.1,0, 0,0,0, 0,-.1,0),)
L,FB,eleres = Muscade.lagrangian(instrumentedbeam,Λ,X,U,A,0.,nothing,(;),@request (εₐₓ,κ,ε))
@testset "strain1" begin
    @test eleres.εₐₓ  ≈ 0.               atol = 1e-12
    @test eleres.κ*20 ≈ [0, 0, 1]        atol = 1e-12
    @test eleres.ε*40 ≈ [0, -1, 0, 1,0]  atol = 1e-12
end
X   = (SVector(0,0,0, 0,0,.1, 0,0,0, 0,0,-.1),)
L,FB,eleres = Muscade.lagrangian(instrumentedbeam,Λ,X,U,A,0.,nothing,(;),@request (εₐₓ,κ,ε))
@testset "strain2" begin
    @test eleres.εₐₓ  ≈ 0.               atol = 1e-12
    @test eleres.κ*20 ≈ [0,-1, 0]        atol = 1e-12
    @test eleres.ε*40 ≈ [1, 0,-1, 0,.5]  atol = 1e-12
end
X   = (SVector(0,0,0, 0,0,0, 0,0,0, 1.,0,0),)
L,FB,eleres = Muscade.lagrangian(instrumentedbeam,Λ,X,U,A,0.,nothing,(;),@request (εₐₓ,κ,ε))
@testset "strain3" begin
    @test eleres.εₐₓ  ≈ 0.               atol = 1e-12
    @test eleres.κ ≈ [.25,0, 0]          atol = 1e-12
    @test eleres.ε ≈ [0, 0,0, 0,0.0625]  atol = 1e-12
end


X   = (SVector(0,0,0, 0,.1,0, 0,0,0, 0,-.1,0),)
out = diffed_lagrangian(instrumentedbeam;Λ,X,U,A,t=0.,SP=nothing)
i    = SVector(1,4,5,6)
@testset "strain cost" begin
    @test out.HL[2,2][1,1][i,i] ≈ [  2.65625e9    8.20345e7    4.33681e-9  -7.42155e7;
                                     8.20345e7    4.38524e7    2.27692e-9  -5.46027e7;
                                     4.33681e-9   2.27692e-9   3.125e8     -2.05989e-9;
                                    -7.42155e7   -5.46027e7   -2.05989e-9   3.47751e8] rtol = 1e-4
    @test out.∇L[2][1][i] ≈ [       62500.00000117949,   45324.23048065669,       6.2500000000000015e7, -279686.1966132973]                                
end



#@testset "strain values" begin
