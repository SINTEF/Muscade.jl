module TestStrainGaugeOnEulerBeam3D

using Test, Muscade, StaticArrays, LinearAlgebra
include("../examples/BeamElement.jl")
include("../examples/StrainGaugeOnBeamElement.jl")

model            = Model(:TestModel)
node1            = addnode!(model,𝕣[0,0,0])
node2            = addnode!(model,𝕣[4,0,0])
elnod            = [model.nod[n.inod] for n∈[node1,node2]]
mat              = BeamCrossSection(EA=10.,EI₂=3.,EI₃=3.,GJ=4.,μ=1.,ι₁=1.0)
P                = SMatrix{3,5}(0.,.5,0.,  0.,0,.5,   0.,-.5,0.,  0.,0,-.5,  0.,.5,0.   )
D                = SMatrix{3,5}(1.,0.,0.,  1.,0.,0.,  1.,0.,0.,   1.,0.,0.,  1/√2,0,1/√2)
L                = 0.1
instrumentedbeam = StrainGaugeOnEulerBeam3D(elnod;P,D,elementkwargs=(mat=mat,orient2=SVector(0.,1.,0.)))
@testset "constructor" begin
    @test instrumentedbeam.eleobj.cₘ    ≈ [2.0, 0.0, 0.0]
    @test instrumentedbeam.eleobj.rₘ    ≈ I
    @test instrumentedbeam.eleobj.ζnod  ≈ [-0.5, 0.5]
    @test instrumentedbeam.eleobj.tgₘ   ≈ [4.0, 0.0, 0.0]
    @test instrumentedbeam.eleobj.tgₑ   ≈ [4.0, 0.0, 0.0]
    @test instrumentedbeam.P            ≈ P
    @test instrumentedbeam.D            ≈ D
    @test instrumentedbeam.E            ≈ SVector(1.,1.,1.,1.,.5) 
    @test instrumentedbeam.K1           ≈ SVector(0.,0.,0.,0.,.25) 
    @test instrumentedbeam.K2           ≈ SVector(-0.5, 0, 0.5, 0,-.25)
    @test instrumentedbeam.K3           ≈ SVector(0, -0.5, 0, 0.5,0) 
    @test typeof(instrumentedbeam)      ==  StrainGaugeOnEulerBeam3D{5, EulerBeam3D{BeamCrossSection, false}, @NamedTuple{strain::@NamedTuple{ε::Nothing, κ::Nothing}}}
end

Λ   =  SVector(0,0,0, 0,.1,0, 0,0,0, 0,-.1,0)
X   = (SVector(0,0,0, 0,.1,0, 0,0,0, 0,-.1,0),)
U   = (SVector{0,𝕣}(),)
A   = SVector{0,𝕣}()

X   = (SVector(0,0,0, 0,.1,0, 0,0,0, 0,-.1,0),)
R,FB,eleres = Muscade.residual(instrumentedbeam,X,U,A,0.,nothing,(;),@request (εₐₓ,κ,ε))
@testset "strain1" begin
    @test eleres.εₐₓ  ≈ 0.               atol = 1e-12
    @test eleres.κ*20 ≈ [0, 0, 1]        atol = 1e-12
    @test eleres.ε*40 ≈ [0, -1, 0, 1,0]  atol = 1e-12
end

X   = (SVector(0,0,0, 0,0,.1, 0,0,0, 0,0,-.1),)
R,FB,eleres = Muscade.residual(instrumentedbeam,X,U,A,0.,nothing,(;),@request (εₐₓ,κ,ε))
@testset "strain2" begin
    @test eleres.εₐₓ  ≈ 0.               atol = 1e-12
    @test eleres.κ*20 ≈ [0,-1, 0]        atol = 1e-12
    @test eleres.ε*40 ≈ [1, 0,-1, 0,.5]  atol = 1e-12
end

X   = (SVector(0,0,0, 0,0,0, 0,0,0, 1.,0,0),)
R,FB,eleres = Muscade.residual(instrumentedbeam,X,U,A,0.,nothing,(;),@request (εₐₓ,κ,ε))
@testset "strain3" begin
    @test eleres.εₐₓ  ≈ 0.               atol = 1e-12
    @test eleres.κ ≈ [.25,0, 0]          atol = 1e-12
    @test eleres.ε ≈ [0, 0,0, 0,0.0625]  atol = 1e-12
end

if true
    axis    = Muscade.SpyAxis()
else
    using GLMakie
    fig      = Figure(size = (500,500))
    display(fig) # open interactive window (gets closed down by "save")
    axis      = Axis3(fig[1,1],title="Test",xlabel="1",ylabel="2",zlabel="3",aspect=:data,viewmode=:free,perspectiveness=.5,clip=false)
end


X       = (SVector(0,0,0, 0,.1,0, 0,0,0, 0,-.1,0),)

Λm      =         reshape(Λ ,(length(Λ ),1))
Xm      = map(Xᵢ->reshape(Xᵢ,(length(Xᵢ),1)),X)
Um      = map(Uᵢ->reshape(Uᵢ,(length(Uᵢ),1)),U)
Am      = map(Aᵢ->reshape(Aᵢ,(length(Aᵢ),1)),A)

α       = 2π*(0:19)/20
circle  = 0.5*[cos.(α) sin.(α)]'
mut,opt = Muscade.allocate_drawing(axis,[instrumentedbeam];EulerBeam3D=(;style=:solid, nseg=1, section=circle, marking=true, Udof=false, Uscale=0.1),
                                                           StrainGaugeOnEulerBeam3D = (;L=L)) 
mut     = Muscade.update_drawing(  axis,[instrumentedbeam],mut,opt, Λm,Xm,Um,Am,0.,nothing,(;)) 
_       = Muscade.display_drawing!(axis,typeof(instrumentedbeam),mut,opt)                          

@testset "draw" begin
    @test axis.call[1].fun == :lines!
    @test axis.call[1].args[1][:,1:2] ≈ [2.1          1.9;
                                         0.51         0.51;
                                         1.7632e-18  -1.7632e-18]
    @test axis.call[2].fun == :scatter! # call to EulerBeam3D took place
end

@functor (;) function straincost(eleres,X,U,A,t) 
    σ  = 15e-6
    ε  = eleres.ε
    εₘ = SVector(cos(t),0.,-cos(t),0.,cos(t)/2)*0.001  
    Δε = ε-εₘ
    cost = (Δε⋅Δε)/(2σ^2)
    return cost
end

costedbeam =  ElementCost(elnod;
                            req = @request(ε),
                            cost=straincost,
                            ElementType=StrainGaugeOnEulerBeam3D,
                            elementkwargs = (P,D,
                                              elementkwargs=(mat=mat,orient2=SVector(0.,1.,0.))))

out = Muscade.diffed_lagrangian(costedbeam;Λ,X,U,A,t=0.)

@testset "costedbeam" begin
    @test costedbeam.eleobj == instrumentedbeam
    @test  out.∇L[2][1] ≈ [  277777.7777783019,       0.0,       1.1102230246251565e-16,  201441.02435832855,       2.7777777927777793e7,      -1.2430497627256343e6, -277777.7777783019,       0.0,      -1.1102230246251565e-16,  -76336.75341947998,      -2.7777777927777793e7,       1.2569502372743965e6]
    @test  out.HL[2,2][1,1] ≈ [  1.18056e9   0.0  0.0   3.64598e7    1.92747e-9   -3.29847e7    -1.18056e9   0.0  0.0  -3.29847e7     9.63735e-8   3.64598e7;
                                 0.0         0.0  0.0   0.0          0.0           0.0           0.0         0.0  0.0   0.0           0.0          0.0;
                                 0.0         0.0  0.0   0.0          0.0           0.0           0.0         0.0  0.0   0.0           0.0          0.0;
                                 3.64598e7   0.0  0.0   1.94899e7    1.01196e-9   -2.42679e7    -3.64598e7   0.0  0.0  -1.69698e7     2.02392e-9   2.60924e7;
                                 1.92747e-9  0.0  0.0   1.01196e-9   1.38889e8    -9.15508e-10  -1.92747e-9  0.0  0.0  -9.15508e-10  -1.38889e8    1.01196e-9;
                                -3.29847e7   0.0  0.0  -2.42679e7   -9.15508e-10   1.54556e8     3.29847e7   0.0  0.0   8.71679e6    -1.83102e-9  -1.56207e8;
                                -1.18056e9   0.0  0.0  -3.64598e7   -1.92747e-9    3.29847e7     1.18056e9   0.0  0.0   3.29847e7    -9.63735e-8  -3.64598e7;
                                 0.0         0.0  0.0   0.0          0.0           0.0           0.0         0.0  0.0   0.0           0.0          0.0;
                                 0.0         0.0  0.0   0.0          0.0           0.0           0.0         0.0  0.0   0.0           0.0          0.0;
                                -3.29847e7   0.0  0.0  -1.69698e7   -9.15508e-10   8.71679e6     3.29847e7   0.0  0.0   1.60148e7    -1.83102e-9  -1.03674e7;
                                 9.63735e-8  0.0  0.0   2.02392e-9  -1.38889e8    -1.83102e-9   -9.63735e-8  0.0  0.0  -1.83102e-9    1.38889e8    2.02392e-9;
                                 3.64598e7   0.0  0.0   2.60924e7    1.01196e-9   -1.56207e8    -3.64598e7   0.0  0.0  -1.03674e7     2.02392e-9   1.58031e8] rtol = 1e-4
end

end