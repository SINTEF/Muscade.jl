module TestPositionElement
using Test, Muscade, Muscade.Toolbox, StaticArrays, LinearAlgebra

model            = Model(:TestModel)
xₘ               = SVector{3,𝕣}(0,0,1)
node             = addnode!(model,xₘ)
elnod            = [model.nod[node]]
P                = SMatrix{3,3,𝕣}(1,0,0,  0,1,0,   0,0,1)
D                = SMatrix{3,3,𝕣}(0,1,0,  0,0,1,   1,0,0)
position         = Position3D(elnod;P,D)
@testset "constructor" begin
    @test position.xₘ           ≈ xₘ 
    @test position.P            ≈ P
    @test position.D            ≈ D
end

Λ   =  SVector(0,0,0, 0,.1,0 )
X   = (SVector(0,1,0, 0,.1,0),SVector(0,0,0, .5,0,0),SVector(0,0,0, .25,0,0))
U   = (SVector{0,𝕣}(),)
A   = SVector{0,𝕣}()

R,FB,eleres = Muscade.residual(position,X,U,A,0.,nothing,(;),@request (R,x,rᵢ,rₑ,a))
@testset "requestables" begin
    @test eleres.R[1]  ≈ [0.995004   0.0  0.0998334;
                          0.0        1.0  0.0      ;
                         -0.0998334  0.0  0.995004]                        rtol = 1e-4
    @test eleres.R[2]  ≈ [0.0        0.0249792   0.0;
                          0.0249792  0.0        -0.499167;
                          0.0        0.499167    0.0]                      rtol = 1e-4
    @test eleres.R[3]  ≈ [0.000208194   0.0124896  -0.008325;
                          0.0124896    -0.249792   -0.249584;
                          0.008325      0.249584   -0.249584]              rtol = 1e-4
    @test eleres.x[1]  ≈ [0.995004  0.0  0.0998334;
                          1.0       2.0  1.0;
                          0.900167  1.0  1.995]                            rtol = 1e-4
    @test eleres.x[2]  ≈ [0.0        0.0249792   0.0;
                          0.0249792  0.0        -0.499167;
                          0.0        0.499167    0.0]                      rtol = 1e-4
    @test eleres.x[3]  ≈ [0.000208194   0.0124896  -0.008325;
                          0.0124896    -0.249792   -0.249584;
                          0.008325      0.249584   -0.249584]              rtol = 1e-4
    @test eleres.rᵢ[1] ≈ [0.0, 0.0, 0.0]                                   rtol = 1e-12
    @test eleres.rᵢ[2] ≈ [0.49916708323414083, 0.0, 0.024979173609871178]  rtol = 1e-12
    @test eleres.rᵢ[3] ≈ [0.24958354161707041, 0.004164583829296608, 0.012489586804935589]  rtol = 1e-12
    @test eleres.rₑ[1] ≈ [0.0, 0.1, 0.0]                                   rtol = 1e-12
    @test eleres.rₑ[2] ≈ [0.5, 0.0, 0.0]                                   rtol = 1e-12
    @test eleres.rₑ[3] ≈ [0.25, 0.0, 0.0]                                  rtol = 1e-12
    @test eleres.a  ≈ [0.012489586804935589 -0.24979173609871178 -0.24958354161707041; 
                        0.008325002975638984 0.24958354161707041 -0.2495835416170705; 
                        0.00020819448164127675 0.012489586804935589 -0.008325002975638984]                    rtol = 1e-12
end

if true
    axis    = Muscade.SpyAxis()
else
    using GLMakie
    fig      = Figure(size = (500,500))
    display(fig) # open interactive window (gets closed down by "save")
    axis      = Axis3(fig[1,1],title="Test",xlabel="1",ylabel="2",zlabel="3",aspect=:data,viewmode=:free,perspectiveness=.5,clip=false)
end

X       = ([0,0,0, 0,.1,0],)

Λm      =         reshape(Λ ,(length(Λ ),1))
Xm      = map(Xᵢ->reshape(Xᵢ,(length(Xᵢ),1)),X)
Um      = map(Uᵢ->reshape(Uᵢ,(length(Uᵢ),1)),U)
Am      = map(Aᵢ->reshape(Aᵢ,(length(Aᵢ),1)),A)

mut,opt = Muscade.allocate_drawing(axis,[position];Position3D=(;L=.1)) 
mut     = Muscade.update_drawing(  axis,[position],mut,opt, Λm,Xm,Um,Am,0.,nothing,(;)) 
_       = Muscade.display_drawing!(axis,Position3D,mut,opt)                          

@testset "draw" begin
    @test axis.call[1].fun == :scatter!
    @test axis.call[1].args[1][:,1:2] ≈ [0.995004  0.0;
                                         0.0       1.0;
                                         0.900167  1.0] rtol = 1e-4
    @test axis.call[2].fun == :lines!
    @test axis.call[2].args[1][:,1:2] ≈ [0.0 0.995004;
                                         0.0       0.0;
                                         1.0  0.900167] rtol = 1e-4
    @test axis.call[3].fun == :lines!
    @test axis.call[3].args[1][:,1:2] ≈ [0.995004  0.995004;
                                         0.0       0.1;
                                         0.900167  0.900167] rtol = 1e-4
end

end