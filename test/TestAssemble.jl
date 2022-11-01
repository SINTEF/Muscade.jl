module TestElementAPI

using Test,StaticArrays
using Muscade
using Muscade.ElTest

include("SomeElements.jl")

### Turbine
sea(t,x) = SVector(1.,0.)
sky(t,x) = SVector(0.,1.)
turbine  = Turbine(SVector(0.,0.),-10., 2.,sea, 3.,sky)
δX       = @SVector [1.,1.]
X        = @SVector [1.,2.]
U        = @SVector []
A        = @SVector [0.,0.]  # [Δseadrag,Δskydrag]


L,Lδx,Lx,Lu,La   = gradient(SeverΛXUAstatic,turbine,δX,[X],[U],A, 0.,0.,())

@testset "Turbine gradient" begin
    @test Lδx           ≈ [2, 3]
    @test Lx            ≈ [0, 0]
    @test length(Lu)    == 0
    @test La            ≈ [1, 1]
end

Lδx,Lx,Lu,La   = test_static_element(turbine;δX,X,U,A,verbose=false)

@testset "test_static_element" begin
    @test Lδx           ≈ [2, 3]
    @test Lx            ≈ [0, 0]
    @test length(Lu)    == 0
    @test La            ≈ [1, 1]
end

L,Ly,Lyy   = hessian(JointΛXAstatic,turbine,δX,[X],[U],A, 0.,0.,())

@testset "hessian" begin
    @test L           ≈ 5.
    @test Ly          ≈ [2,3,0,0,1,1]
    @test Lyy         ≈ [0 0 0 0 1 0;0 0 0 0 0 1;0 0 0 0 0 0;0 0 0 0 0 0;1 0 0 0 0 0;0 1 0 0 0 0]
end

# ###  AnchorLine

anchorline      = AnchorLine(SVector(0.,0.,100.), SVector(0,2.,0), SVector(94.,0.), 170., -1.)

δX       = @SVector [1.,1.,1.]
X        = @SVector [0.,0.,0.]
U        = @SVector []
A        = @SVector [0.,0.]  # [Δseadrag,Δskydrag]
L,Lδx,Lx,Lu,La   = gradient(SeverΛXUAstatic,anchorline,δX,[X],[U],A, 0.,0.,())
@testset "anchorline1" begin
    @test Lδx           ≈ [-12.256289016934003, 0.26077210674327667, -0.5215442134865533]
    @test Lx            ≈ [0.9150974560878556, -0.14708204066347275, 22.682383121692297]
    @test length(Lu)    == 0
    @test La            ≈ [0.9180190940688681, 12.51706112367728]
end

end
