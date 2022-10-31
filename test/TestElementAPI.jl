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

@testset "Turbine" begin
    @test Lδx           ≈ [2, 3]
    @test Lx            ≈ [0, 0]
    @test length(Lu)    == 0
    @test La            ≈ [1, 1]
end
T = typeof(turbine)
@testset "Element utility functions" begin
    @test Muscade.getdoflist(T)  == ([1, 1, 2, 2],[:X, :X, :A, :A],[:tx1, :tx2, :Δseadrag, :Δskydrag])
    @test Muscade.getidof(T,:X) == [1,2]
    @test Muscade.getidof(T,:U) == []
    @test Muscade.getidof(T,:A) == [3,4]
    @test Muscade.getndof(T)     == 4
    @test Muscade.getndof(T,:X)  == 2
    @test Muscade.getndofs(T)    == (2,0,2)
    @test Muscade.getnnod(T)     == 2
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

X               = [0,-1,45/180*π]
L,Lδx,Lx,Lu,La   = gradient(SeverΛXUAstatic,anchorline,δX,[X],[U],A, 0.,0.,())

@testset "anchorline2" begin
    @test Lδx           ≈ [-13.554507845277369, 0.05884302528088421, 19.08575222168048]
    @test Lx            ≈ [-0.3964022106464097, -0.05737760846863315, 19.893927005010767]
    @test length(Lu)    == 0
    @test La            ≈ [-0.39614938921718024, -5.590087401683997]
end


end
