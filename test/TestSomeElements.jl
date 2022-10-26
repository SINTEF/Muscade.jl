module TestSomeElements

using Test,StaticArrays
using Muscade.Tools.ElementTestBench

include("SomeElements.jl")

@testset "SomeElements" begin
### Turbine
sea(t,x) = SVector(1.,0.)
sky(t,x) = SVector(0.,1.)
turbine  = Turbine(SVector(0.,0.),-10., 2.,sea, 3.,sky)
δX       = @SVector [1.,1.]
X        = @SVector [1.,2.]
U        = @SVector []
A        = @SVector [0.,0.]  # [Δseadrag,Δskydrag]
Lδx,Lx,Lu,La  = testStaticElement(turbine; δX,X,U,A,verbose=false)
@testset "Turbine" begin
    @test Lδx           ≈ [2, 3]
    @test Lx            ≈ [0, 0]
    @test length(Lu)    == 0
    @test La            ≈ [1, 1]
end

# ###  AnchorLine

anchorline      = AnchorLine(SVector(0.,0.,100.), SVector(0,2.,0), SVector(94.,0.), 170., -1.)

δX       = @SVector [1.,1.,1.]
X        = @SVector [0.,0.,0.]
U        = @SVector []
A        = @SVector [0.,0.]  # [Δseadrag,Δskydrag]
Lδx,Lx,Lu,La  = testStaticElement(anchorline; δX,X,U,A,verbose=false)
@testset "anchorline1" begin
    @test Lδx           ≈ [-12.256289016934003, 0.26077210674327667, -0.5215442134865533]
    @test Lx            ≈ [0.9150974560878556, -0.14708204066347275, 22.682383121692297]
    @test length(Lu)    == 0
    @test La            ≈ [0.9180190940688681, 12.51706112367728]
end

X               = [0,-1,45/180*π]
Lδx,Lx,Lu,La  = testStaticElement(anchorline; δX,X,U,A,verbose=false)

@testset "anchorline2" begin
    @test Lδx           ≈ [-13.554507845277369, 0.05884302528088421, 19.08575222168048]
    @test Lx            ≈ [-0.3964022106464097, -0.05737760846863315, 19.893927005010767]
    @test length(Lu)    == 0
    @test La            ≈ [-0.39614938921718024, -5.590087401683997]
end



end # @testset





end
;