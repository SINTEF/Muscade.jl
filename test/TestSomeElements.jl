module TestOptimist

include("Optimist.jl")
using Muscade,Test,StaticArrays
using Muscade.Tools.Unit


node   = nodesforelementtest([0.,0.,0.])

### Turbine

turbine  = Turbine(node, sea=sea, sky=sky, seadrag=2., skydrag=3.)
#req      = @request (:X,)
δX       = @SVector [1.,1.]
X        = @SVector [1.,2.]
U        = @SVector []
A        = @SVector [0.,0.]  # [Δseadrag,Δskydrag]
Lδx,Lx,Lu,La,χo,χn  = testStaticElement(turbine; δX,X,U,A)
@testset "Turbine" begin
    @test r           ≈ [2, 3, 0]
    @test out[key.X]  ≈ [0, 0, 0]
end

###  AnchorLine

el              = AnchorLine(node, ΔX₀top=[0,2.,0]←m, X₀bot=[100.,0]←m, L=170←m, buoyancy=-1←N/m)
req             = @request (:Xtop,:ΔXchain,:xaf,:cr,:Fh,:ltf)
y               = [[0←m,0←m,0←deg]]
r,χo,χn,out,key = testelement(el,y=y,req=req,verbose=false)
@testset "AnchorChain1" begin
    @test r            ≈ [-18.49678868176869, 0.36993577363537383, 36.99357736353738]
    @test out[key.cr]  ≈ 18.50048766964324
    @test out[key.ltf] ≈ 117.04741575074884
end
y               = [[0←m,-1←m,45←deg]]
r,χo,χn,out,key = testelement(el,y=y,req=req,verbose=false)
@testset "AnchorChain2" begin
    @test r            ≈ [-20.227795576207416, 0.08261787939047101, 28.489583515254516]
    @test out[key.cr]  ≈ 20.22796429665702
    @test out[key.ltf] ≈ 118.51410405235067
end





#end # @testset





end
;