module TestOptimist
include("Optimist.jl")
using Muscade,Test,StaticArrays

using Muscade.Tools.Unit

#nodes    = nodesforelementtest([0.,0.,0.])

nodes = [Node(zeros(3))]
sea(t,x) = SVector(10,0)
el       = Ballast(nodes,seadrag=3e3, sea=sea, buoyancy=-1e4)
δX       = @SVector [1.,1,1]
X        = @SVector [1.,2,3]
U        = @SVector []
A        = @SVector [0.,0]
Lδx,Lx,Lu,La,χo,χn  = testStaticElement(el; δX,X,U,A)

@testset "Ballast" begin
    @test Lδx ≈ [3e4,0,1e4]
    @test Lx  ≈ [0.,0.,0.]
    @test La  ≈ [10.,-1.]
end

# ###############

# sea(t,x) = SVector(1.,0.)
# sky(t,x) = SVector(0.,1.)

# nod             = nodesforelementtest([0. 0. 100.;50. 0. 80.; 60. 0. 80.]←m)
# nod1            = nod[1:1]

# @testset "Optimist unit tests" begin
# ### Ballast
# el              = Ballast(nod1, sea=sea, seadrag=2. ←N/(m/s),buoyancy=-100←N)
# req             = @request (:X,)
# y               = [[1,2,3]←m]
# r,χo,χn,out,key = testelement(el,y=y,req=req,verbose=false)
# @testset "Ballast" begin
#     @test r           ≈ [2, 0, 100]
#     @test out[key.X]  ≈ [1, 2, 103]
# end


# ### Turbine
# el              = Turbine(nod1, sea=sea, sky=sky, seadrag=2., skydrag=3.)
# req             = @request (:X,)
# y               = [[0←m,0←m,0←deg]]
# r,χo,χn,out,key = testelement(el,y=y,req=req,verbose=false)
# @testset "Turbine" begin
#     @test r           ≈ [2, 3, 0]
#     @test out[key.X]  ≈ [0, 0, 0]
# end


# ### Anchor chain
# el              = AnchorLine(nod1, ΔX₀top=[0,2.,0]←m, X₀bot=[100.,0]←m, L=170←m, buoyancy=-1←N/m)
# req             = @request (:Xtop,:ΔXchain,:xaf,:cr,:Fh,:ltf)
# y               = [[0←m,0←m,0←deg]]
# r,χo,χn,out,key = testelement(el,y=y,req=req,verbose=false)
# @testset "AnchorChain1" begin
#     @test r            ≈ [-18.49678868176869, 0.36993577363537383, 36.99357736353738]
#     @test out[key.cr]  ≈ 18.50048766964324
#     @test out[key.ltf] ≈ 117.04741575074884
# end
# y               = [[0←m,-1←m,45←deg]]
# r,χo,χn,out,key = testelement(el,y=y,req=req,verbose=false)
# @testset "AnchorChain2" begin
#     @test r            ≈ [-20.227795576207416, 0.08261787939047101, 28.489583515254516]
#     @test out[key.cr]  ≈ 20.22796429665702
#     @test out[key.ltf] ≈ 118.51410405235067
# end


# ### HawserTop
# el              = HawserTop(nod[1:2], ΔX₀top = [0.,2.,0.],L₀=52←m ,EA=10000←N)
# req             = @request (:Xtop,:Xbot,:ΔXtop,:ΔX,:L,:T)
# y               = [[0←m,0←m,0←deg, 0←m,0←m,0←m]]
# r,χo,χn,out,key = testelement(el,y=y,req=req,verbose=false)
# @testset "HawserTop1" begin
#     @test r              ≈ [-337.01437756954743, 13.480575102781897, 674.0287551390949, 337.01437756954743, -13.480575102781897, -134.80575102781896]
#     @test out[key.Xtop]  ≈ [0, 0, 100]
#     @test out[key.Xbot]  ≈ [50, 0, 80]
#     @test out[key.ΔXtop] ≈ [0,2, 0]
#     @test out[key.L]     ≈ 53.88877434122992
#     @test out[key.T]     ≈ 363.22583485190796
# end
# y               = [[0←m,-1←m,45←deg,0←m,0←m,0←m]]
# r,χo,χn,out,key = testelement(el,y=y,req=req,verbose=false)
# @testset "HawserTop2" begin
#     @test r              ≈ [-567.9066737579479, 4.575284345200771, 796.6708910179865, 567.9066737579479, -4.575284345200771, -220.9142703579401]
#     @test out[key.Xtop]  ≈ [0.0, -1.0, 100.0]
#     @test out[key.Xbot]  ≈ [50, 0, 80]
#     @test out[key.ΔXtop] ≈ [-1.414213562373095, 1.4142135623730951, 0.0]
#     @test out[key.L]     ≈ 55.168767696157246
#     @test out[key.T]     ≈ 609.3784031071636
# end


# ### Hawser
# el              = Hawser(nod[1:2], L₀=52←m ,EA=10000←N)
# req             = @request (:X1,:X2,:L,:T)
# y               = [[0←m,0←m,0←deg, 0←m,0←m,0←m]]
# r,χo,χn,out,key = testelement(el,y=y,req=req,verbose=false)
# @testset "Hawser1" begin
#     @test r              ≈ [-330.61770653202234, -0.0, 132.24708261280895, 330.61770653202234, 0.0, -132.24708261280895]
#     @test out[key.X1]    ≈ [0, 0, 100]
#     @test out[key.X2]    ≈ [50, 0, 80]
#     @test out[key.L]     ≈ 53.85164807134504
#     @test out[key.T]     ≈ 356.08616756635405
# end
# y               = [[0←m,-1←m,0←m,0←m,0←m,0←m]]
# r,χo,χn,out,key = testelement(el,y=y,req=req,verbose=false)
# @testset "Hawser2" begin
#     @test r              ≈ [-332.2181145261312, -6.644362290522625, 132.8872458104525, 332.2181145261312, 6.644362290522625, -132.8872458104525]
#     @test out[key.X1]    ≈ [0.0, -1.0, 100.0]
#     @test out[key.X2]    ≈ [50, 0, 80]
#     @test out[key.L]     ≈ 53.86093203798092
#     @test out[key.T]     ≈ 357.87154576556236
# end

#end # @testset





end
;