module TestElementAPI

using Test,StaticArrays
using Muscade
using Muscade.ElTest

include("SomeElements.jl")

### Turbine
sea(t,x) = SVector(1.,0.)
sky(t,x) = SVector(0.,1.)
turbine  = Turbine(SVector(0.,0.),-10., 2.,sea, 3.,sky)
X        = @SVector [1.,2.]
U        = @SVector []
A        = @SVector [0.,0.]  # [Δseadrag,Δskydrag]

Re = zeros(2)
residual(turbine, Re,[X],[U],A, 0.,0.,())
@testset "Turbine" begin
    @test Re            ≈ [2, 3]
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



# ###  AnchorLine

anchorline      = AnchorLine(SVector(0.,0.,100.), SVector(0,2.,0), SVector(94.,0.), 170., -1.)

δX       = @SVector [1.,1.,1.]
X        = @SVector [0.,0.,0.]
U        = @SVector []
A        = @SVector [0.,0.]  # [Δseadrag,Δskydrag]
L1       = lagrangian(anchorline, δX,[X],[U],A, 0.,0.,())
X               = [0,-1,45/180*π]
L2       = lagrangian(anchorline, δX,[X],[U],A, 0.,0.,())

@testset "Lagrangian" begin
   @test L1 ≈ -12.517061123678818
   @test L2 ≈ 5.590087401683872
end

end
