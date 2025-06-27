module TestElementAPI

using Test,StaticArrays
using Muscade

include("SomeElements.jl")

### Turbine
sea(t,x) = SVector(1.,0.)
sky(t,x) = SVector(0.,1.)
turbine  = Turbine(SVector(0.,0.),-10., 2.,sea, 3.,sky)
X        = @SVector [1.,2.]
U        = @SVector []
A        = @SVector [0.,0.]  # [Δseadrag,Δskydrag]

R,FB = Muscade.residual(turbine, [X],[U],A, 0.,nothing,(dbg=true,))
@testset "Turbine" begin
    @test R             ≈ [-2, -3]
end
T = typeof(turbine)
@testset "Element utility functions" begin
    @test Muscade.getdoflist(T)  == ((1, 1, 2, 2),(:X, :X, :A, :A),(:tx1, :tx2, :Δseadrag, :Δskydrag))
    @test Muscade.getidof(T,:X)  == [1,2]
    @test Muscade.getidof(T,:U)  == []
    @test Muscade.getidof(T,:A)  == [3,4]
    @test Muscade.getndof(T)     == 4
    @test Muscade.getndof(T,:X)  == 2
    @test Muscade.getndof(T,(:X,:U,:A))    == (2,0,2)
    @test Muscade.getnnod(T)     == 2
end



# ###  AnchorLine

anchorline      = AnchorLine(SVector(0.,0.,100.), SVector(0,2.,0), SVector(94.,0.), 170., -1.)

δX       = @SVector [1.,1.,1.]
X        = @SVector [0.,0.,0.]
U        = @SVector []
A        = @SVector [0.,0.]  # [Δseadrag,Δskydrag]
L1,FB    = Muscade.lagrangian(anchorline, δX,[X],[U],A, 0.,nothing,(dbg=true,))
X        = [0,-1,45/180*π]
L2,FB    = Muscade.lagrangian(anchorline, δX,[X],[U],A, 0.,nothing,(dbg=true,))

@testset "Lagrangian" begin
   @test L1 ≈ 12.517061123678818
   @test L2 ≈ 5.590087401683872
end


### Spring

T = Spring{2}
@testset "Element utility functions" begin
    @test Muscade.getdoflist(T)  == ((1, 1, 2, 2, 3, 3), (:X, :X, :X, :X, :A, :A), (:tx1, :tx2, :tx1, :tx2, :ΞL₀, :ΞEI))
    @test Muscade.getidof(T,:X)  == [1,2,3,4]
    @test Muscade.getidof(T,:U)  == []
    @test Muscade.getidof(T,:A)  == [5,6]
    @test Muscade.getndof(T)     == 6
    @test Muscade.getndof(T,:X)  == 4
    @test Muscade.getndof(T,(:X,:U,:A))    == (4,0,2)
    @test Muscade.getnnod(T)     == 3
end


end
