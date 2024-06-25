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

R,χ,FB =residual(turbine, [X],[U],A, 0.,nothing,identity,nothing,(dbg=true,))
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
L1,χ,FB  = lagrangian(anchorline, δX,[X],[U],A, 0.,nothing,nothing,(dbg=true,))
X        = [0,-1,45/180*π]
L2,χ,FB  = lagrangian(anchorline, δX,[X],[U],A, 0.,nothing,nothing,(dbg=true,))

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


const X1 = (SVector{3,𝕣}(1,2,3),)
const X2 = (SVector{3,𝕣}(1,2,3),SVector{3,𝕣}(4,5,6))
const X3 = (SVector{3,𝕣}(1,2,3),SVector{3,𝕣}(4,5,6),SVector{3,𝕣}(7,8,9))
const X4 = (variate{1,3}(SVector{3,𝕣}(1,2,3)),)
const X5 = (variate{1,3}(SVector{3,𝕣}(1,2,3)),variate{1,3}(SVector{3,𝕣}(4,5,6)))
const X6 = (variate{1,3}(SVector{3,𝕣}(1,2,3)),variate{1,3}(SVector{3,𝕣}(4,5,6)),variate{1,3}(SVector{3,𝕣}(7,8,9)))
const P  = 2
Y1=Muscade.motion{P,1}(X1)
Y2=Muscade.motion{P,2}(X2)
Y3=Muscade.motion{P,3}(X3)
Y4=Muscade.motion{P,1}(X4)
Y5=Muscade.motion{P,2}(X5)
Y6=Muscade.motion{P,3}(X6)
@testset "motion" begin
    @test Muscade.position{P,1}(Y1) === SVector{3,𝕣}(1,2,3)
    @test Muscade.velocity{P,1}(Y1) === SVector{3,𝕣}(0,0,0)
    @test Muscade.acceleration{P,1}(Y1) === SVector{3,𝕣}(0,0,0)

    @test Muscade.position{P,2}(Y2) === SVector{3,𝕣}(1,2,3)
    @test Muscade.velocity{P,2}(Y2) === SVector{3,𝕣}(4,5,6)
    @test Muscade.acceleration{P,2}(Y2) === SVector{3,𝕣}(0,0,0)

    @test Muscade.position{P,3}(Y3) === SVector{3,𝕣}(1,2,3)
    @test Muscade.velocity{P,3}(Y3) === SVector{3,𝕣}(4,5,6)
    @test Muscade.acceleration{P,3}(Y3) === SVector{3,𝕣}(7,8,9)

    @test Muscade.position{P,1}(Y4) === variate{1,3}(SVector{3,𝕣}(1,2,3))
    @test Muscade.velocity{P,1}(Y4) === SVector{3,𝕣}(0,0,0)
    @test Muscade.acceleration{P,1}(Y4) === SVector{3,𝕣}(0,0,0)

    @test Muscade.position{P,2}(Y5) === variate{1,3}(SVector{3,𝕣}(1,2,3))
    @test Muscade.velocity{P,2}(Y5) === variate{1,3}(SVector{3,𝕣}(4,5,6))
    @test Muscade.acceleration{P,2}(Y5) === SVector{3,𝕣}(0,0,0)

    @test Muscade.position{P,3}(Y6) === variate{1,3}(SVector{3,𝕣}(1,2,3))
    @test Muscade.velocity{P,3}(Y6) === variate{1,3}(SVector{3,𝕣}(4,5,6))
    @test Muscade.acceleration{P,3}(Y6) === variate{1,3}(SVector{3,𝕣}(7,8,9))
end

end
