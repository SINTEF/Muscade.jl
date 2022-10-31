module TestUnit
include("../src/core/Unit.jl")

using Test

imp   = pound/foot^3
dou   = 2m
@testset "Unit" begin
    @test ( 4. ← pound/foot^3 )≈ 64.07385349584058
    @test ( 4. ← imp )≈ 64.07385349584058
    @test ( 2^2←imp )≈ 64.07385349584058
    @test ( 2+2←imp )≈ 64.07385349584058
    @test ( [1,2,3]←dou )≈ [2.,4.,6.]
    @test ( (1,2,3)←dou )== (2.,4.,6.)
    @test ( ([1,2],3)←dou )== ([2.,4.],6.)
    @test ( (1:3)←dou)   ≈ [2.,4.,6.]
    @test ( [[1,2],3]←dou)≈ [[2.,4.],6.]
    @test ( (2←imp)→pound/foot^3 )≈ 2.
    @test (2+2→imp )≈ (4→imp)
    @test string(imp )== "16.018463373960145 m^-3.0*kg"
end
end
