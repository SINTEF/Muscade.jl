#module TestFunctors
using Muscade
using Test

f  = FunctionFromVector(0:.1:.9,[0.,1.,1.,1.,1.,2.,2.,2.,3.,4.])

q1 = QuadraticFunction(1.,1.)
q2 = QuadraticFunction(f,1.)

@testset "functors" begin
    @test f(.45)        ≈ 1.5
    @test q1(1.5)       ≈ 0.125
    @test q2(1.5,.45)   ≈ 0.
    @test f  isa Function
    @test q1 isa Function
    @test q2 isa Function
end

# end

