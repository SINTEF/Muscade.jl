module TestMultiplex

using Test
using Muscade
using Muscade: Multiplex,demultiplex
using StaticArrays

t  = Multiplex(1.:3.)
x  = 1. +sin(t)
y  = 3. *t
r  = SVector(x,y)
dr = demultiplex(r)
dx = demultiplex(x)

@testset "Multiplex Vector" begin
    @test dr ≈ [1.8414709848078965 1.9092974268256817 1.1411200080598671; 3. 6. 9.]
    @test dx ≈ [1.8414709848078965,1.9092974268256817,1.1411200080598671]
  #  @test 
end

st  = Multiplex(SVector(1.,2.,3.))
sx  = 1. +sin(st)
sy  = 3. *st
sr  = SVector(sx,sy)
se  = SMatrix{2,2}(sx,sy,-sy,sx)
sdr = demultiplex(sr)
sdx = demultiplex(sx)
@testset "Multiplex SVector" begin
    @test sdr ≈ [1.8414709848078965 1.9092974268256817 1.1411200080598671; 3. 6. 9.]
    @test sdx ≈ [1.8414709848078965,1.9092974268256817,1.1411200080598671]
    @test zero(sx) ≈ zero(𝕣)
end


end