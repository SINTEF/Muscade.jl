module TestFunctors
using Muscade
using Test

a = 3
@functor with(a,e=2)  g(x::Real)=a*x^e
@functor with(a,e=2) function f(x::Real)
    return a*x^e
end
@functor with() h(x) = 2x
@functor with() fu(u,t) = u^2
fukwargs = (;) 
a = :a

a = 2 # changing a: test capture by value, not reference
b = 1
@functor with()        cost1(x) = x
@functor with(a)       cost2(x) = a*x
@functor with(a=2)     cost3(x) = a*x
@functor with(a,b=1)   cost4(x) = a*x+b
@functor with(a=2,b=1) cost5(x) = a*x+b
@functor with(a,b)     cost6(x) = a*x+b



@testset "functors" begin
    @test typeof(f)    == Functor{:f, @NamedTuple{a::Int64, e::Int64}}
    @test typeof(g)    == Functor{:g, @NamedTuple{a::Int64, e::Int64}}
    @test typeof(h)    == Functor{:h, @NamedTuple{}}
    @test f isa Functor
    @test f isa Function
    @test @inferred f(2.) ≈ 12.
    @test @inferred g(2.) ≈ 12.
    @test @inferred h(2.) ≈ 4.
    @test @inferred f(2) == 12
    @test @inferred fu(1.,0,fukwargs...) == 1.
    @test @inferred cost1(3) == 3
    @test @inferred cost2(3) == 6
    @test @inferred cost3(3) == 6
    @test @inferred cost4(3) == 7
    @test @inferred cost5(3) == 7
    @test @inferred cost6(3) == 7

end


n  = 5
X  = (0:n-1)#
Xd = [0.1,.95,2.,3.03,3.99]
Y1 = X
Y2 = [X';X'.+1]

A1 = interpolator(X ,Y1)                  # range
B1 = interpolator(Xd,Y1;quasirange=true ) # quasirange
C1 = interpolator(Xd,Y1;quasirange=false) # vector

A2 = interpolator(X ,Y2)
B2 = interpolator(Xd,Y2;quasirange=true )
C2 = interpolator(Xd,Y2;quasirange=false)

using Test
@testset "interpolators" begin
    @test A1(3.5) ≈ 3.5
    @test B1(3.5) ≈ 3.4895833333333335
    @test C1(3.5) ≈ 3.4895833333333335

    @test A1(-1.) ≈ -1.0
    @test B1(-1.) ≈ -1.2941176470588236
    @test C1(-1.) ≈ -1.2941176470588236

    @test A1(4.3) ≈ 4.3
    @test B1(4.3) ≈ 4.322916666666666
    @test C1(4.3) ≈ 4.322916666666666

    @test A2(4.3) ≈ 4.3 .+ [0,1]
    @test B2(4.3) ≈ 4.322916666666666  .+ [0,1]
    @test C2(4.3) ≈ 4.322916666666666  .+ [0,1]
end


end


 
