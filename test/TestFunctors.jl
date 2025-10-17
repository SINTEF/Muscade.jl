module TestFunctors
using Muscade
using Test

a = 3
@functor (a,e=2)  g(x::Real)=a*x^e
@functor (a,e=2) function f(x::Real)
    return a*x^e
end
@functor (;) h(x) = 2x
@functor (;) fu(u,t) = u^2
fukwargs = (;) 
a = :a

a = 2 # changing a: test capture by value, not reference
b = 1
@functor ()        cost1(x) = x
@functor (a)       cost2(x) = a*x
@functor (a=2)     cost3(x) = a*x
@functor (a,b=1)   cost4(x) = a*x+b
@functor (a=2,b=1) cost5(x) = a*x+b



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

end

end

