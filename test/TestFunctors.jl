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

@testset "functors" begin
    @test typeof(f)    == Functor{:f, @NamedTuple{a::Int64, e::Int64}}
    @test typeof(g)    == Functor{:g, @NamedTuple{a::Int64, e::Int64}}
    @test typeof(h)    == Functor{:h, @NamedTuple{}}
    @test f isa Functor
    @test f isa Function
    @test f(2.) ≈ 12.
    @test g(2.) ≈ 12.
    @test h(2.) ≈ 4.
    @test f(2) == 12
    @test fu(1.,0,fukwargs...) ≈ 1.
    @inferred f(2.)
    @inferred f(2)
end

end

