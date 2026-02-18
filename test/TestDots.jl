module TestDots
using Muscade
using Test,StaticArrays


@testset "SArrays" begin
    @testset "doubledot" begin
        a = SArray{Tuple{2,3,4}}(1:24)
        b = SArray{Tuple{3,4,2,2}}(1:48)
        c = dots(a,b,Val(2))
        @test size(c)==(2,2,2)
    end
    @testset "matmat" begin
        a = SArray{Tuple{2,3}}(1:6)
        b = SArray{Tuple{3,4}}(1:12)
        c = dots(a,b,Val(1))
        @test c==a*b
    end
    @testset "external" begin
        a = SArray{Tuple{2,3}}(1:6)
        b = SArray{Tuple{3,4}}(1:12)
        c = dots(a,b,Val(0))
        d = dots(a',b,Val(0))
        e = dots(conj(a),b',Val(0))
        @test size(c)==(2,3,3,4)
        @test size(d)==(3,2,3,4)
        @test size(e)==(2,3,4,3)
    end
    @testset "doubledot" begin
        a = SArray{Tuple{2,3}}(1:6)
        c = dots(a,a,Val(2))
        @test c==91
    end
end

@testset "Arrays" begin
    @testset "doubledot" begin
        a = reshape(1:24,2,3,4)
        b = reshape(1:48,3,4,2,2)
        c = dots(a,b,Val(2))
        @test size(c)==(2,2,2)
    end
    @testset "matmat" begin
        a = reshape(1:6,2,3)
        b = reshape(1:12,3,4)
        c = dots(a,b,Val(1))
        @test c==a*b
    end
    @testset "external" begin
        a = reshape(1:6,2,3)
        b = reshape(1:12,3,4)
        c = dots(a,b,Val(0))
        d = dots(a',b,Val(0))
        e = dots(conj(a),b',Val(0))
        @test size(c)==(2,3,3,4)
        @test size(d)==(3,2,3,4)
        @test size(e)==(2,3,4,3)
    end
    @testset "doubledot" begin
        a = reshape(1:6,2,3)
        c = dots(a,a,Val(2))
        @test c==91
    end
end

end