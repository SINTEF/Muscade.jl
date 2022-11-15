using Test,Documenter,Muscade
@testset "Muscade.jl package" begin
    @testset "TestElementAPI" begin
        include("TestElementAPI.jl")
    end
    @testset "TestAdiff" begin
        include("TestAdiff.jl")
    end
    @testset "TestModelDescription" begin
        include("TestModelDescription.jl")
    end
    @testset "TestAssemble" begin
        include("TestAssemble.jl")
    end
    @testset "TestStaticX" begin
        include("TestStaticX.jl")
    end
    @testset "TestUnit" begin
        include("TestUnit.jl")
    end
    #doctest(Muscade)
end
;