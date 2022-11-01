using Test,Documenter,Muscade
@testset "Muscade.jl package" begin
    @testset "SomeElements" begin
        include("TestElementAPI.jl")
        include("TestModelDescription.jl")
        include("TestAssemble.jl")
    end
    @testset "TestUnit" begin
        include("TestUnit.jl")
    end
    #doctest(Muscade)
end
;