using Test,Documenter,Muscade
println("\nMuscade test suite\n")

@testset "Muscade.jl package" begin
    @testset "SomeElements" begin
        include("TestElementAPI.jl")
        include("TestModelDescription.jl")
    end
    @testset "TestUnit" begin
        include("TestUnit.jl")
    end
    #doctest(Muscade)
end
;