using Test,Documenter,Muscade
println("\nMuscade test suite\n")

@testset "Muscade.jl package" begin
    @testset "SomeElements" begin
        include("TestSomeElements.jl")
        include("TestAssemble.jl")
    end
    @testset "TestUnit" begin
        include("TestUnit.jl")
    end
    #doctest(Muscade)
end
;