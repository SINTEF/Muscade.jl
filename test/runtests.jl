using Test,Documenter,Muscade
println("\nMuscade test suite\n")

@testset "Muscade.jl package" begin
    include("TestSomeElements.jl")
    include("TestUnit.jl")
    #doctest(Muscade)
end
