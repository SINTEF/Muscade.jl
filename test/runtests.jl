using Test,Documenter,Muscade
println("\nMuscade test suite\n")

@testset "Muscade.jl package" begin
    include("TestOptimist.jl")
    #doctest(Muscade)
end
;