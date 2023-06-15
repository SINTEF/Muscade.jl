using Muscade
module Runtest
    using Test,Documenter,Muscade
    @testset "Muscade.jl package" begin
        @testset "TestEspy" begin
            include("TestEspy.jl")
        end
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
        @testset "TestStaticXUA" begin
            include("TestStaticXUA.jl")
        end
        @testset "TestScale" begin
            include("TestScale.jl")
        end
        @testset "TestConstraints" begin
            include("TestConstraints.jl")
        end
        @testset "ElementCost" begin
            include("TestElementCost.jl")
        end
        @testset "TestUnit" begin
            include("TestUnit.jl")
        end
        doctest(Muscade)
    end
end