test    = @__DIR__
muscade = normpath(joinpath(test,".."))
Pkg.activate(test)

module Runtest
    using Test,Literate, DocumenterCitations,Printf,Documenter,Muscade

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
        @testset "TestMultiplex" begin
            include("TestMultiplex.jl")
        end
        @testset "TestModelDescription" begin
            include("TestModelDescription.jl")
        end
        @testset "TestAssemble" begin
            include("TestAssemble.jl")
        end
        @testset "TestNewmarkX" begin
            include("TestNewmarkX.jl")
        end
        @testset "TestStaticX" begin
            include("TestStaticX.jl")
        end
        @testset "TestStaticXUA" begin
            include("TestStaticXUA.jl")
        end
        @testset "TestDirectXUA" begin
            include("TestDirectXUA.jl")
        end
        @testset "TestStaticXUAwithineq" begin
            include("TestStaticXUAwithineq.jl")
        end
        @testset "TestScale" begin
            include("TestScale.jl")
        end
        @testset "TestDofConstraints" begin
            include("TestDofConstraints.jl")
        end
        @testset "ElementCost" begin
            include("TestElementCost.jl")
        end
        @testset "BeamElement" begin
            include("TestBeamElement.jl")
        end
        @testset "Rotations" begin
            include("TestRotations.jl")
        end
        @testset "TestUnit" begin
            include("TestUnit.jl")
        end
        @testset "TestBlockSparse" begin
            include("TestBlockSparse.jl")
        end
        @testset "TestFiniteDifferences" begin
            include("TestFiniteDifferences.jl")
        end
        doctest(Muscade)
    end
end

Pkg.activate(muscade) 
