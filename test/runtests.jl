module Runtest
    using Test,Muscade

    @testset "Muscade.jl package" begin
        @testset "TestEspy" begin
            include("TestEspy.jl")
        end
        @testset "TestElementAPI" begin
            include("TestElementAPI.jl")
        end
        @testset "TestDots" begin
            include("TestDots.jl")
        end
        @testset "TestAdiff" begin
            include("TestAdiff.jl")
        end
        @testset "TestTaylor" begin
            include("TestTaylor.jl")
        end
        @testset "TestFunctors" begin
            include("TestFunctors.jl")
        end
        @testset "TestModelDescription" begin
            include("TestModelDescription.jl")
        end
        @testset "TestAssemble" begin
            include("TestAssemble.jl")
        end
        @testset "TestSweepX2" begin
            include("TestSweepX2.jl")
        end
        @testset "TestSweepX1" begin
            include("TestSweepX1.jl")
        end
        @testset "TestSweepX0" begin
            include("TestSweepX0.jl")
        end
        #@testset "TestNewmarkSweep" begin
        #    include("TestNewmarkSweep.jl")
        #end
        @testset "TestDirectXUA" begin
            include("TestDirectXUA.jl")
        end
        @testset "TestDirectXUA001" begin
            include("TestDirectXUA001.jl")
        end
        @testset "TestFreqXU" begin
            include("TestFreqXU.jl")
        end
        @testset "TestScale" begin
            include("TestScale.jl")
        end
        @testset "TestElementCost" begin
            include("TestElementCost.jl")
        end
        @testset "TestRotations" begin
            include("TestRotations.jl")
        end
        @testset "TestUnit" begin
            include("TestUnit.jl")
        end
        @testset "TestSparseTools" begin
            include("TestSparseTools.jl")
        end
        @testset "TestFFT" begin
            include("TestFFT.jl")
        end
        @testset "TestEigenmodes" begin
            include("TestEigenmodes.jl")
        end
        @testset "TestEigX" begin
            include("TestEigX.jl")
        end
        @testset "TestFiniteDifferences" begin
            include("TestFiniteDifferences.jl")
        end
        @testset "TestElementTestTools" begin
            include("TestElementTestTools.jl")
        end
        @testset "TestBarElement" begin
            include("TestBarElement.jl")
        end
        @testset "TestBeamElement" begin
            include("TestBeamElement.jl")
        end
        @testset "TestBeamElementStrainGauge" begin
            include("TestBeamElementStrainGauge.jl")
        end
        @testset "TestDrawBeamElement" begin
            include("TestDrawBeamElement.jl")
        end
        @testset "TestPositionElement" begin
            include("TestPositionElement.jl")
        end
    end
end
