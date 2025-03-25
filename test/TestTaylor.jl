#module TestTaylor
using Muscade
using Test,StaticArrays


f(X) = SVector(cos(X[1])*X[3],sin(X[2])*X[3])
X₀   = SVector(0.,0.,1.)
X₁   = SVector(1.,1.,1.)

te1 = Taylor{1}(f,X₀)
t1  = te1(X₁)
te2 = Taylor{2}(f,X₀)
t2  = te2(X₁)

@testset "Taylor" begin
    @test typeof(te1) == Taylor{1, 3, Tuple{SVector{2, Float64}, SMatrix{2, 3, Float64, 6}}}
    @test te1.x ≈ X₀
    @test te1.A[1] ≈ [1,0]
    @test te1.A[2] ≈ [0 0 1;0 1 0]
    @test t1 ≈ [1,1]

    @test typeof(te2) == Taylor{2, 3, Tuple{SVector{2, Float64}, SMatrix{2, 3, Float64, 6}, SArray{Tuple{2, 3, 3}, Float64, 3, 18}}}
    @test te2.x ≈ X₀
    @test te2.A[1] ≈ [1,0]
    @test te2.A[2] ≈ [0 0 1;0 1 0]
    @test te2.A[3] ≈ [-.5;0 ;; 0;0 ;; 0;0 ;;; 0;0 ;; 0;0 ;; 0;.5 ;;; 0;0 ;; 0;.5 ;; 0;0]
    @test t2 ≈ [.5,1]
end

#end