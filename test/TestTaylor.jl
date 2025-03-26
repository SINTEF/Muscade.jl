module TestTaylor
using Muscade
using Test,StaticArrays


f(X) = SVector(cos(X[1])*X[3],sin(X[2])*X[3])
X₀   = SVector(0.,0.,1.)
X₁   = SVector(1.,1.,1.)

te1 = Taylor{1}(f,X₀)
t1  = te1(X₁)
te2 = Taylor{2}(f,X₀)
t2  = te2(X₁)

te2′ = ∂(te2)
t2′  = te2′(X₁)

@testset "Taylor" begin
    @test typeof(te1) == Taylor{1, 3, SVector{2, ∂ℝ{1, 3, Float64}}}

    @test te1.x ≈ X₀
    @test value{1}(te1.y) ≈ [1,0]
    @test ∂{1,3}(te1.y) ≈ [0 0 1;0 1 0]
    @test t1 ≈ [1,1]

    @test typeof(te2) ==  Taylor{2, 3, SVector{2, ∂ℝ{2, 3, ∂ℝ{1, 3, Float64}}}}
    @test te2.x ≈ X₀
    @test value{1}(value{2}(te2.y)) ≈ [1,0]
    @test ∂{1,3}(value{2}(te2.y)) ≈ [0 0 1;0 1 0]
    @test ∂{1,3}(∂{2,3}(te2.y)) ≈ [-1.0 0.0 0.0; 0.0 0.0 0.0;;; 0.0 0.0 0.0; 0.0 0.0 1.0;;; 0.0 0.0 0.0; 0.0 1.0 0.0]
    @test t2 ≈ [.5,1]

    @test typeof(te2′) == Taylor{1, 3, SMatrix{2, 3, ∂ℝ{1, 3, Float64}, 6}}
    @test te2′.x ≈ X₀
    @test value{1}(te2′.y) ≈ [0 0 1;0 1 0]
    @test ∂{1,3}(te2′.y) ≈ [-1.0 0.0 0.0; 0.0 0.0 0.0;;; 0.0 0.0 0.0; 0.0 0.0 1.0;;; 0.0 0.0 0.0; 0.0 1.0 0.0]
    @test t2′ ≈ [-1 0 1;0 1 1]
end

end