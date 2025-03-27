module TestTaylor
using Muscade
using Test,StaticArrays


f(X) = SVector(cos(X[1])*X[3],sin(X[2])*X[3])
g(X) = 4.5*X[1]^2/2 + X[1] + 2.
h(X) = SVector(3.,4.)
k(X) = 5.

w(X) = (f(X),g(X),h(X),k(X))
X₀   = SVector(0.,0.,1.)
X₁   = SVector(1.,1.,1.)

tf0,tg0,th0,tk0 = Taylor{0}(w,X₀)
tf1,tg1,th1,tk1 = Taylor{1}(w,X₀)
tf2,tg2,th2,tk2 = Taylor{2}(w,X₀)
tf2′            = ∂(tf2)

f0  = tf0(X₁)
g0  = tg0(X₁)
h0  = th0(X₁)
k0  = tk0(X₁)

f1  = tf1(X₁)
g1  = tg1(X₁)
h1  = th1(X₁)
k1  = tk1(X₁)

f2  = tf2(X₁)
g2  = tg2(X₁)
h2  = th2(X₁)
k2  = tk1(X₁)
f2′ = tf2′(X₁)

@testset "Taylor" begin
    @test typeof(tf1) == Taylor{1, 3, SVector{2, ∂ℝ{1, 3, Float64}}}

    @test tf1.x ≈ X₀
    @test tf1.y ≈ [∂ℝ{1, 3, Float64}(1.0, [0.0, 0.0, 1.0]),    ∂ℝ{1, 3, Float64}(0.0, [0.0, 1.0, 0.0])]
    @test f1 ≈ [1,1]

    @test typeof(tf2) ==  Taylor{2, 3, SVector{2, ∂ℝ{2, 3, ∂ℝ{1, 3, Float64}}}}
    @test tf2.x ≈ X₀
    @test tf2.y ≈ [∂ℝ{2, 3, ∂ℝ{1, 3, Float64}}(∂ℝ{1, 3, Float64}(1.0, [0.0, 0.0, 1.0]), ∂ℝ{1, 3, Float64}[∂ℝ{1, 3, Float64}(0.0, [-1.0, 0.0, 0.0]), ∂ℝ{1, 3, Float64}(0.0, [0.0, 0.0, 0.0]), ∂ℝ{1, 3, Float64}(1.0, [0.0, 0.0, 0.0])]),
                   ∂ℝ{2, 3, ∂ℝ{1, 3, Float64}}(∂ℝ{1, 3, Float64}(0.0, [0.0, 1.0, 0.0]), ∂ℝ{1, 3, Float64}[∂ℝ{1, 3, Float64}(0.0, [ 0.0, 0.0, 0.0]), ∂ℝ{1, 3, Float64}(1.0, [0.0, 0.0, 1.0]), ∂ℝ{1, 3, Float64}(0.0, [0.0, 1.0, 0.0])])]
    @test f2 ≈ [.5,1]

    @test typeof(tf2′) == Taylor{1, 3, SMatrix{2, 3, ∂ℝ{1, 3, Float64}, 6}}
    @test tf2′.x ≈ X₀
    @test tf2′.y ≈ [∂ℝ{1, 3, Float64}(0.0, [-1.0, 0.0, 0.0])  ∂ℝ{1, 3, Float64}(0.0, [0.0, 0.0, 0.0])  ∂ℝ{1, 3, Float64}(1.0, [0.0, 0.0, 0.0]);
                    ∂ℝ{1, 3, Float64}(0.0, [0.0, 0.0, 0.0])   ∂ℝ{1, 3, Float64}(1.0, [0.0, 0.0, 1.0])  ∂ℝ{1, 3, Float64}(0.0, [0.0, 1.0, 0.0])]
    @test f2′ ≈ [-1 0 1;0 1 1]

    @test tg0.y ≈ 2.0
    @test g0  ≈ 2.
    @test tg1.y ≈ ∂ℝ{1, 3, Float64}(2.0, [1.0, 0.0, 0.0])
    @test g1  ≈ 3.
    @test tg2.y ≈ ∂ℝ{2, 3, ∂ℝ{1, 3, Float64}}(∂ℝ{1, 3, Float64}(2.0, [1.0, 0.0, 0.0]), ∂ℝ{1, 3, Float64}[∂ℝ{1, 3, Float64}(1.0, [4.5, 0.0, 0.0]), ∂ℝ{1, 3, Float64}(0.0, [0.0, 0.0, 0.0]), ∂ℝ{1, 3, Float64}(0.0, [0.0, 0.0, 0.0])])
    @test g2  ≈ 5.25

    @test th0.y ≈ [3.0, 4.0]
    @test th1.y ≈ [3.0, 4.0]
    @test th2.y ≈ [3.0, 4.0]
    @test h0    ≈ [3.0, 4.0]
    @test h1    ≈ [3.0, 4.0]
    @test h2    ≈ [3.0, 4.0]

    @test tk0.y ≈ 5.
    @test tk1.y ≈ 5.
    @test tk2.y ≈ 5.
    @test k0    ≈ 5.
    @test k1    ≈ 5.
    @test k2    ≈ 5.

end
end # module