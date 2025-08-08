module TestTaylor
using Muscade
using Test,StaticArrays

const X1 = (SVector{3,𝕣}(1,2,3),)
const X2 = (SVector{3,𝕣}(1,2,3),SVector{3,𝕣}(4,5,6))
const X3 = (SVector{3,𝕣}(1,2,3),SVector{3,𝕣}(4,5,6),SVector{3,𝕣}(7,8,9))
const X4 = (variate{1,3}(SVector{3,𝕣}(1,2,3)),)
const X5 = (variate{1,3}(SVector{3,𝕣}(1,2,3)),variate{1,3}(SVector{3,𝕣}(4,5,6)))
const X6 = (variate{1,3}(SVector{3,𝕣}(1,2,3)),variate{1,3}(SVector{3,𝕣}(4,5,6)),variate{1,3}(SVector{3,𝕣}(7,8,9)))

Y1=motion{2}(X1)
Y2=motion{2}(X2)
Y3=motion{2}(X3)
Y4=motion{3}(X4)
Y5=motion{3}(X5)
Y6=motion{3}(X6)

@testset "motion" begin
    @test typeof(X1) == Tuple{SVector{3, Float64}}
    @test typeof(Y1) == SVector{3, Float64}

    @test typeof(X2) == Tuple{SVector{3, Float64}, SVector{3, Float64}}
    @test typeof(Y2) == SVector{3, ∂ℝ{2, 1, Float64}} 

    @test typeof(X3) == Tuple{SVector{3, Float64}, SVector{3, Float64}, SVector{3, Float64}}
    @test typeof(Y3) == SVector{3, ∂ℝ{3, 1, ∂ℝ{2, 1, Float64}}}

    @test typeof(X4) == Tuple{SVector{3, ∂ℝ{1, 3, Float64}}}
    @test typeof(Y4) == SVector{3, ∂ℝ{1, 3, Float64}}

    @test typeof(X5) == Tuple{SVector{3, ∂ℝ{1, 3, Float64}}, SVector{3, ∂ℝ{1, 3, Float64}}}
    @test typeof(Y5) == SVector{3, ∂ℝ{3, 1, ∂ℝ{1, 3, Float64}}}

    @test typeof(X6) == Tuple{SVector{3, ∂ℝ{1, 3, Float64}}, SVector{3, ∂ℝ{1, 3, Float64}}, SVector{3, ∂ℝ{1, 3, Float64}}}
    @test typeof(Y6) == SVector{3, ∂ℝ{4, 1, ∂ℝ{3, 1, ∂ℝ{1, 3, Float64}}}}
end

@testset "motion⁻¹" begin
    @test motion⁻¹{2,1,0}(Y1) === SVector{3,𝕣}(1,2,3)
    @test motion⁻¹{2,1,1}(Y1) === SVector{3,𝕣}(0,0,0)
    @test motion⁻¹{2,1,2}(Y1) === SVector{3,𝕣}(0,0,0)

    @test motion⁻¹{2,2,0}(Y2) === SVector{3,𝕣}(1,2,3)
    @test motion⁻¹{2,2,1}(Y2) === SVector{3,𝕣}(4,5,6)
    @test motion⁻¹{2,2,2}(Y2) === SVector{3,𝕣}(0,0,0)

    @test motion⁻¹{2,3,0}(Y3) === SVector{3,𝕣}(1,2,3)
    @test motion⁻¹{2,3,1}(Y3) === SVector{3,𝕣}(4,5,6)
    @test motion⁻¹{2,3,2}(Y3) === SVector{3,𝕣}(7,8,9)

    @test motion⁻¹{3,1,0}(Y4) === variate{1,3}(SVector{3,𝕣}(1,2,3))
    @test motion⁻¹{3,1,1}(Y4) === SVector{3,𝕣}(0,0,0)
    @test motion⁻¹{3,1,2}(Y4) === SVector{3,𝕣}(0,0,0)

    @test motion⁻¹{3,2,0}(Y5) === variate{1,3}(SVector{3,𝕣}(1,2,3))
    @test motion⁻¹{3,2,1}(Y5) === variate{1,3}(SVector{3,𝕣}(4,5,6))
    @test motion⁻¹{3,2,2}(Y5) === SVector{3,𝕣}(0,0,0)

    @test motion⁻¹{3,3,0}(Y6) === variate{1,3}(SVector{3,𝕣}(1,2,3))
    @test motion⁻¹{3,3,1}(Y6) === variate{1,3}(SVector{3,𝕣}(4,5,6))
    @test motion⁻¹{3,3,2}(Y6) === variate{1,3}(SVector{3,𝕣}(7,8,9))

    @test motion⁻¹{2,1}(Y1) === (SVector{3,𝕣}(1,2,3),)
    @test motion⁻¹{2,2}(Y2) === (SVector{3,𝕣}(1,2,3),SVector{3,𝕣}(4,5,6))
    @test motion⁻¹{2,3}(Y3) === (SVector{3,𝕣}(1,2,3),SVector{3,𝕣}(4,5,6),SVector{3,𝕣}(7,8,9))
end

a = SVector(3.,4.)
@testset "revariate" begin
    @test Muscade.revariate{0}(a)    === a
    @test Muscade.revariate{1}(a)[1] ===                                                                   ∂ℝ{1, 2, Float64}(3.0, [1.0, 0.0])
    @test Muscade.revariate{2}(a)[1] ===                                       ∂ℝ{2, 2, ∂ℝ{1, 2, Float64}}(∂ℝ{1, 2, Float64}(3.0, [1.0, 0.0]), ∂ℝ{1, 2, Float64}[∂ℝ{1, 2, Float64}(1.0, [0.0, 0.0]), ∂ℝ{1, 2, Float64}(0.0, [0.0, 0.0])])
    @test Muscade.revariate{3}(a)[1] === ∂ℝ{3, 2, ∂ℝ{2, 2, ∂ℝ{1, 2, Float64}}}(∂ℝ{2, 2, ∂ℝ{1, 2, Float64}}(∂ℝ{1, 2, Float64}(3.0, [1.0, 0.0]), ∂ℝ{1, 2, Float64}[∂ℝ{1, 2, Float64}(1.0, [0.0, 0.0]), ∂ℝ{1, 2, Float64}(0.0, [0.0, 0.0])]), ∂ℝ{2, 2, ∂ℝ{1, 2, Float64}}[∂ℝ{2, 2, ∂ℝ{1, 2, Float64}}(∂ℝ{1, 2, Float64}(1.0, [0.0, 0.0]), ∂ℝ{1, 2, Float64}[∂ℝ{1, 2, Float64}(0.0, [0.0, 0.0]), ∂ℝ{1, 2, Float64}(0.0, [0.0, 0.0])]), ∂ℝ{2, 2, ∂ℝ{1, 2, Float64}}(∂ℝ{1, 2, Float64}(0.0, [0.0, 0.0]), ∂ℝ{1, 2, Float64}[∂ℝ{1, 2, Float64}(0.0, [0.0, 0.0]), ∂ℝ{1, 2, Float64}(0.0, [0.0, 0.0])])])
end

f(X) = SVector(cos(X[1])*X[3],sin(X[2])*X[3])
g(X) = 4.5*X[1]^2/2 + X[1] + 2.
h(X) = SVector(3.,4.)
k(X) = 5.
w(X) = (f(X),g(X),h(X),k(X))
X₀   = SVector(0.,0.,1.)
vX₀  = variate{1,3}(X₀)

@testset "fast" begin
    @test fast(w, X₀) === w( X₀)
    @test fast(w,vX₀) === w(vX₀)
end

yy    = Muscade.firstorderonly(variate{2,3}(variate{1,3}(SVector{3,𝕣}(1,2,3))),
                               variate{2,3}(variate{1,3}(SVector{3,𝕣}(4,5,6))),
                               3.)
fooyy =                       (variate{1,3}(SVector{3,𝕣}(1,2,3)) ,
                               variate{1,3}(SVector{3,𝕣}(4,5,6)) ,
                               3.)                    
@testset "firstorderonly" begin
    @test yy === fooyy
end

end # module
