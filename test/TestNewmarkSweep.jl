#module TestNewmarkSweep
using Muscade,Test,StaticArrays

K  = 1.
C  = .4
M  = .3
Re = 0.
f(x,x′,x″) = K*x+C*x′+M*x″ -Re

Δt  = 0.01
β,γ = 1/4,1/2
a₁,a₂,a₃,b₁,b₂,b₃ = γ/(β*Δt),γ/β ,(γ/2β-1)*Δt,1/(β*Δt^2),1/(β*Δt),1/2β

#n  = 3
n  = 1000
x  = randn(n) # randomness gets overwriten
x′ = randn(n)
x″ = randn(n)
# y  = randn(n)  # randomness corrected in calculations (in this linear problem!)
# y′ = randn(n)
# y″ = randn(n)
y  = zeros(n)  # randomness corrected in calculations (in this linear problem!)
y′ = zeros(n)
y″ = zeros(n)

x[ 1] = y[ 1] = 0.
x′[1] = y′[1] = 1.
x″[1] = y″[1] = 0.

# std Newmark-β step
for i = 2:n
    i⁻         = i-1
    δX         = Muscade.∂ℝ{1,2,𝕣}(0.,SVector(1.,0.))
    δr         = Muscade.∂ℝ{1,2,𝕣}(0.,SVector(0.,1.))
    a          = a₂*x′[i⁻] + a₃*x″[i⁻]
    b          = b₂*x′[i⁻] + b₃*x″[i⁻]
    vx         = x[ i⁻] +    δX
    vx′        = x′[i⁻] + a₁*δX + a*δr 
    vx″        = x″[i⁻] + b₁*δX + b*δr
    vr         = f(vx,vx′,vx″)
    r          = value{1}(vr)
    B          = r - ∂{1,2}(vr)[2]
    A          = ∂{1,2}(vr)[1] 
    Δx         = -A\B
    Δx′        = a₁*Δx - a
    Δx″        = b₁*Δx - b
    x[i]       = x[i⁻] + Δx
    x′[i]      = x′[i⁻] + Δx′
    x″[i]      = x″[i⁻] + Δx″
end

# Newmark-β sweep 
for i = 2:n
    i⁻         = i-1
    #@show y[i], y′[i], y″[i], y[i⁻], y′[i⁻], y″[i⁻]
    δX         = Muscade.∂ℝ{1,2,𝕣}(0.,SVector(1.,0.))
    δr         = Muscade.∂ℝ{1,2,𝕣}(0.,SVector(0.,1.))
    a          = a₁*(y[i⁻].-y[i]) + (a₂-1)*y′[i⁻] +     a₃*y″[i⁻] + y′[i]      
    b          = b₁*(y[i⁻].-y[i]) +     b₂*y′[i⁻] + (b₃-1)*y″[i⁻] + y″[i]       
    vy         = y[ i] +    δX
    vy′        = y′[i] + a₁*δX + a*δr 
    vy″        = y″[i] + b₁*δX + b*δr
    vr         = f(vy,vy′,vy″)
    r          = value{1}(vr)
    B          = r - ∂{1,2}(vr)[2]
    A          = ∂{1,2}(vr)[1] 
    dy         = -A\B
    #@show A,B
    dy′        = a₁*dy - a
    dy″        = b₁*dy - b
    y[i]       = y[i]  + dy    
    y′[i]      = y′[i] + dy′  
    y″[i]      = y″[i] + dy″  
end

@testset "Newmark sweep" begin
    @test x  ≈ y   
    @test x′ ≈ y′
    @test x″ ≈ y″
end

using GLMakie
fig      = Figure(size = (1000,750))
axe      = Axis(fig[1,1],title="Test",xlabel="time",ylabel="x")
display(fig)
lines!(  axe,Δt*(1:n),x , linewidth = 1,color=:black)
lines!(  axe,Δt*(1:n),y , linewidth = 1,color=:red  )

#end



