module TestDirectXUA_UQ

using Test
using Muscade
using StaticArrays,SparseArrays

# Classic spring element, with unknown force u
struct El1 <: AbstractElement
    K :: 𝕣
end
El1(nod::Vector{Node};K::𝕣) = El1(K)
@espy function Muscade.residual(o::El1, X,U,A, t,SP,dbg) 
    x,u = ∂0(X)[1], ∂0(U)[1]
    r         = -u + o.K*x
    return SVector(r),noFB
end
Muscade.doflist( ::Type{El1})  = (inod =(1 ,1), class=(:X,:U), field=(:tx1,:utx1))

# Parameters
model           = Model(:NoConfidenceX)
model2           = Model(:SameConfidenceXU)
k = 5. # Stiffness (N/m)
σₓ = 1e10
σᵤ = 0.1
σₓ₂ = 2e-2

n1              = addnode!(model,𝕣[0])
e1              = addelement!(model,El1,[n1], K=k)
n2              = addnode!(model2,𝕣[0])
e2              = addelement!(model2,El1,[n2], K=k)

# Costs
@functor with() fu(u,t)   = 0.5 * (u / σᵤ)^2
@functor with() l1(x,t) = 0.5 * ((x-1.) / σₓ)^2
@functor with() l2(x,t) = 0.5 * ((x-1.) / σₓ₂)^2
e11            = addelement!(model,SingleDofCost,[n1];class=:X,field=:tx1,     cost=l1)
e12             = addelement!(model,SingleDofCost,[n1];class=:U,field=:utx1,     cost=fu)
e21            = addelement!(model2,SingleDofCost,[n2];class=:X,field=:tx1,     cost=l2)
e22             = addelement!(model2,SingleDofCost,[n2];class=:U,field=:utx1,     cost=fu)

# Resolutions
initialstate    = initialize!(model)
stateX          = solve(SweepX{0};  initialstate=initialstate,time=[0.],verbose=false)
output        = solve(DirectXUA{0,0,0};primerstate=[stateX[1]],time = [1:0.5:1],verbose=false,maxiter=60, saveLvv=true)
stateXUA = output.state
Lvv_out = output.Lvv
covarMatrix = Matrix(Lvv_out)^-1

initialstate2    = initialize!(model2)
stateX2          = solve(SweepX{0};  initialstate=initialstate2,time=[0.],verbose=false)
output2        = solve(DirectXUA{0,0,0};primerstate=[stateX2[1]],time = [1:0.5:1],verbose=false,maxiter=60, saveLvv=true)
stateXUA2 = output2.state
Lvv_out2 = output2.Lvv
covarMatrix2 = Matrix(Lvv_out2)^-1

@testset "uncertain X" begin
    # This test considers a standard deviation of σₓ = 1e10 m (infinite) for the position x, but a std for the force u of σᵤ = 0.01 N.
    # It can be then compared to a analytical problem where you express the probability distribution of x, knowing the Gaussian pdf of u with mean 0 and variance σᵤ² = 0.0001.
    # The solution is analytically derived: x follows a Gaussian pdf with mean 0 and standard deviation σₓ = σᵤ/K. 
    @test covarMatrix[2,2] ≈ (σᵤ/k)^2
    @test covarMatrix[3,3] ≈ (σᵤ)^2
end

@testset "equally certain X and U" begin
    @test covarMatrix2[2,2] ≈ (σₓ₂)^2 * 0.5
    @test covarMatrix2[3,3] ≈ (σᵤ)^2  * 0.5
end

# display(Lvv_out)
# println("\n1/σₓ² = ", 1/σₓ^2)
# println("1/σᵤ² = ", 1/σᵤ^2)
# println("-----\n")
# covarMatrix = Matrix(Lvv_out)^-1
# display(covarMatrix)
# println("\nσₓ² = ", σₓ^2)
# println("σᵤ² = ", σᵤ^2)
# println("(σᵤ/K)² = ", (σᵤ/k)^2)

end 