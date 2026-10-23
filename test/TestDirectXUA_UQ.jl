# module TestDirectXUA

using Test
using Muscade
using StaticArrays,SparseArrays

struct El1 <: AbstractElement
    K :: 𝕣
end
El1(nod::Vector{Node};K::𝕣) = El1(K)
@espy function Muscade.residual(o::El1, X,U,A, t,SP,dbg) 
    x,u = ∂0(X)[1], ∂0(U)[1]
    r         = -u + o.K*x
    return SVector(r),noFB
end
Muscade.doflist( ::Type{El1})  = (inod =(1 ,1), class=(:X,:U), field=(:tx1,:u))

model           = Model(:TrueModel)
n1              = addnode!(model,𝕣[0])
e1              = addelement!(model,El1,[n1], K=1.)
σₓ = 0.5
σᵤ = 0.1
@functor with() fu(u,t)   = 0.5 * (u/ σᵤ)^2
@functor with() l1(tx1,t) = 0.5 * ((tx1-0.1*sin(t))/σₓ)^2
e10             = addelement!(model,SingleUdof   ,[n1];Xfield=:tx1,Ufield=:utx1,cost=fu)
e12             = addelement!(model,SingleDofCost,[n1];class=:X,field=:tx1,     cost=l1)

initialstate    = initialize!(model)
stateX          = solve(SweepX{0};  initialstate,time=[0.],verbose=false)
output        = solve(DirectXUA{0,0,0};primerstate=[stateX[1]],time = [1:0.5:1],verbose=false,maxiter=50, saveLvv=true)
# output        = solve(DirectXUA{0,0,0};primerstate=[stateX[1]],time = [0:.1:1],verbose=false,maxiter=50, saveLvv=true)
stateXUA = output.state
Lvv_out = output.Lvv
display(Lvv_out)
println("1/σₓ² = ", 1/σₓ^2)
println("1/σᵤ² = ", 1/σᵤ^2)

# covarMatrix = Matrix(Lvv_out)^-1

# @testset "solution" begin
#     @test 1==1
# end

# @testset "saveLvv" begin
#     @test 1==1
# end

# @testset "UQstatic" begin
#     @test 1==1
# end

# end 