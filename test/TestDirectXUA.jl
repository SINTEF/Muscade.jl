#=
Verify that UU, XX XA and UA blocks are as expected (and set regression test for ΛX block)
Move on to solving this and line search
Profile and optimise
=#

# module TestDirectXUA

# cd("C:\\Users\\philippem\\.julia\\dev\\Muscade")
# using Pkg 
# Pkg.activate(".")



using Test
using Muscade
using StaticArrays,SparseArrays

###

using StaticArrays
struct El1 <: AbstractElement
    K :: 𝕣
    C :: 𝕣
    M :: 𝕣
end
El1(nod::Vector{Node};K::𝕣,C::𝕣,M::𝕣) = El1(K,C,M)
@espy function Muscade.residual(o::El1, X,U,A, t,SP,dbg) 
    x,x′,x″,u,ΞC,ΞM = ∂0(X)[1], ∂1(X)[1], ∂2(X)[1], ∂0(U)[1], A[1], A[2]
    r         = -u + o.K*x + o.C*exp10(ΞC)*x′ + o.M*exp10(ΞM)*x″
    return SVector(r),noFB
end
Muscade.doflist( ::Type{El1})  = (inod =(1 ,1 ,1 ,1), class=(:X,:U,:A,:A), field=(:tx1,:u,:ΞC,:ΞM))

include("SomeElements.jl")

model1          = Model(:TrueModel)
n1              = addnode!(model1,𝕣[0])  
n2              = addnode!(model1,𝕣[1])  
n3              = addnode!(model1,𝕣[]) # anode for spring
e1              = addelement!(model1,El1,[n1], K=1.,C=0.05,M=1.)
e2              = addelement!(model1,El1,[n2], K=0.,C=0.0,M=1.)
e3              = addelement!(model1,Spring{1},[n1,n2,n3], EI=1.1)
e4              = addelement!(model1,SingleDofCost,[n3];class=:A,field=:ΞL₀,cost=a->a^2)
e5              = addelement!(model1,SingleDofCost,[n3];class=:A,field=:ΞEI,cost=a->a^2)
e6              = addelement!(model1,SingleDofCost,[n1];class=:A,field=:ΞC ,cost=a->a^2)
e7              = addelement!(model1,SingleDofCost,[n1];class=:A,field=:ΞM ,cost=a->a^2)
e8              = addelement!(model1,SingleDofCost,[n2];class=:A,field=:ΞC ,cost=a->a^2)
e9              = addelement!(model1,SingleDofCost,[n2];class=:A,field=:ΞM ,cost=a->a^2)
e10             = addelement!(model1,SingleUdof   ,[n1];Xfield=:tx1,Ufield=:utx1,cost=u->u^2)
e11             = addelement!(model1,SingleUdof   ,[n2];Xfield=:tx1,Ufield=:utx1,cost=u->u^2)
e12             = addelement!(model1,SingleDofCost,[n1];class=:X,field=:tx1,cost=(tx1,t)->(tx1-0.1*sin(t))^2)
e13             = addelement!(model1,SingleDofCost,[n2];class=:X,field=:tx1,cost=(tx1,t)->(tx1-0.1*cos(t))^2)

state0          = initialize!(model1;nXder=2,time=0.)  # make space for 1st-order derivatives, 

dis             = state0.dis
out1,asm1       = Muscade.prepare(Muscade.AssemblyDirect    ,model1,dis,(1,3,1,1))
out2,asm2       = Muscade.prepare(Muscade.AssemblyDirectLine,model1)

zero!(out1)
zero!(out2)

s2 = sparse([1, 2, 1, 2], [1, 1, 2, 2], [0.0, 0.0, 0.0, 0.0], 2, 2)

Muscade.assemble!(out1,asm1,dis,model1,state0,(;))


nstep            = 6
ND               = 3
NA               = 1
Lv,Lvv,bigasm    = Muscade.preparebig(ND,NA,nstep,out1)

state            = [Muscade.State{1,ND,ND,@NamedTuple{γ::Float64}}(copy(state0)) for i = 1:nstep]
γ = 9.
Δt = 1.

Muscade.assemblebig!(Lvv,Lv,bigasm,asm1,model1,dis,out1,state,nstep,Δt,NA,γ,(caller=:TestDirectXUA,))

# using Spy
# fig = spy(Lvv,title="bigsparse Lvv sparsity",size=1000)
# # save("C:\\Users\\philippem\\.julia\\dev\\Muscade\\spy.jpg",fig)

stateXUA         = solve(DirectXUA{NA,ND};initialstate=state0,time=0:1.:5)

@testset "prepare_out" begin
    @test out1.L1[1] ≈ [[0.0, 0.0]]
    @test out1.L1[2] ≈ [[0.0, -0.2], [0.0, 0.0], [0.0, 0.0]]
    @test out1.L1[3] ≈ [[0.0, 0.0, 0.0, 0.0]]
    @test out1.L1[4] ≈ [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
    @test out1.L2[1,1] ≈ fill!(Matrix{Any}(undef,1,1),s2)
    @test out1.L2[2,2][1,1] ≈ sparse([1, 2, 1, 2], [1, 1, 2, 2], [2.0, 0.0, 0.0, 2.0], 2, 2)  
    @test out1.L2[2,2][1,2] ≈ sparse([1, 2, 1, 2], [1, 1, 2, 2], [0.0, 0.0, 0.0, 0.0], 2, 2)  
    @test out1.L2[2,2][1,3] ≈ sparse([1, 2, 1, 2], [1, 1, 2, 2], [0.0, 0.0, 0.0, 0.0], 2, 2)
    @test out1.L2[2,1][1,1] ≈ sparse([1, 2, 1, 2], [1, 1, 2, 2], [2.1, -1.1, -1.1, 1.1], 2, 2)
    @test out1.L2[2,1][2,1] ≈ sparse([1, 2, 1, 2], [1, 1, 2, 2], [0.05, 0.0, 0.0, 0.0], 2, 2)
    @test out1.L2[2,1][3,1] ≈ sparse([1, 2, 1, 2], [1, 1, 2, 2], [1.0, 0.0, 0.0, 1.0], 2, 2)
    @test out1.L2[3,3][1,1] ≈ sparse([1, 2, 3, 4], [1, 2, 3, 4], [0.0, 0.0, 2.0, 2.0], 4, 4)
    @test out1.L2[3,4][1,1] ≈ sparse([1, 1, 2, 2], [1, 2, 3, 4], [0.0, 0.0, 0.0, 0.0], 4, 6)
    @test out1.L2[4,4][1,1] ≈ sparse([1, 2, 1, 2, 3, 4, 3, 4, 5, 6, 5, 6], [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6], [2.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 2.0], 6, 6)
    @test typeof(out2) ==  Muscade.AssemblyDirectLine
    @test out2.ming ≈ Inf
    @test out2.minλ ≈ Inf
    @test out2.Σλg ≈ 0.
    @test out2.npos ≈ 0
end
@testset "prepare_asm" begin
    @test asm1[1,1]  ≈ [1 2]      # asm1[iarray,ieletyp][ieledof,iele] -> idof|inz
    @test asm1[1,2]  ≈ [1; 2;;]
    @test asm1[2,1]  ≈ [1 2] 
    @test asm1[2,2]  ≈ [1; 2;;]
    @test asm1[3,1]  ≈ [1 2] 
    @test asm1[3,2]  ≈ Matrix{Int64}(undef, 0, 1)
    @test asm1[4,1]  ≈ [1 3; 2 4]
    @test asm1[4,2]  ≈ [5; 6;;]
    @test asm1[5,1]  ≈ [1 4]                   
    @test asm1[20,1] ≈ [1 5; 2 6; 3 7; 4 8]    
end

@testset "preparebig Lvv" begin
    @test size(Lv)         == (54,)
    @test size(Lvv)        == (54,54)
    @test Lvv.colptr       == [1,29,57,85,113,127,141,153,165,193,221,249,277,291,305,317,329,369,409,449,489,509,529,547,565,605,645,685,725,745,765,783,801,829,857,885,913,927,941,953,965,993,1021,1049,1077,1091,1105,1117,1129,1149,1169,1189,1209,1235,1261]
    @test Lvv.rowval[1:60] == [1,2,3,4,5,7,9,10,11,12,13,15,17,18,19,20,21,23,25,26,27,28,29,31,49,50,53,54,1,2,3,4,6,8,9,10,11,12,14,16,17,18,19,20,22,24,25,26,27,28,30,32,51,52,53,54,1,2,3,4]
end

@testset "preparebig ,bigasm" begin
    @test bigasm.pgc      == [1, 3, 5, 9, 11, 13, 17, 19, 21, 25, 27, 29, 33, 35, 37, 41, 43, 45, 49, 55]
    @test bigasm.pigr[1]' == [  1   3   5   7   9  11  13  15  17  19  21  23  0  0  0  0  0  0  25;    29  31  33  35  37  39  41  43  45  47  49  51  0  0  0  0  0  0  53]
end



#end 

;