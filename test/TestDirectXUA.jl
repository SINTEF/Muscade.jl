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

### SDOF oscillator

# model1          = Model(:TrueModel)
# n               = addnode!(model1,𝕣[ 0, 0])  
# e               = addelement!(model1,El1,[n], K=1.,C=0.05,M=1.)
# state0          = initialize!(model1;nXder=2,time=0.)  # make space for 1st-order derivatives, 
# setdof!(state0,[1.];field=:x,nodID=[n],order=1)
# time            = 1.:.1:100
# state1          = solve(SweepX{2};  initialstate=state0,time,verbose=false)
# x               = getdof(state1;field=:x,nodID=[n],order=0 )
# x               = reshape(x,length(x))
# using GLMakie
# fig             = Figure(size = (1000,800))
# axe             = Axis(fig[1,1],title="Test",xlabel="time",ylabel="x")
# oedge           = lines!(  axe,time,x , linewidth = 1)

### Test out1, asm1

include("SomeElements.jl")

model1          = Model(:TrueModel)
n1              = addnode!(model1,𝕣[0])  
n2              = addnode!(model1,𝕣[1])  
n3              = addnode!(model1,𝕣[]) # anode
e1              = addelement!(model1,El1,[n1], K=1.,C=0.05,M=1.)
e2              = addelement!(model1,El1,[n2], K=0.,C=0.0,M=1.)
e3              = addelement!(model1,Spring{1},[n1,n2,n3], EI=1.1)
state0          = initialize!(model1;nXder=2,time=0.)  # make space for 1st-order derivatives, 

dis             = state0.dis
out1,asm1       = Muscade.prepare(Muscade.AssemblyDirect    ,model1,dis,(1,3,1,1))
out2,asm2       = Muscade.prepare(Muscade.AssemblyDirectLine,model1)

zero!(out1)
zero!(out2)
# Dofs of class :X
# 1. field= :tx1             NodID(1)
# 2. field= :tx1             NodID(2)
#
# Dofs of class :U
# 1. field= :u               NodID(1)
# 2. field= :u               NodID(2)
#
# Dofs of class :A
# 1. field= :ΞC              NodID(1)
# 2. field= :ΞM              NodID(1)
# 3. field= :ΞC              NodID(2)
# 4. field= :ΞM              NodID(2)
# 5. field= :ΞL₀             NodID(3)
# 6. field= :ΞEI             NodID(3)
s2 = sparse([1, 2, 1, 2], [1, 1, 2, 2], [0.0, 0.0, 0.0, 0.0], 2, 2)
@testset "prepare_out" begin
    @test out1.L1[1] ≈ [[0.0, 0.0]]
    @test out1.L1[2] ≈ [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
    @test out1.L1[3] ≈ [[0.0, 0.0]]
    @test out1.L1[4] ≈ [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
    @test out1.L2[1,1] ≈ fill!(Matrix{Any}(undef,1,1),s2)
    @test out1.L2[2,2] ≈ fill!(Matrix{Any}(undef,3,3),s2)
    @test out1.L2[1,2] ≈ fill!(Matrix{Any}(undef,1,3),s2)
    @test out1.L2[2,1] ≈ fill!(Matrix{Any}(undef,3,1),s2)
    @test out1.L2[3,3] ≈ fill!(Matrix{Any}(undef,1,1),s2)
    @test out1.L2[3,4] ≈ [sparse([1, 1, 2, 2], [1, 2, 3, 4], [0.0, 0.0, 0.0, 0.0], 2, 6)]
    @test out1.L2[4,4] ≈ [sparse([1, 2, 1, 2, 3, 4, 3, 4, 5, 6, 5, 6], [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 6, 6)]
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
Muscade.assemble!(out1,asm1,dis,model1,state0,(;))
@testset "prepare_out" begin
    @test out1.L1[1] ≈ [[0.0, 0.0]]
    @test out1.L1[2] ≈ [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
    @test out1.L1[3] ≈ [[0.0, 0.0]]
    @test out1.L1[4] ≈ [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
    @test out1.L2[1,1] ≈ fill!(Matrix{Any}(undef,1,1),s2)
    @test out1.L2[2,2] ≈ fill!(Matrix{Any}(undef,3,3),s2)
    @test out1.L2[1,2][1] ≈ sparse([1, 2, 1, 2], [1, 1, 2, 2], [2.1, -1.1, -1.1, 1.1], 2, 2)
    @test out1.L2[1,2][2] ≈ sparse([1, 2, 1, 2], [1, 1, 2, 2], [0.05, 0.0, 0.0, 0.0], 2, 2) 
    @test out1.L2[1,2][3] ≈ sparse([1, 2, 1, 2], [1, 1, 2, 2], [1.0, 0.0, 0.0, 1.0], 2, 2)
    @test out1.L2[2,1][1] ≈ sparse([1, 2, 1, 2], [1, 1, 2, 2], [2.1, -1.1, -1.1, 1.1], 2, 2)
    @test out1.L2[2,1][2] ≈ sparse([1, 2, 1, 2], [1, 1, 2, 2], [0.05, 0.0, 0.0, 0.0], 2, 2) 
    @test out1.L2[2,1][3] ≈ sparse([1, 2, 1, 2], [1, 1, 2, 2], [1.0, 0.0, 0.0, 1.0], 2, 2)
    @test out1.L2[3,1][1] ≈ sparse([1, 2], [1, 2], [-1.0, -1.0], 2, 2)
    @test out1.L2[1,3][1] ≈ sparse([1, 2], [1, 2], [-1.0, -1.0], 2, 2)
    @test out1.L2[3,3] ≈ fill!(Matrix{Any}(undef,1,1),s2)
    @test out1.L2[3,4] ≈ [sparse([1, 1, 2, 2], [1, 2, 3, 4], [0.0, 0.0, 0.0, 0.0], 2, 6)]
    @test out1.L2[4,4] ≈ [sparse([1, 2, 1, 2, 3, 4, 3, 4, 5, 6, 5, 6], [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 6, 6)]
    @test out1.L2[1,4][1] ≈ sparse([1, 1, 2, 2, 1, 2, 1, 2], [1, 2, 3, 4, 5, 5, 6, 6], [0.0, 0.0, 0.0, 0.0, 2.532843602293451, -2.532843602293451, 0.0, 0.0], 2, 6)
    @test out1.L2[4,1][1] ≈ sparse([1, 2, 5, 6, 3, 4, 5, 6], [1, 1, 1, 1, 2, 2, 2, 2], [0.0, 0.0, 2.532843602293451, 0.0, 0.0, 0.0, -2.532843602293451, 0.0], 6, 2)
    @test typeof(out2) ==  Muscade.AssemblyDirectLine
    @test out2.ming ≈ Inf
    @test out2.minλ ≈ Inf
    @test out2.Σλg  ≈ 0.
    @test out2.npos ≈ 0
end


nstep            = 6
ND               = 2
NA               = 1
Lv,Lvv,bigasm    = Muscade.preparebig(ND,NA,nstep,out1)

using GLMakie
fig      = Figure(size = (500,500))
display(fig) # open interactive window (gets closed down by "save")
axe      = Axis(fig[1,1],title="bigsparse Lvv sparsity",xlabel="i",ylabel="j")
(i,j,v)  = findnz(Lvv)
scatter!(axe,i,-j)

# TODO repair tests AFTER fixing FDsparsity.  Test both 1sd and 2rd order solvers.
@testset "preparebig Lvv" begin
#    @test size(Lv)         == (42,)
#    @test size(Lvv)        == (42,42)
#    @test Lvv.colptr       == [1,20,39,58,77,88,99,123,147,171,195,209,223,252,281,310,339,356,373,402,431,460,489,506,523,547,571,595,619,633,647,666,685,704,723,734,745,765,785,805,825,851,877]
#    @test Lvv.rowval[1:60] == [1,2,3,4,5,7,8,9,10,11,13,14,15,16,17,37,38,41,42,1,2,3,4,6,7,8,9,10,12,13,14,15,16,18,39,40,41,42,1,2,3,4,5,7,8,9,10,11,13,14,15,16,17,37,38,41,42,1,2,3]
end

@testset "preparebig ,bigasm" begin
#    @test bigasm.pgc      == [1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,43]
#    @test bigasm.pigr[1]' == [1 3 5 6 8 10 11 13 15 0 0 0 0 0 0 0 0 0 16; 20 22 24 25 27 29 30 32 34 0 0 0 0 0 0 0 0 0 35]
end

# state                 = [Muscade.State{1,ND,ND,@NamedTuple{γ::Float64}}(copy(state0)) for i = 1:nstep]
# γ = 9.
# Δt = 1.

# Muscade.assemblebig!(Lvv,Lv,bigasm,asm1,model1,dis,out1,state,nstep,Δt,NA,γ,(caller=:TestDirectXUA,))


stateXUA         = solve(DirectXUA{NA,ND};initialstate=state0,time=0:1.:5)

#end 

;