#module TestDirectXUA

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




model           = Model(:TrueModel)
n1              = addnode!(model,𝕣[0])  
n2              = addnode!(model,𝕣[1])  
n3              = addnode!(model,𝕣[]) # anode for spring
e1              = addelement!(model,El1,[n1], K=1.,C=0.05,M=1.)
e2              = addelement!(model,El1,[n2], K=0.,C=0.0 ,M=1.)
e3              = addelement!(model,Spring{1},[n1,n2,n3], EI=1.1)
@once f(x) = x^2
@once l1(tx1,t) = (tx1-0.1*sin(t))^2
@once l2(tx1,t) = (tx1-0.1*cos(t))^2
e4              = addelement!(model,SingleDofCost,[n3];class=:A,field=:ΞL₀,     cost=f)
e5              = addelement!(model,SingleDofCost,[n3];class=:A,field=:ΞEI,     cost=f)
e6              = addelement!(model,SingleDofCost,[n1];class=:A,field=:ΞC ,     cost=f)
e7              = addelement!(model,SingleDofCost,[n1];class=:A,field=:ΞM ,     cost=f)
e8              = addelement!(model,SingleDofCost,[n2];class=:A,field=:ΞC ,     cost=f)
e9              = addelement!(model,SingleDofCost,[n2];class=:A,field=:ΞM ,     cost=f)
e10             = addelement!(model,SingleUdof   ,[n1];Xfield=:tx1,Ufield=:utx1,cost=f)
e11             = addelement!(model,SingleUdof   ,[n2];Xfield=:tx1,Ufield=:utx1,cost=f)
e12             = addelement!(model,SingleDofCost,[n1];class=:X,field=:tx1,     cost=l1)
e13             = addelement!(model,SingleDofCost,[n2];class=:X,field=:tx1,     cost=l2)

state0          = initialize!(model)   

nstep            = 6
OX               = 2
OU               = 0
IA               = 1
FAST             = true

Δt = 1.

dis             = state0.dis
out,asm,dofgr = Muscade.prepare(Muscade.AssemblyDirect{OX,OU,IA},model,dis,fastresidual=true)#;Uwhite=true,Xwhite=true,XUindep=true,UAindep=true,XAindep=true)
zero!(out)
state           = [Muscade.State{1,OX+1,OU+1,@NamedTuple{γ::Float64,iter::Int64}}(copy(state0)) for i = 1:nstep]
for i=1:nstep
    state[i].time = Δt*i
end

Muscade.assemble!(out,asm,dis,model,state[1],(;))

pattern    = Muscade.makepattern(OX,OU,IA,nstep,out)
# using Spy,GLMakie
# fig = spypattern(pattern)
# save("C:\\Users\\philippem\\.julia\\dev\\Muscade\\spypattern.jpg",fig)

Lvv,Lv,Lvvasm,Lvasm,Lvdis  = Muscade.preparebig(OX,OU,IA,nstep,out)
Muscade.assemblebig!(Lvv,Lv,Lvvasm,Lvasm,asm,model,dis,out,state,nstep,Δt,(γ=0.,iter=1),(caller=:TestDirectXUA,))

# using Spy,GLMakie
# fig = Spy.spy(Lvv,title="bigsparse Lvv sparsity",size=500)
# save("C:\\Users\\philippem\\.julia\\dev\\Muscade\\spy.jpg",fig)

stateXUA         = solve(DirectXUA{OX,OU,IA};initialstate=state0,time=0:1.:5,maxiter=10,verbose=false)

@testset "prepare_out" begin
    @test out.L1[1] ≈ [[0.0, 0.0]]
    @test out.L1[2] ≈ [[0.055883099639785175, -0.1920340573300732], [0.0, 0.0], [0.0, 0.0]]
    @test out.L1[3] ≈ [[0.0, 0.0, 0.0, 0.0]]
    @test out.L1[4] ≈ [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
    @test size(out.L2[1,1]) == (0,0)
    @test out.L2[2,2][1,1] ≈ sparse([1, 2, 1, 2], [1, 1, 2, 2], [2.0, 0.0, 0.0, 2.0], 2, 2)  
    @test out.L2[2,2][1,2] ≈ sparse([1, 2, 1, 2], [1, 1, 2, 2], [0.0, 0.0, 0.0, 0.0], 2, 2)  
    @test out.L2[2,2][1,3] ≈ sparse([1, 2, 1, 2], [1, 1, 2, 2], [0.0, 0.0, 0.0, 0.0], 2, 2)
    @test out.L2[2,1][1,1] ≈ sparse([1, 2, 1, 2], [1, 1, 2, 2], [2.1,-1.1,-1.1, 1.1], 2, 2)
    @test out.L2[2,1][2,1] ≈ sparse([1, 2, 1, 2], [1, 1, 2, 2], [.05, 0.0, 0.0, 0.0], 2, 2)
    @test out.L2[2,1][3,1] ≈ sparse([1, 2, 1, 2], [1, 1, 2, 2], [1.0, 0.0, 0.0, 1.0], 2, 2)
    @test out.L2[3,3][1,1] ≈ sparse([1, 2, 3, 4], [1, 2, 3, 4], [0.0, 0.0, 2.0, 2.0], 4, 4)
    @test out.L2[3,4][1,1] ≈ sparse([1, 1, 2, 2], [1, 2, 3, 4], [0.0, 0.0, 0.0, 0.0], 4, 6)
    @test out.L2[4,4][1,1] ≈ sparse([1, 2, 1, 2, 3, 4, 3, 4, 5, 6, 5, 6], [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6], [2.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 2.0], 6, 6)
end
@testset "prepare_asm" begin
    @test asm[1,1]  ≈ [1 2]      # asm[iarray,ieletyp][ieledof,iele] -> idof|inz
    @test asm[1,2]  ≈ [1; 2;;]
    @test asm[2,1]  ≈ [1 2] 
    @test asm[2,2]  ≈ [1; 2;;]
    @test asm[3,1]  ≈ [1 2] 
    @test asm[3,2]  ≈ Matrix{Int64}(undef, 0, 1)
    @test asm[4,1]  ≈ [1 3; 2 4]
    @test asm[4,2]  ≈ [5; 6;;]
    @test asm[5,1]  ≈ [1 4]                   
    @test asm[20,1] ≈ [1 5; 2 6; 3 7; 4 8]    
end

@testset "preparebig Lvv" begin
    @test size(Lv)         == (54,)
    @test size(Lvv)        == (54,54)
    @test Lvv.colptr[1:50] == [1,11,21,45,69,75,81,85,89,101,113,137,161,168,175,180,185,199,213,241,269,277,285,291,297,311,325,353,381,389,397,403,409,421,433,457,481,488,495,500,505,515,525,549,573,579,585,589,593,613]
    @test Lvv.rowval[1:60] == [3,4,5,7,11,12,49,50,53,54,3,4,6,8,11,12,51,52,53,54,1,2,3,4,5,7,9,10,11,12,13,15,17,18,19,20,21,23,27,28,49,50,53,54,1,2,3,4,6,8,9,10,11,12,14,16,17,18,19,20]
end

@testset "preparebig ,Lvvasm" begin
    @test Lvasm      == [1, 3, 5, 9, 11, 13, 17, 19, 21, 25, 27, 29, 33, 35, 37, 41, 43, 45, 49, 55]
    @test Lvdis      == [1, 3, 5, 9, 11, 13, 17, 19, 21, 25, 27, 29, 33, 35, 37, 41, 43, 45, 49, 55]
    @test Lvvasm.colptr' == [1  5  16  21  26  37  43  49  62  69  75  88  95  100  111  117  121  132  137  156]
    @test Lvvasm.nzval[1]' == [1  2  11  12]
    @test Lvvasm.nzval[4]' == [7  8  9  10  17  18  19  20]
    @test Lvvasm.nzval[155]' == [611  612  631  632  651  652  671  672  697  698  723  724]
end

#end 

