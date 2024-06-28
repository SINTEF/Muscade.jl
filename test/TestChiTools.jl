module TestChiTools
using Test,StaticArrays,SparseArrays
using Muscade
using Muscade: DofID,EleID,NodID

include("SomeElements.jl")


c0  = 3.
c1  = ∂ℝ{1,1,𝕣}(3.,SVector{1}(1.))
c2  = ∂ℝ{2,1,𝕣}(3.,SVector{1}(2.))
c12 = ∂ℝ{2, 1, ∂ℝ{1, 1, Float64}}(∂ℝ{1, 1, Float64}(3.0, [1.0]), ∂ℝ{1, 1, Float64}[∂ℝ{1, 1, Float64}(2.0, [12.0])])
T0  = typeof(c0)
T1  = typeof(c1)
T2  = typeof(c2)
T12 = typeof(c12)

@testset "χcast adiffs" begin
    @test Muscade.χcast(T0,c0) === c0
    @test Muscade.χcast(T1,c0) === ∂ℝ{1, 1, Float64}(3.0, [0.0])
    @test Muscade.χcast(T2,c0) === ∂ℝ{2, 1, Float64}(3.0, [0.0])
    @test Muscade.χcast(T12,c0) === ∂ℝ{2, 1, ∂ℝ{1, 1, Float64}}(∂ℝ{1, 1, Float64}(3.0, [0.0]), ∂ℝ{1, 1, Float64}[∂ℝ{1, 1, Float64}(0.0, [0.0])])
    @test Muscade.χcast(T0,c1) === c0
    @test Muscade.χcast(T1,c1) === c1
    @test Muscade.χcast(T2,c1) === ∂ℝ{2, 1, Float64}(3.0, [0.0])
    @test Muscade.χcast(T12,c1) === ∂ℝ{2, 1, ∂ℝ{1, 1, Float64}}(∂ℝ{1, 1, Float64}(3.0, [1.0]), ∂ℝ{1, 1, Float64}[∂ℝ{1, 1, Float64}(0.0, [0.0])])
    @test Muscade.χcast(T0,c2) === c0
    @test Muscade.χcast(T1,c2) === ∂ℝ{1, 1, Float64}(3.0, [0.0])
    @test Muscade.χcast(T2,c2) === c2
    @test Muscade.χcast(T12,c2) ===  ∂ℝ{2, 1, ∂ℝ{1, 1, Float64}}(∂ℝ{1, 1, Float64}(3.0, [0.0]), ∂ℝ{1, 1, Float64}[∂ℝ{1, 1, Float64}(2.0, [0.0])])
    @test Muscade.χcast(T0,c12) === c0
    @test Muscade.χcast(T1,c12) === c1
    @test Muscade.χcast(T2,c12) === c2  
    @test Muscade.χcast(T12,c12) === c12
end


@testset "χcasttype adiffs" begin
    @test Muscade.χcasttype(T0,T0) === T0
    @test Muscade.χcasttype(T1,T0) === T1
    @test Muscade.χcasttype(T2,T0) === T2
    @test Muscade.χcasttype(T12,T0) === T12
    @test Muscade.χcasttype(T0,T1) === T0
    @test Muscade.χcasttype(T1,T1) === T1
    @test Muscade.χcasttype(T2,T1) === T2
    @test Muscade.χcasttype(T12,T1) === T12
    @test Muscade.χcasttype(T0,T2) === T0
    @test Muscade.χcasttype(T1,T2) === T1
    @test Muscade.χcasttype(T2,T2) === T2
    @test Muscade.χcasttype(T12,T2) ===  T12
    @test Muscade.χcasttype(T0,T12) === T0
    @test Muscade.χcasttype(T1,T12) === T1
    @test Muscade.χcasttype(T2,T12) === T2  
    @test Muscade.χcasttype(T12,T12) === T12
end


T∂  = ∂ℝ{1,1,𝕣}
y = (3,290.,:a)
x = (a=y,b=:helloworld)
X = [x,x] 

χ2 = Muscade.χcast(T∂,X)

@testset "χcast χ-structures" begin
    @test χ2[1].a[1] ≈ 3
    @test χ2[1].a[2] ≈ variate{1}(290.)
    @test χ2[1].b   ==  :helloworld
    @test χ2[1]     ==  χ2[2]
end

@testset "typecase χ-structures" begin
    @test Muscade.χcasttype(T∂,typeof(y)) == Tuple{Int64,T∂, Symbol}
    @test Muscade.χcasttype(T∂,typeof(x)) == @NamedTuple{a::Tuple{Int64, T∂, Symbol}, b::Symbol}
    @test Muscade.χcasttype(T∂,SVector{2,𝕣}) == SVector{2,T∂}
    @test Muscade.χcasttype(T∂,SMatrix{2,2,𝕣}) == SMatrix{2,2,T∂}
    @test Muscade.χcasttype(T∂,SArray{Tuple{2,2,3},𝕣}) == SArray{Tuple{2,2,3},T∂}
end
χmem = Muscade.χalloc_(∂ℝ{1,2,𝕣},X)
@testset "χalloc" begin
    @test eltype(χmem) == @NamedTuple{a::Tuple{Int64, ∂ℝ{1, 2, Float64},Symbol}, b::Symbol}
    @test length(χmem) == 2
end

###

model           = Model(:TestModel)
n1              = addnode!(model,𝕣[0,0,+100]) # turbine
n2              = addnode!(model,𝕣[])  # Anod for turbine 
n3              = addnode!(model,𝕣[])  # Anod for anchor
@once sea(t,x)  = SVector(1.,0.)*t
@once sky(t,x)  = SVector(0.,10.)
α(i)            = SVector(cos(i*2π/3),sin(i*2π/3))
e1              =  addelement!(model,Turbine   ,[n1,n2], seadrag=1e6, sea=sea, skydrag=1e5, sky=sky)
e2              = [addelement!(model,AnchorLine,[n1,n3], Δxₘtop=vcat(5*α(i),[0.]), xₘbot=250*α(i), L=290., buoyancy=-5e3) for i∈0:2]
initialstate    = initialize!(model)
state           = solve(StaticX;initialstate,time=[0.,1.],verbose=false)
@testset "initalstate" begin
    @test initialstate.χ[1]        == [nothing]
    @test initialstate.χ[2][1].a[1] ≈ 3
    @test initialstate.χ[2][1].a[2] ≈ 290.
    @test initialstate.χ[2][1].b   ==  :helloworld
    @test initialstate.χ[2][1]     ==  initialstate.χ[2][2]
end

end