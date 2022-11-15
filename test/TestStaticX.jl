module TestStaticX

using Test,StaticArrays,SparseArrays
using Muscade
using Muscade.ElTest

include("SomeElements.jl")


model           = Model(:TestModel)
n1              = addnode!(model,𝕣[0,0,+100]) # turbine
n2              = addnode!(model,𝕣[])  # Anod for turbine 
n3              = addnode!(model,𝕣[])  # Anod for anchor
sea(t,x)        = SVector(1.,0.)
sky(t,x)        = SVector(0.,10.)
α(i)            = SVector(cos(i*2π/3),sin(i*2π/3))
e1              = addelement!(model,Turbine   ,[n1,n2], seadrag=1e6, sea=sea, skydrag=1e5, sky=sky)
e2              = [addelement!(model,AnchorLine,[n1,n3], Δxₘtop=vcat(5*α(i),[0.]), xₘbot=150*α(i), L=180., buoyancy=-5e3) for i∈0:2]
state           = step!(StaticX;model,time = [0.],verbose=false)
@testset "StaticX" begin
    @test  state[1].Λ ≈ [0.0, 0.0, 0.0]
    @test  state[1].X[1] ≈  [-9.696622088552829, -12.438188104101995, 0.005763499207626334]
    @test  state[1].U[1] ≈  Float64[]
    @test  state[1].A ≈ [0.0, 0.0, 0.0, 0.0]
    @test  state[1].t ≈ 0.
end

dis         = Muscade.Disassembler(model)
dofgr       = Muscade.AllXdofs(model,dis)
s           = deepcopy(state[1])
s[dofgr]    = [1.,1.,1.]
@testset "AllXdofs construction" begin
    @test  dofgr.scale ≈ [1.0, 1.0, 1.0]
    @test  state[1][dofgr] ≈ [-9.696622088552829,-12.438188104101995,0.005763499207626334]
    @test  s[dofgr] ≈ [1.,1.,1.]
end
end
