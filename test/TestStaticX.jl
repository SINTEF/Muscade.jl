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
e1              =  addelement!(model,Turbine   ,[n1,n2], seadrag=1e6, sea=sea, skydrag=1e5, sky=sky)
e2              = [addelement!(model,AnchorLine,[n1,n3], Δxₘtop=vcat(5*α(i),[0.]), xₘbot=250*α(i), L=290., buoyancy=-5e3) for i∈0:2]
state           = solve(StaticX;model,time=[0.],verbose=false)

@testset "StaticX" begin
    @test  state[1].Λ ≈ [0.0, 0.0, 0.0]
    @test  state[1].X[1] ≈  [-17.46832446885514, -24.570658899684172, 0.011313890183180228]
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
    @test  state[1][dofgr] ≈ [-17.46832446885514, -24.570658899684172, 0.011313890183180228]
    @test  s[dofgr] ≈ [1.,1.,1.]
end

#using GLMakie
#fig      = Figure(resolution = (2000,1500))
#display(fig) # open interactive window (gets closed down by "save")
#axe      = Axis3(fig[1,1],title="Muscade made this drawing",xlabel="X",ylabel="Y",zlabel="Z",aspect=:data,viewmode=:fit,perspectiveness=.5)
#draw(axe,state[1])
#save("C:\\Users\\philippem\\C home\\GIT\\Muscade.jl\\test\\first_light.jpg",fig)

include("GLMakieTester.jl")
axe = SpyAxe()
draw(axe,state[1])
@testset "drawing" begin
    @test  axe.data[1].fun == :lines!
    @test  axe.data[1].args[1] ≈ [-17.46832446885514 -17.46832446885514; -24.570658899684172 -24.570658899684172; 90.0 110.0]
    @test  axe.data[1].kwargs[:color] == :orange
    @test  axe.data[1].kwargs[:linewidth] == 5
    @test  axe.data[2].fun == :lines!
    @test  axe.data[2].args[1][1:5]≈[ 220.71174584912032,197.39370681663635,174.0756677841524,150.75762875166845,127.43958971918447]
    @test  axe.data[2].kwargs[:color] == :blue
    @test  axe.data[2].kwargs[:linewidth] == 2
end
end
