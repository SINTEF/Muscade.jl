#module TestStaticX

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

step = 1
@testset "StaticX" begin
    @test  state[step].Λ ≈ [0.0, 0.0, 0.0]
    @test  state[step].X[1] ≈  [20.184170880401076, 10.987078031136829, -0.016856115935358795]
    @test  state[step].U[1] ≈  Float64[]
    @test  state[step].A ≈ [0.0, 0.0, 0.0, 0.0]
    @test  state[step].t ≈ 0.
end

dis         = Muscade.Disassembler(model)
dofgr       = Muscade.AllXdofs(model,dis)
s           = deepcopy(state[step])
s[dofgr]    = [1.,1.,1.]
@testset "AllXdofs construction" begin
    @test  dofgr.scale ≈ [1., 1., 1.]
    @test  state[step][dofgr] ≈ [20.184170880401076, 10.987078031136829, -0.016856115935358795]
    @test  s[dofgr] ≈ [1.,1.,1.]
end

#using GLMakie
#fig      = Figure(resolution = (2000,1500))
#display(fig) # open interactive window (gets closed down by "save")
#axe      = Axis3(fig[1,1],title="Muscade made this drawing",xlabel="X",ylabel="Y",zlabel="Z",aspect=:data,viewmode=:fit,perspectiveness=.5)
#draw(axe,state[step])
#save("C:\\Users\\philippem\\C home\\GIT\\Muscade.jl\\test\\first_light.jpg",fig)

include("GLMakieTester.jl")
axe = SpyAxe()
draw(axe,state[step])
@testset "drawing" begin
    @test  axe.call[1].fun == :lines!
    @test  axe.call[1].args[1] ≈ [20.184170880401076 20.184170880401076; 10.987078031136829 10.987078031136829; 90.0 110.0]
    @test  axe.call[1].kwargs[:color] == :orange
    @test  axe.call[1].kwargs[:linewidth] == 5
    @test  axe.call[2].fun == :lines!
    @test  axe.call[2].args[1][1:5]≈[82.97477173657569, 77.19564062047895, 71.41650950438223, 65.6373783882855, 59.85824727218878]
    @test  axe.call[2].kwargs[:color] == :blue
    @test  axe.call[2].kwargs[:linewidth] == 2
end

out1,dofid1 = getdof(state[1],field=:tx1)
out2,dofid2 = getdof(state   ,field=:tx1)
out3,dofid3 = getdof(state[1],class=:A,field=:ΔL)
out4,dofid4 = getdof(state[1],field=:tx1,nodID=[n1])
out5,dofid5 = getdof(state   ,field=:tx1,nodID=[n1])
@testset "getdof" begin
    @test  out1 ≈ [20.184170880401076;;]
    @test  dofid1 == DofID[DofID(:X, 1)]
    @test  out2 ≈ [20.184170880401076;;;]
    @test  dofid2 == DofID[DofID(:X, 1)]
    @test  out3 ≈ [0.0;;]
    @test  dofid3 == DofID[DofID(:A, 3)]
    @test  out4 ≈ [20.184170880401076;;]
    @test  dofid4 == DofID[DofID(:X, 1)]
    @test  out5 ≈ [20.184170880401076;;;]
    @test  dofid5 == DofID[DofID(:X, 1)]
end
req     = @request cr,ltf
out,key = getresult(state,req, eleID=e2)
cr      = out[key.cr ,:,1] # out[ikey,iele,istep]
ltf     = out[key.ltf,:,1]
@testset "getdof" begin
    @test  cr  ≈ [25.37276764221335, 89.41914791181317, 322.4125480714144]
    @test  ltf ≈ [122.77847339188848,166.98451899012264,272.9148394907886]
end

out,key = getresult(state[1],req, eleID=e2)
cr      = out[key.cr ,:] # out[ikey,iele]
ltf     = out[key.ltf,:]
@testset "getdof" begin
    @test  cr  ≈ [25.37276764221335, 89.41914791181317, 322.4125480714144]
    @test  ltf ≈ [122.77847339188848,166.98451899012264,272.9148394907886]
end

#end
