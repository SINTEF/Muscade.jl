module TestStaticX

using Test,StaticArrays,SparseArrays
using Muscade

include("SomeElements.jl")


model           = Model(:TestModel)
n1              = addnode!(model,π£[0,0,+100]) # turbine
n2              = addnode!(model,π£[])  # Anod for turbine 
n3              = addnode!(model,π£[])  # Anod for anchor
@once sea sea(t,x)        = SVector(1.,0.)*t
@once sky sky(t,x)        = SVector(0.,10.)
Ξ±(i)            = SVector(cos(i*2Ο/3),sin(i*2Ο/3))
e1              =  addelement!(model,Turbine   ,[n1,n2], seadrag=1e6, sea=sea, skydrag=1e5, sky=sky)
e2              = [addelement!(model,AnchorLine,[n1,n3], Ξxβtop=vcat(5*Ξ±(i),[0.]), xβbot=250*Ξ±(i), L=290., buoyancy=-5e3) for iβ0:2]
state           = solve(staticX;model,time=[0.,1.],verbose=false)
step = 1
@testset "StaticX" begin
    @test  state[step].Ξ β [0.0, 0.0, 0.0]
    @test  state[step].X[1] β  [-5.332268523655259, 21.09778288272267, 0.011304253608808651]
    @test  state[step].U[1] β  Float64[]
    @test  state[step].A β [0.0, 0.0, 0.0, 0.0]
    @test  state[step].time β 0.
end

dis         = Muscade.Disassembler(model)
dofgr       = Muscade.allXdofs(model,dis)
s           = deepcopy(state[step])
Muscade.decrement!(s,0,[1.,1.,-1.],dofgr)
#s[dofgr]    = [1.,1.,1.]
# @testset "AllXdofs construction" begin
# #    @test  dofgr.scale β [1., 1., 1.]
# #    @test  state[step][dofgr] β [-5.332268523655259, 21.09778288272267, 0.011304253608808651]
# #    @test  s[dofgr] β [-6.332268523655259, 20.09778288272267, 1.011304253608808651]
# end

#using GLMakie
#fig      = Figure(resolution = (2000,1500))
#display(fig) # open interactive window (gets closed down by "save")
#axe      = Axis3(fig[1,1],title="Muscade made this drawing",xlabel="X",ylabel="Y",zlabel="Z",aspect=:data,viewmode=:fit,perspectiveness=.5)
#draw(axe,state[step])
#save("C:\\Users\\philippem\\C home\\GIT\\Muscade.jl\\test\\first_light.jpg",fig)

include("GLMakieTester.jl")
axe = SpyAxe()
draw(axe,state[step],ieletyp=[1,2])
@testset "drawing" begin
    @test  axe.call[1].fun == :lines!
    @test  axe.call[1].args[1] β [-5.332268523655259 -5.332268523655259; 21.09778288272267 21.09778288272267; 90.0 110.0]
    @test  axe.call[1].kwargs[:color] == :orange
    @test  axe.call[1].kwargs[:linewidth] == 5
    @test  axe.call[2].fun == :lines!
    @test  axe.call[2].args[1][1:5]β[144.05988106384137, 129.6206341588945, 115.1813872539476, 100.74214034900072, 86.30289344405384]
    @test  axe.call[2].kwargs[:color] == :blue
    @test  axe.call[2].kwargs[:linewidth] == 2
end

out1,dofid1 = getdof(state[1],field=:tx1)
out2,dofid2 = getdof(state   ,field=:tx1)
out3,dofid3 = getdof(state[1],class=:A,field=:ΞL)
out4,dofid4 = getdof(state[1],field=:tx1,nodID=[n1])
out5,dofid5 = getdof(state   ,field=:tx1,nodID=[n1])
@testset "getdof" begin
    @test  out1 β [-5.332268523655259;;]
    @test  dofid1 == DofID[DofID(:X, 1)]
    @test  out2 β [-5.332268523655259;;;20.184170880401076]
    @test  dofid2 == DofID[DofID(:X, 1)]
    @test  out3 β [0.0;;]
    @test  dofid3 == DofID[DofID(:A, 3)]
    @test  out4 β [-5.332268523655259;;]
    @test  dofid4 == DofID[DofID(:X, 1)]
    @test  out5 β [-5.332268523655259;;; 20.18417088040054]
    @test  dofid5 == DofID[DofID(:X, 1)]
end
req     = @request cr,ltf
out,key = getresult(state,req, eleID=e2)
cr      = out[key.cr ,:,1] # out[ikey,iele,istep]
ltf     = out[key.ltf,:,1]
@testset "getdof" begin
    @test  cr  β [118.69592125130082, 23.961941539907585, 235.51441727552435]
    @test  ltf β [183.68229160771097, 121.62396272109176, 238.96209627282917]
end

out,key = getresult(state[1],req, eleID=e2)
cr      = out[key.cr ,:] # out[ikey,iele]
ltf     = out[key.ltf,:]
@testset "getdof" begin
    @test  cr  β [118.69592125130082, 23.961941539907585, 235.51441727552435]
    @test  ltf β [183.68229160771097, 121.62396272109176, 238.96209627282917]
end

end
