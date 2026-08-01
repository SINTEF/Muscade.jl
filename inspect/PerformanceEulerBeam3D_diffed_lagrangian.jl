#module Performance

using SnoopCompileCore, SnoopCompile#, Cthulu
using Profile,ProfileView
using BenchmarkTools
using Muscade
using Muscade.Toolbox
using StaticArrays

# module TmpModule
#     export BeamCrossSection, EulerBeam3D,EulerBeam3DwithStrainGauge
#     include("../toolbox/BeamElement.jl")
#     include("../toolbox/StrainGaugeOnBeamElement.jl")
# end


model            = Model(:TestModel)
node1            = addnode!(model,𝕣[0,0,0])
node2            = addnode!(model,𝕣[4,0,0])
elnod            = [model.nod[n.inod] for n∈[node1,node2]]
mat              = BeamCrossSection(EA=10.,EI₂=3.,EI₃=3.,GJ=4.,μ=1.,ι₁=1.0)
P                = SMatrix{3,5}(0.,.5,0.,  0.,0,.5,   0.,-.5,0.,  0.,0,-.5,  0.,.5,0.   )
D                = SMatrix{3,5}(1.,0.,0.,  1.,0.,0.,  1.,0.,0.,   1.,0.,0.,  1/√2,0,1/√2)
L                = 0.1

const beam      = EulerBeam3D(elnod;mat,orient2=SVector(0.,1.,0.))

const t         = 0.
const SP        = (;)
const dbg       = (status=:testing,)
const displacement =  SVector(0,0,0, 0,.1,0, 0,0,0, 0,-.1,0); 
const velocity     =  SVector(0.,0.,0.,0.,0.,0.,  0.,0.,0.,0.,0.,0.); 
const acceleration =  SVector(0.,0.,0.,0.,0.,0.,  0.,0.,0.,0.,0.,0.); 
Λ   =  SVector(0,0,0, 0,.1,0, 0,0,0, 0,-.1,0)
X   = (displacement,)
U   = (SVector{0,𝕣}(),)
A   = SVector{0,𝕣}()

# function straincost(eleres,X,U,A,t) 
straincost = Functor{:straincost}(;)  
# function (o::Functor{:straincost})(eleres,X,U,A,t) 
#     σ  = 15e-6
#     ε  = eleres.ε
#     εₘ = SVector(cos(t),0.,-cos(t),0.,cos(t)/2)*0.001  
#     Δε = ε-εₘ
#     cost = (Δε⋅Δε)/(2σ^2)
#     return cost
# end

#invs = @snoop_invalidations ElementCostAndConstraint(elnod;
const costedbeam =  ElementCostAndConstraint(elnod;
                            req = @request(ε),
                            cost=straincost,
                            TargetElement=EulerBeam3DwithStrainGauge,
                            elementkwargs = (P,D,
                                              elementkwargs=(mat=mat,orient2=SVector(0.,1.,0.))))

display("tic")
out = Muscade.diffed_lagrangian(costedbeam;Λ,X,U,A,t=0.)
display("tac")
out = Muscade.diffed_lagrangian(costedbeam;Λ,X,U,A,t=0.)
display("toc")
# invs = @snoop_invalidations Muscade.diffed_lagrangian(costedbeam;Λ,X,U,A,t=0.)

# 2025-07-22
# Julia 1.10: ttfx=2:05s @time=522μs secondrun=20s  
# Julia 1.11: ttfx=1:52s @time=548μs secondrun=20s  
# Julia 1.12: ttfx=1;33s @time=610μs secondrun=20s second@time=200μs

# mission = :time
# if mission == :report
#     out = Muscade.diffed_lagrangian(costedbeam;Λ,X,U,A,t=0.)
# elseif mission == :time
#     out = Muscade.diffed_lagrangian(costedbeam;Λ,X,U,A,t=0.)
#     @btime out = Muscade.diffed_lagrangian(costedbeam;Λ,X,U,A,t=0.)
# elseif mission == :profile
#     out = Muscade.diffed_lagrangian(costedbeam;Λ,X,U,A,t=0.)
#     Profile.clear()
#     Profile.@profile for i=1:100000
#         local out = Muscade.diffed_lagrangian(costedbeam;Λ,X,U,A,t=0.)
#     end
#     ProfileView.view(fontsize=30);
#     # After clicking on a bar in the flame diagram, you can type warntype_last() and see the result of 
#     # code_warntype for the call represented by that bar.
# end

;
#end # module


