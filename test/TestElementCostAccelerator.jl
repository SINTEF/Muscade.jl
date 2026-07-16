#module TestElementCostAccelerator

using Test,StaticArrays
using Muscade



model           = Model()
n1              = addnode!(model,[0.,0.,0.]) 
n2              = addnode!(model,[1.,0.,0.])
n3              = addnode!(model,[2.,0.,0.])

req             = @request(gp(κ))
@functor with() cost(κ) = sum(κ.*κ)
@functor with() gap( κ) = κ[1]
material        = BeamCrossSection(EA=1. ,EI₂=2.,EI₃=3.,GJ=4.,μ=5.,ι₁=6.,   Ca₁=7., Ca₂=8.,Ca₃=9.)
elementkwargs = (mat=material,orient2=SVector(0.,1.,0.))
el1             = addelement!(model,ElementCost      ,[n1,n2]; req, ElementType=EulerBeam3D, elementkwargs, cost                                         )
el3             = addelement!(model,ElementConstraint,[n2,n3]; req, ElementType=EulerBeam3D, elementkwargs, gap ,λinod=1, λfield=:λκ₁, mode=Muscade.equal)

# TODO test the addin! accelerator with constraint on cost on element, or constraint on strain on element.  I don't think this works...

initialstate    = initialize!(model;time=0.)

#end
