using Test, Muscade, StaticArrays, LinearAlgebra
using SnoopCompileCore, SnoopCompile, AbstractTrees, ProfileView
using Profile
using BenchmarkTools

include("../toolbox/BeamElement.jl")
include("../examples/StrainGaugeOnBeamElement.jl")

model        = Model(:TestModel)
node1        = addnode!(model,𝕣[0,0,0])
node2        = addnode!(model,𝕣[4,3,0])
node3        = addnode!(model,𝕣[]) # Unod
elnod        = [model.nod[n.inod] for n∈[node1,node2,node3]]
mat          = BeamCrossSection(EA=10.,EI₂=3.,EI₃=3.,GJ=4.,μ=1.,ι₁=1.0)
beam         = EulerBeam3D{false}(elnod;mat,orient2=SVector(0.,1.,0.))
strainbeam = StrainGaugeOnEulerBeam3D(elnod;
                                         P             = SMatrix{3,4}(0.,0.,.05, 0.,0.05,0.,  0.,0.,-.05,  0.,-.05,0.),
                                         D             = SMatrix{3,4}(1.,0.,0.,  1.,0.,0.,    1.,0.,0.,    1.,0.,0.  ),
                                         TargetElement   = EulerBeam3D{true},
                                         elementkwargs = (mat     = mat,
                                                          orient2 = SVector(0.,0.,1.)))
@functor with(σε=1.) cost(eleres,X,U,A,t) = sum((eleres.ε/σε).^2)/2
coststrainbeam = ElementCostAndConstraint(elnod;
                        req           = @request(ε),
                        cost          = cost,
                        TargetElement   = StrainGaugeOnEulerBeam3D,
                        elementkwargs = (P             = SMatrix{3,4}(0.,0.,.05, 0.,0.05,0.,  0.,0.,-.05,  0.,-.05,0.),
                                         D             = SMatrix{3,4}(1.,0.,0.,  1.,0.,0.,    1.,0.,0.,    1.,0.,0.  ),
                                         TargetElement   = EulerBeam3D{true},
                                         elementkwargs = (mat     = mat,
                                                          orient2 = SVector(0.,0.,1.))))


displacement =  SVector(0.,0.,0.,0.,0.,0.,  0.,0.,0.,0.,0.,0.); 
velocity     =  SVector(0.,0.,0.,0.,0.,0.,  0.,0.,0.,0.,0.,0.); 
acceleration =  SVector(0.,0.,0.,0.,0.,0.,  0.,0.,0.,0.,0.,0.); 
Λ            = displacement
#X            = (displacement,)
#X            = (displacement,velocity)
X            = (displacement,velocity,acceleration)
U            = (SVector{3,𝕣}(1,2,3),)
A            = SVector{0,𝕣}()
t,SP         = 0.,(;)

# Study https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-captured 

## Invalidation

# invalidation = @snoop_invalidation diffed_lagrangian(beam; Λ,X,U,A,t,SP) # 

## Inference

# display("First compilation")
# #out            = diffed_lagrangian(coststrainbeam; Λ,X,U,A,t,SP)
# #out            = diffed_residual(strainbeam; X,U,A,t,SP)
# out            = diffed_residual(beam; X,U=(SVector{0,𝕣}(),),A,t,SP)
# display(out.inftyp)
# display(out.rettyp)

# display("Snoop compilation")
# #inference     = @snoop_inference diffed_lagrangian(coststrainbeam; Λ,X,U,A,t,SP) # 
# #inference     = @snoop_inference diffed_residual(strainbeam; X,U,A,t,SP) # 
# inference     = @snoop_inference diffed_residual(beam; X,U=(SVector{0,𝕣}(),),A,t,SP) # 

# display("Snooping done")
# ProfileView.view(flamegraph(inference))
# trig          = inference_triggers(inference)
# flat          = flatten(inference,tmin = 0.0, sortby=exclusive)
# source        = accumulate_by_source(flat; tmin = 1., by=inclusive) # Accumulate by source (aka trigger), i.e. methods triggering the need for inference
# # MuscadeSource = filtermod(Muscade, source)
# # summary(MuscadeSource[1])
# # suggest(MuscadeSource[1])
# # Trigger trees

# ## Bleeding edge
# # Julia 1.12: @trace_compile solve(EigXU{OX,OU};Δω, p, nmod,initialstate,verbose=true,verbosity=1,tol=1e-20,σₓᵤ)


mission = :profile
if mission == :report
    out = Muscade.diffed_residual(strainbeam; X,U,A,t,SP)
elseif mission == :time
    out = Muscade.diffed_residual(strainbeam; X,U,A,t,SP)
    @btime out =diffed_residual(strainbeam; X,U,A,t,SP)
elseif mission == :profile
    out = Muscade.diffed_residual(strainbeam; X,U,A,t,SP)
    Profile.clear()
    Profile.@profile for i=1:30000
        local out = Muscade.diffed_residual(strainbeam; X,U,A,t,SP)
    end
    ProfileView.view(fontsize=30);
    # After clicking on a bar in the flame diagram, you can type warntype_last() and see the result of 
    # code_warntype for the call represented by that bar.
end

;