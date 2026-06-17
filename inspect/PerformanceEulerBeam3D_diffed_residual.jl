module Performance

using Profile,ProfileView
using BenchmarkTools
using Muscade
using StaticArrays
using Muscade.Toolbox


model           = Model(:TestModel)
node1           = addnode!(model,𝕣[0,0,0])
node2           = addnode!(model,𝕣[4,3,0])
elnod           = [model.nod[n.inod] for n∈[node1,node2]]
mat             = BeamCrossSection(EA=10.,EI₂=3.,EI₃=3.,GJ=4.,μ=1.0,ι₁=1.)

const beam      = EulerBeam3D(elnod;mat,orient2=SVector(0.,1.,0.))

const t         = 0.
const SP        = (;)
const dbg       = (status=:testing,)
const displacement =  SVector(0.,0.,0.,0.,0.,0.,  0.,0.,0.,0.,0.,0.); 
const velocity     =  SVector(0.,0.,0.,0.,0.,0.,  0.,0.,0.,0.,0.,0.); 
const acceleration =  SVector(0.,0.,0.,0.,0.,0.,  0.,0.,0.,0.,0.,0.); 
const X = (displacement,)#,velocity,acceleration)
const U            = (SVector{0,𝕣}(),)
const A            =  SVector{0,𝕣}()

mission = :time
if mission == :eval
    out = Muscade.diffed_residual(beam; X,U,A,t,SP)
elseif mission == :time
    # 121.500 μs baseline 
    # 104.400 μs if ...∂X₀ are computed with kinematics{:direct}
    # 101.800 μs with improvements in Adiff.jl
    #  87.300 μs with improvements in Rotation.jl
    out =  Muscade.diffed_residual(beam; X,U,A,t,SP)
    @btime Muscade.diffed_residual(beam; X,U,A,t,SP)
elseif mission == :profile
    Muscade.diffed_residual(beam; X,U,A,t,SP)
    Profile.clear()
    Profile.@profile for i=1:250000
        local out = Muscade.diffed_residual(beam; X,U,A,t,SP)
    end
    ProfileView.view(fontsize=30);
    # After clicking on a bar in the flame diagram, you can type warntype_last() and see the result of 
    # code_warntype for the call represented by that bar.
end


end # module


