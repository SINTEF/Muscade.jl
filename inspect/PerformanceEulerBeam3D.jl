using Profile,ProfileView,BenchmarkTools
using Muscade
using StaticArrays

model           = Model(:TestModel)
node1           = addnode!(model,𝕣[0,0,0])
node2           = addnode!(model,𝕣[4,3,0])
elnod           = [model.nod[n.inod] for n∈[node1,node2]]
mat             = Muscade.BeamCrossSection(EA=10.,EI=3.,GJ=4.)

const beam      = EulerBeam3D(elnod;mat,orient2=SVector(0.,1.,0.))

const t= 0.
const χo = (nothing,nothing)
const χcv = identity
const SP = (;)
const dbg  = (status=:testing,)
const x = SVector(1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.)
const xv = variate{1,12}(x)
const X = (xv,)
const U = (SVector{0,𝕣}(),)
const A = SVector{0,𝕣}()

mission = :time
if mission == :report
    R,χ,FB=residual(beam,   X,U,A,t,χo,χcv,SP,dbg)
    r,K = value_∂{1,12}(R)
    display(r')
    display(K)
elseif mission == :time
    R,χ,FB=residual(beam,   X,U,A,t,χo,χcv,SP,dbg)
    @btime residual(beam,   X,U,A,t,χo,χcv,SP,dbg)
elseif mission == :profile
    Profile.clear()
    Profile.@profile for i=1:10000
        local R,χ,FB=residual(beam,   X,U,A,t,χo,χcv,SP,dbg)
    end
    ProfileView.view(fontsize=30);
    # After clicking on a bar in the flame diagram, you can type warntype_last() and see the result of 
    # code_warntype for the call represented by that bar.
end
;


