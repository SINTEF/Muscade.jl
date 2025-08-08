# TODO
# - normalize X dofs to maximum(abs.(X)) = 1 with the exception of reaction forces?
# - model a circle to avoid "inertial singularities"
# - what does λ mean? Report it, as a function of ω and imod.
# - do we indeed have convergence of geneig?  Introduce an option to get geneig toself check.


#invs = @snoop_invalidations 
# using Muscade, Test, StaticArrays,SparseArrays;


q          = 6
nel        = 2^q
Δi         = 2^3-1 # 8-1=7 sensor stations
istrain    = Δi:Δi:(nel-Δi) # leaving a gap in the set of sensors
iacc       = Δi:2Δi:(nel-Δi)
iels       = 1:nel
inaked     = setdiff(iels,istrain)  # complement
xn         = 1.
λn         = 1e3
un         = 1e3

## beam in space
include("../examples/BeamElement.jl")
include("../examples/StrainGaugeOnBeamElement.jl")
include("../examples/PositionElement.jl")

L    = 1;    # Beam length [m]
q    = 0.0;  # Uniform lateral load [N/m]
EI₂  = 1;    # Bending stiffness [Nm²]
EI₃  = 10;   # Bending stiffness [Nm²]
EA   = 1e6;  # Axial stiffness [N]
GJ   = 1e6;  # Torsional stiffness [Nm²]
μ    = 1;
ι₁   = 1;
hasU = true
const σε   = 100e-6*100 # precision of strain measurements
σx   = 1e-1
σu   = 1e-0
const σa   = 1e5

α           = iels/nel*2π
XnodeCoord  = hcat(cos.(α),sin.(α),zeros(nel,1))
UnodeCoord  = zeros(nel,3)
mat         = BeamCrossSection(;EA,EI₂,EI₃,GJ,μ,ι₁)
model       = Model(:TestModel)
Xnod        = addnode!(model,XnodeCoord)
Unod        = addnode!(model,UnodeCoord)
mesh        = hcat(Xnod,Xnod[mod_onebased.(iels.+1,nel)],Unod)
nakedmesh   = mesh[inaked,:]
strainmesh  = mesh[istrain,:]
accmesh     = reshape(Xnod[iacc],(length(iacc),1))
addelement!(model,EulerBeam3D{hasU},nakedmesh;mat=mat,orient2=SVector(0.,0.,1.))
addelement!(model,ElementCost,strainmesh;
                        req           = @request(ε),
                        cost          = (eleres,X,U,A,t) -> sum((eleres.ε/σε).^2)/2,
                        ElementType   = StrainGaugeOnEulerBeam3D,
                        elementkwargs = (P             = SMatrix{3,4}(0.,0.,.05, 0.,0.05,0.,  0.,0.,-.05,  0.,-.05,0.),
                                         D             = SMatrix{3,4}(1.,0.,0.,  1.,0.,0.,    1.,0.,0.,    1.,0.,0.  ),
                                         ElementType   = EulerBeam3D{true},
                                         elementkwargs = (mat     = mat,
                                                          orient2 = SVector(0.,0.,1.))))
addelement!(model,ElementCost,accmesh;
                        req           = @request(a),
                        cost          = (eleres,X,U,A,t) -> sum((eleres.a/σa).^2)/2,
                        ElementType   = Position3D,
                        elementkwargs = (P             = SMatrix{3,3}(0.,0.,.1,  0.,0.,.1,  0.,0.,.1),
                                         D             = SMatrix{3,3}(1.,0.,0.,  0.,1.,0.,    0.,0.,1.) ))
# Ucost
addelement!( model, SingleDofCost, Muscade.columnmatrix(Unod         )    ,class=:U, field=:t1,cost=QuadraticFunction(0.,0.5(σu^-2)))
addelement!( model, SingleDofCost, Muscade.columnmatrix(Unod         )    ,class=:U, field=:t2,cost=QuadraticFunction(0.,0.5(σu^-2)))
addelement!( model, SingleDofCost, Muscade.columnmatrix(Unod         )    ,class=:U, field=:t3,cost=QuadraticFunction(0.,0.5(σu^-2)))
# disp meas
addelement!( model, SingleDofCost, Muscade.columnmatrix(Xnod[istrain])    ,class=:X, field=:t1,cost=QuadraticFunction(0.,0.5(σx^-2)))
addelement!( model, SingleDofCost, Muscade.columnmatrix(Xnod[istrain])    ,class=:X, field=:t2,cost=QuadraticFunction(0.,0.5(σx^-2)))
addelement!( model, SingleDofCost, Muscade.columnmatrix(Xnod[istrain])    ,class=:X, field=:t3,cost=QuadraticFunction(0.,0.5(σx^-2)))

initialstate      = initialize!(model)   
initialstate.time     = 0.
OX,OU                 = 2,0

## N
σₓᵤ     = (X=(t1 =xn,t2 =xn,t3 =xn,r1 =xn,r2 =xn,r3 =xn,
             λt1=λn,λt2=λn,λt3=λn,λr1=λn,λr2=λn,λr3=λn),
          U=(t1 =un,t2 =un,t3 =un,r1 =un,r2 =un,r3 =un))

# ## testing makeXUnorm
# dis             = Muscade.Disassembler(model)
# dofgr           = Muscade.allXdofs(model,dis)
# N               = [zeros(getndof(dofgr)),zeros(getndof(dofgr)),zeros(getndof(dofgr))]
# Muscade.makeXUnorm!(N,dofgr,σₓᵤ)
# @show N

## EigXU analysis
Δω                = 2^-6 
p                 = 11
nmod              = 5



# using SnoopCompileCore, SnoopCompile, AbstractTrees, ProfileView
# @trace_compile solve(EigXU{OX,OU};Δω, p, nmod,initialstate,verbose=true,verbosity=1,tol=1e-20,σₓᵤ)
# #inference = @snoop_inference solve(EigXU{OX,OU};Δω, p, nmod,initialstate,verbose=true,verbosity=1,tol=1e-20,σₓᵤ)
# # io=open("inference.txt","w")
# # print_tree(io,inference,maxdepth=1)
# # close(io)
# # ProfileView.view(flamegraph(inference))
# inf_exc      = flatten(inference,tmin = 0.0, sortby=exclusive)[end:-1:max(1,end-100)]      # by instance
# # inf_exc      = flatten(inference,tmin = 0.0, sortby=exclusive)[end:-1:1]      # by instance
# #inf_inc      = flatten(inference,tmin = 1.0, sortby=inclusive)[end:-1:max(1,end-49)]      # by instance
# # acced = accumulate_by_source(flat; tmin = 1., by=inclusive) # by method (buggy - repeats methods)



eigincXU          = solve(EigXU{OX,OU};Δω, p, nmod,initialstate,verbose=true,verbosity=1,tol=1e-20,σₓᵤ)

nα                = 32
α                 = 2π*(1:nα)/nα
circle            = 0.05*[cos.(α) sin.(α)]'
GUI(initialstate,eigincXU;shadow = (;EulerBeam3D              = (;style=:shape,line_color=:grey,Udof=false),
                                    StrainGaugeOnEulerBeam3D  = (;gauge_color=:transparent)         ),
                          model  = (;EulerBeam3D              = (;style=:solid,section=circle),
                                     StrainGaugeOnEulerBeam3D = (;L=0.03),
                                     Position3D               = (;L=.03)) ) 


# using SnoopCompileCore, SnoopCompile, AbstractTrees, ProfileView
# using Profile
# using BenchmarkTools

# mission = :profile
# if mission == :report
#     out = solve(EigXU{OX,OU};Δω, p, nmod,initialstate,verbose=true,verbosity=1,tol=1e-20,σₓᵤ)
# elseif mission == :time
#     out = solve(EigXU{OX,OU};Δω, p, nmod,initialstate,verbose=true,verbosity=1,tol=1e-20,σₓᵤ)
#     @btime out =solve(EigXU{OX,OU};Δω, p, nmod,initialstate,verbose=false,verbosity=1,tol=1e-20,σₓᵤ)
# elseif mission == :profile
#     out = solve(EigXU{OX,OU};Δω, p, nmod,initialstate,verbose=true,verbosity=1,tol=1e-20,σₓᵤ)
#     Profile.clear()
#     Profile.@profile for i=1:30
#         local out = solve(EigXU{OX,OU};Δω, p, nmod,initialstate,verbose=false,verbosity=1,tol=1e-20,σₓᵤ)
#     end
#     ProfileView.view(fontsize=30);
#     # After clicking on a bar in the flame diagram, you can type warntype_last() and see the result of 
#     # code_warntype for the call represented by that bar.
# end
;



