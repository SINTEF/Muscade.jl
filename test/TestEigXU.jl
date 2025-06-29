# TODO
# - normalize X dofs to maximum(abs.(X)) = 1 with the exception of reaction forces?
# - what does λ mean? Report it, as a function of ω and imod.
# - do we indeed have convergence of geneig?  Introduce an option to get geneig toself check.


using Test
using Muscade
using StaticArrays,SparseArrays

q          = 6
nel        = 2^q
Δi         = 2^(q-2) # 8-1=7 sensor stations
isensor    = Δi:Δi:(nel-Δi)
xn         = 1.
λn         = 1e3
un         = 1e3


## simply supported beam
include("../examples/BeamElements.jl")

L    = 1;    # Beam length [m]
q    = 0.0;  # Uniform lateral load [N/m]
EI₂  = 1;    # Bending stiffness [Nm²]
EI₃  = 10;    # Bending stiffness [Nm²]
EA   = 1e6;  # Axial stiffness [N]
GJ   = 1e6;  # Torsional stiffness [Nm²]
μ    = 1;
ι₁   = 1;
hasU = true
column(v::Vector) = reshape(v,(length(v),1))
row(   v::Vector) = reshape(v,(1,length(v)))
vec(   v::Muscade.NodID ) = SVector{1}(v)

XnodeCoord  = hcat(range(0,L,length=nel+1),zeros(nel+1,2))
UnodeCoord  = zeros(nel,3)
mat         = BeamCrossSection(;EA,EI₂,EI₃,GJ,μ,ι₁)
model       = Model(:TestModel)
Xnod        = addnode!(model,XnodeCoord)
Unod        = addnode!(model,UnodeCoord)
mesh        = hcat(Xnod[1:end-1],Xnod[2:end],Unod)
eleid       = addelement!(model,EulerBeam3D{hasU},mesh;mat=mat,orient2=SVector(0.,1.,0.))
[addelement!(model, Hold         , vec(Xnod[  1])           ; field) for field∈(:t1,:t2,:t3,:r1)] 
[addelement!(model, Hold         , vec(Xnod[end])           ; field) for field∈(:t1,:t2,:t3    )] 

# Ucost
addelement!( model, SingleDofCost, column(Unod)    ,class=:U, field=:t1,cost=QuadraticFunction(0.,.1 ))
addelement!( model, SingleDofCost, column(Unod)    ,class=:U, field=:t2,cost=QuadraticFunction(0.,.5 ))
addelement!( model, SingleDofCost, column(Unod)    ,class=:U, field=:t3,cost=QuadraticFunction(0.,1. ))
# disp meas
addelement!( model, SingleDofCost, column(Xnod[isensor])    ,class=:X, field=:t2,cost=QuadraticFunction(0.,1.01))
addelement!( model, SingleDofCost, column(Xnod[isensor])    ,class=:X, field=:t3,cost=QuadraticFunction(0.,1.01))

initialstate      = initialize!(model)   
initialstate.time     = 0.
OX,OU                 = 2,0

## N
σₓᵤ     = (X=(t1 =xn,t2 =xn,t3 =xn,r1 =xn,r2 =xn,r3 =xn,
             λt1=λn,λt2=λn,λt3=λn,λr1=λn,λr2=λn,λr3=λn),
          U=(t1 =un,t2 =un,t3 =un,r1 =un,r2 =un,r3 =un))

# dis             = Muscade.Disassembler(model)
# dofgr           = Muscade.allXdofs(model,dis)
# N               = [zeros(getndof(dofgr)),zeros(getndof(dofgr)),zeros(getndof(dofgr))]
# Muscade.makeXUnorm!(N,dofgr,σₓᵤ)
# @show N


## EigXU analysis
# Δω                = 2^-6 *4
# p                 = 13   -2
Δω                = 10.
p                 = 3
nmod              = 5
eiginc            = solve(EigXU{OX,OU};Δω, p, nmod,initialstate,verbose=true,verbosity=1,tol=1e-20,σₓᵤ)

# ## EigX analysis
# eigincX = solve(EigX{ℝ};state=initialstate,nmod=10)
# ωₚ                = eigincX.ω


# ## draw EigX vector
# jmod              = [10]
# A                 = [100] 
# state             = increment{OX}(initialstate,eigincX,jmod,A);
# using GLMakie
# fig      = Figure(size = (500,500))
# display(fig) # open interactive window (gets closed down by "save")
# axe      = Axis3(fig[1,1],title="Test",xlabel="X",ylabel="Y",zlabel="Z",aspect=:data,viewmode=:fit,perspectiveness=.5)
# draw(axe,state;EulerBeam3D=(;style=:shape,nseg=10,marking=true,Uscale))




## draw EigXU vector
jmod              = [1]
A                 = [1000000] 
iω                = 5
Uscale            = 0.0000002
state             = increment{OX}(initialstate,eiginc,iω,jmod,A);

using GLMakie
fig      = Figure(size = (500,500))
display(fig) # open interactive window (gets closed down by "save")
axe      = Axis3(fig[1,1],title="Test",xlabel="X",ylabel="Y",zlabel="Z",aspect=:data,viewmode=:fit,perspectiveness=.5)
draw(axe,state;EulerBeam3D=(;style=:shape,nseg=10,marking=true,Uscale));

# ## draw error graph
# using GLMakie
# fig      = Figure(size = (500,500))
# display(fig) # open interactive window (gets closed down by "save")
# axe      = Axis(fig[1,1],title="Information content",xlabel="ω [rad/s]",ylabel="magnitude of error",yscale=log)
# nω = 2^p
# ω   = range(start=0.,step=Δω,length=nω) 
# nor = 𝕣1(undef,nω)
# λ   = 𝕣1(undef,nω)
# for imod = 1:maximum(eiginc.ncv)
#     for iω= 1:nω
#         if imod≤eiginc.ncv[iω]
#             nor[iω] = eiginc.nor[iω][imod]
#         else
#             nor[iω] = NaN
#         end
#     end
#     scatter!(axe,ω,nor,markersize=2,color=:black)
# end
# nωₚ = findlast(ωₚ .< ω[end])
# scatter!(axe,ωₚ[1:nωₚ],ones(nωₚ))
;

