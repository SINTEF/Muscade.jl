# TODO
# - normalize X dofs to maximum(abs.(X)) = 1 with the exception of reaction forces?
# - model a circle to avoid "inertial singularities"
# - what does λ mean? Report it, as a function of ω and imod.
# - do we indeed have convergence of geneig?  Introduce an option to get geneig toself check.


using Test
using Muscade
using StaticArrays,SparseArrays

q          = 6
nel        = 2^q
Δi         = 2^3-1 # 8-1=7 sensor stations
isensor    = Δi:Δi:(nel-Δi) # leaving a gap in the set of sensors
iels       = 1:nel
xn         = 1.
λn         = 1e3
un         = 1e3

## beam in space
include("../examples/BeamElements.jl")

L    = 1;    # Beam length [m]
q    = 0.0;  # Uniform lateral load [N/m]
EI₂  = 1;    # Bending stiffness [Nm²]
EI₃  = 10;   # Bending stiffness [Nm²]
EA   = 1e6;  # Axial stiffness [N]
GJ   = 1e6;  # Torsional stiffness [Nm²]
μ    = 1;
ι₁   = 1;
hasU = true

α           = iels/nel*2π
XnodeCoord  = hcat(cos.(α),sin.(α),zeros(nel,1))
UnodeCoord  = zeros(nel,3)
mat         = BeamCrossSection(;EA,EI₂,EI₃,GJ,μ,ι₁)
model       = Model(:TestModel)
Xnod        = addnode!(model,XnodeCoord)
Unod        = addnode!(model,UnodeCoord)
mesh        = hcat(Xnod,Xnod[mod_onebased.(iels.+1,nel)],Unod)
addelement!(model,EulerBeam3D{hasU},mesh;mat=mat,orient2=SVector(0.,0.,1.))

# Ucost
addelement!( model, SingleDofCost, Muscade.columnmatrix(Unod         )    ,class=:U, field=:t1,cost=QuadraticFunction(0.,.1  ))
addelement!( model, SingleDofCost, Muscade.columnmatrix(Unod         )    ,class=:U, field=:t2,cost=QuadraticFunction(0.,.1  ))
addelement!( model, SingleDofCost, Muscade.columnmatrix(Unod         )    ,class=:U, field=:t3,cost=QuadraticFunction(0.,.1  ))
# disp meas
addelement!( model, SingleDofCost, Muscade.columnmatrix(Xnod[isensor])    ,class=:X, field=:t1,cost=QuadraticFunction(0.,1.))
addelement!( model, SingleDofCost, Muscade.columnmatrix(Xnod[isensor])    ,class=:X, field=:t2,cost=QuadraticFunction(0.,1.))
addelement!( model, SingleDofCost, Muscade.columnmatrix(Xnod[isensor])    ,class=:X, field=:t3,cost=QuadraticFunction(0.,1.))

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
#eigincXU          = solve(EigXU{OX,OU};Δω, p, nmod,initialstate,verbose=true,verbosity=1,tol=1e-20,σₓᵤ)

α                 = 2π*(1:8)/8
circle            = 0.05*[cos.(α) sin.(α)]'
draw(initialstate,eigincXU;shadow = (;EulerBeam3D=(;style=:shape,nseg=1,line_color=:grey,Udof=false)),
                           model  = (;EulerBeam3D=(;style=:solid,nseg=1,section=circle))           )

# ### draw ###########################################
# using GLMakie

# ## Organize the figure
# fig      = Figure(size = (1500,900),title="EigXU results")
# nω = 2^p
# ω   = range(start=0.,step=Δω,length=nω) 


# display(fig) # open interactive window (gets closed down by "save")

# panelFreqs = fig[1,1]        
# panelNorm  = panelFreqs[1,1] 
# axisNorm   = Axis(panelNorm,xlabel="ω [rad/s]",ylabel="magnitude of error",yscale=log10)
# panelSlide = panelFreqs[2,1] 
# panelModel = fig[1,2:3]        
# Box(panelModel, cornerradius = 20,z=1., color = :transparent)
# axisModel  = Axis3(panelModel,title="EigXU mode shape",aspect=:data,viewmode=:free,perspectiveness=.5,clip=false)

# ## sliders
# sg = SliderGrid(panelSlide,
#                 (label="ω"      , range = ω        , startvalue = 0,snap=true,update_while_dragging=true,format = "{:.1f} rad/s"),
#                 (label="mode"   , range = 1:nmod   , startvalue = 1,snap=true,update_while_dragging=true                        ),
#                 (label="X scale", range = -5:0.01:5, startvalue = 0,snap=true,update_while_dragging=true,format = "10^{:.1f}"   ),
#                 (label="U scale", range = -5:0.01:5, startvalue = 0,snap=true,update_while_dragging=true,format = "10^{:.1f}"   ))
# obs = (ω      = sg.sliders[1].value,
#        imode  = sg.sliders[2].value,
#        Xscale = sg.sliders[3].value,
#        Uscale = sg.sliders[4].value)


# ## norm spectre
# nor = 𝕣1(undef,nω)
# λ   = 𝕣1(undef,nω)
# for imod = 1:maximum(eigincXU.ncv)
#     for iω= 1:nω
#         if imod≤eigincXU.ncv[iω]
#             nor[iω] = eigincXU.nor[iω][imod]
#             λ[  iω] = eigincXU.λ[  iω][imod]
#         else
#             nor[iω] = NaN
#             λ[  iω] = NaN
#         end
#     end
#     scatter!(axisNorm,ω,nor,markersize=1,color=:black)
# end

# iω  = map(obs.ω) do ω
#     round(Int64,ω/Δω)+1
# end       
# nor = map(obs.imode,iω) do imode,iω 
#     eigincXU.nor[iω][imode]
# end

# scatter!(axisNorm,obs.ω,nor,color=:red,markersize=10)

# ## Model
#     α      = 2π*(1:8)/8
#     circle = 0.05*[cos.(α) sin.(α)]'


# draw!(          axisModel,initialstate;EulerBeam3D=(;style=:shape,nseg=1,line_color=:grey,Udof=false)); # draw initial state once to keep on screen

# graphic = draw!(axisModel,initialstate;EulerBeam3D=(;style=:solid,nseg=1,section=circle)); # and twice to start the pump
# _ = map(iω,obs.imode,obs.Xscale,obs.Uscale) do iω,imod,Xscale,Uscale # Then observe the sliders
#     state = Muscade.visualincrement(initialstate,eigincXU,iω,imod;Xscale=exp10(Xscale),Uscale=exp10(Uscale))
#     draw!(graphic,state;EulerBeam3D=(;style=:solid,nseg=1,section=circle));
# end

;
