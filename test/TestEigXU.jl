# TODO
# - normalize X dofs to maximum(abs.(X)) = 1 with the exception of reaction forces?
# - what does λ mean? Report it, as a function of ω and imod.
# - do we indeed have convergence of geneig?  Introduce an option to get geneig toself check.

#module TestFreqXU

using Test
using Muscade
using StaticArrays,SparseArrays

include("../examples/BeamElements.jl")

L    = 1;    # Beam length [m]
q    = 0.0;  # Uniform lateral load [N/m]
EI₂  = 1;    # Bending stiffness [Nm²]
EI₃  = 1;    # Bending stiffness [Nm²]
EA   = 1e6;  # Axial stiffness [N]
GJ   = 1e6;  # Torsional stiffness [Nm²]
μ    = 1;
ι₁   = 1;
hasU = true
column(v::Vector) = reshape(v,(length(v),1))
row(   v::Vector) = reshape(v,(1,length(v)))
vec(   v::Muscade.NodID ) = SVector{1}(v)

nel         = 10
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
addelement!( model, SingleDofCost, column(Unod)    ,class=:U, field=:t1,cost=QuadraticFunction(0.,.1 ))
addelement!( model, SingleDofCost, column(Unod)    ,class=:U, field=:t2,cost=QuadraticFunction(0.,.5 ))
addelement!( model, SingleDofCost, column(Unod)    ,class=:U, field=:t3,cost=QuadraticFunction(0.,1. ))
addelement!( model, SingleDofCost, vec(Xnod[4])    ,class=:X, field=:t2,cost=QuadraticFunction(0.,10.))
addelement!( model, SingleDofCost, vec(Xnod[7])    ,class=:X, field=:t3,cost=QuadraticFunction(0.,10.))

initialstate      = initialize!(model)   
initialstate.time = 0.

OX,OU             = 2,0
if true # eigXU
    Δω                = 2^-6
    p                 = 13
    nmod              = 2
    eiginc            = solve(EigXU{OX,OU};Δω, p, nmod,initialstate,verbose=true,verbosity=1,tol=1e-20)

    # α = 2π*(0:19)/20
    # circle = 0.05*[cos.(α) sin.(α)]'
    # jmod              = [3]
    # A                 = [1000] 
    # iω                = 11
    # Uscale            = 0.002
    # state             = increment{OX}(initialstate,eiginc,iω,jmod,A)
    # using GLMakie
    # fig      = Figure(size = (500,500))
    # display(fig) # open interactive window (gets closed down by "save")
    # axe      = Axis3(fig[1,1],title="Test",xlabel="X",ylabel="Y",zlabel="Z",aspect=:data,viewmode=:fit,perspectiveness=.5)
    # draw(axe,state;EulerBeam3D=(;style=:shape,nseg=10,section = circle,marking=true,Uscale))
end
if true # eigX
    eigincX = solve(EigX{ℝ};state=initialstate,nmod=10)
    # jmod              = [1]
    # A                 = [1] 
    # state             = increment(initialstate,eigincX,jmod,A)
    ωₚ                = eigincX.ω
end
using GLMakie
fig      = Figure(size = (500,500))
display(fig) # open interactive window (gets closed down by "save")
axe      = Axis(fig[1,1],title="Information content",xlabel="ω [rad/s]",ylabel="S [bit/s]",yscale=log)


nω = 2^p
ω   = range(start=0.,step=Δω,length=nω) 
S = 𝕣1(undef,nω)
λ = 𝕣1(undef,nω)
for imod = 1:maximum(eiginc.ncv)
    for iω= 1:nω
        if imod≤eiginc.ncv[iω]
            S[iω] = eiginc.S[iω][imod]
        else
            S[iω] = NaN
        end
    end
    scatter!(axe,ω,S,markersize=2,color=:black)
end
nωₚ = findlast(ωₚ .< ω[end])
scatter!(axe,ωₚ[1:nωₚ],ones(nωₚ))


# dof             = getdof(state,field=:t1)
# @testset "output" begin
#     @test dof[1,1:100:end]      == [-0.0008522224876621403, 0.09549942165460336, -0.05008960568585242,  -0.06878034794799491, 0.08677879820843742, 0.02249032282693445,  -0.09877571024378282, 0.030199161274387067, 0.08266670708994063, -0.07429569207304236, -0.043035489287326145]
#     @test dof[2,1:100:end]      == [  0.09808063976613206,    -0.026980709523268566,    -0.08368844089497145,     0.07162225920301411,     0.045483300397735474,    -0.09588421213929227,     0.005663784838222031,     0.09286300479539936,    -0.05519928099099323,    -0.06341829990791986,     0.0890282202766159]
# end

#end 

