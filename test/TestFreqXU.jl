module TestFreqXU

using Test
using Muscade
using StaticArrays,SparseArrays

struct El1 <: AbstractElement
    K :: 𝕣
    C :: 𝕣
    M :: 𝕣
end
El1(nod::Vector{Node};K::𝕣,C::𝕣,M::𝕣) = El1(K,C,M)
@espy function Muscade.residual(o::El1, X,U,A, t,SP,dbg) 
    x,x′,x″,ΞC,ΞM = ∂0(X)[1], ∂1(X)[1], ∂2(X)[1],  A[1], A[2]
    r               = o.K*x + o.C*exp10(ΞC)*x′ + o.M*exp10(ΞM)*x″
    return SVector(r),noFB
end
Muscade.doflist( ::Type{El1})  = (inod =(1 ,1 ,1), class=(:X,:A,:A), field=(:tx1,:ΞC,:ΞM))

include("SomeElements.jl")

model           = Model(:TrueModel)
n1              = addnode!(model,𝕣[0])  
n2              = addnode!(model,𝕣[1])  
n3              = addnode!(model,𝕣[ ]) # anode for spring
e1              = addelement!(model,El1,[n1], K=1.,C=0.05,M=2.)
e2              = addelement!(model,El1,[n2], K=0.,C=0.0 ,M=2.)
e3              = addelement!(model,Spring{1},[n1,n2,n3], EA=1.1)
@functor with() fu(u,t)   = (u/.1)^2/2
@functor with() l1(tx1,t) = ((tx1-0.1*sin(t))/.01)^2/2
@functor with() l2(tx1,t) = ((tx1-0.1*cos(t))/.01)^2/2
e10             = addelement!(model,SingleUdof   ,[n1];Xfield=:tx1,Ufield=:tx1, cost=fu)
e11             = addelement!(model,SingleUdof   ,[n2];Xfield=:tx1,Ufield=:tx1, cost=fu)
e12             = addelement!(model,SingleDofCost,[n1];class=:X,field=:tx1,     cost=l1)
e13             = addelement!(model,SingleDofCost,[n2];class=:X,field=:tx1,     cost=l2)

p               = 10
Δt              = (3*2π)/2^p
t₀              = 0.
OX              = 2
OU              = 2
t               = range(start=t₀,step=Δt,length=2^p)

initialstate    = initialize!(model)   
state           = solve(FreqXU{OX,OU};Δt, p, t₀,initialstate,verbose=false)
dof             = getdof(state,field=:tx1)

# using GLMakie
# fig      = Figure(size = (1000,700))
# display(fig) # open interactive window (gets closed down by "save")
# axe      = Axis(fig[1,1],title="Test",xlabel="Time",ylabel="X-dofs")
# with_theme(Theme(fontsize = 30,font="Arial")) do
#     h1    = lines!(  axe,t,dof[1,:],color = :black, linewidth = 1)
#     h2    = lines!(  axe,t,dof[2,:],color = :red  , linewidth = 1)
# end

@testset "output" begin
    @test dof[1,1:100:end]      ≈ [-0.0008522224876621403, 0.09549942165460336, -0.05008960568585242,  -0.06878034794799491, 0.08677879820843742, 0.02249032282693445,  -0.09877571024378282, 0.030199161274387067, 0.08266670708994063, -0.07429569207304236, -0.043035489287326145] rtol=1e-10
    @test dof[2,1:100:end]      ≈ [  0.09808063976613206,    -0.026980709523268566,    -0.08368844089497145,     0.07162225920301411,     0.045483300397735474,    -0.09588421213929227,     0.005663784838222031,     0.09286300479539936,    -0.05519928099099323,    -0.06341829990791986,     0.0890282202766159] rtol=1e-10
end

end 

