#module TestDirectXUA

using Test
using Muscade

###

using StaticArrays
struct SdofOscillator <: AbstractElement
    K :: 𝕣
    C :: 𝕣
    M :: 𝕣
end
SdofOscillator(nod::Vector{Node};K::𝕣,C::𝕣,M::𝕣) = SdofOscillator(K,C,M)
@espy function Muscade.residual(o::SdofOscillator, X,U,A, t,SP,dbg) 
    x,x′,x″,u,ΞC,ΞM = ∂0(X)[1], ∂1(X)[1], ∂2(X)[1], ∂0(U)[1], A[1], A[2]
    r         = -u + o.K*x + o.C*exp10(ΞC)*x′ + o.M*exp10(ΞM)*x″
    return SVector(r),noFB
end
Muscade.doflist( ::Type{SdofOscillator})  = (inod =(1 ,1 ,1 ,1), class=(:X,:U,:A,:A), field=(:x,:u,:ΞC,:ΞM))

###

model1          = Model(:TrueModel)
n               = addnode!(model1,𝕣[ 0, 0])  
e               = addelement!(model1,SdofOscillator,[n], K=1.,C=0.05,M=1.)
#initialstate    = Muscade.State{1,3,1}(initialize!(model1;time=0.))  # recast to force the state to have 2nd derivatives, 
initialstate    = initialize!(model1;nXder=2,time=0.)  # recast to force the state to have 2nd derivatives, 
setdof!(initialstate,[1.];field=:x,nodID=[n],ider=1)
time            = 0:.1:10
state1          = solve(SweepX{0};  initialstate,time,verbose=true)
x               = getdof(state1;field=:x,nodID=[n] )
x               = reshape(x,length(x))

###

using GLMakie




# stateXUA           = solve(DirectXUA;initialstate=stateX,verbose=false)
# stateXUA           = solve(DirectXUA;initialstate=stateX,saveiter=true,verbose=false)

#end 

;