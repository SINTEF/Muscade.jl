module TestElementCostAndConstraint

iλ,ix,iu,ia=1,2,3,4

using Test,StaticArrays
using Muscade
include("SomeElements.jl")

model           = Model(:TestModel)
n1              = addnode!(model,𝕣[0,0,+100]) # turbine
n2              = addnode!(model,𝕣[])         # Anod for anchor
@functor with() cost(eleres,t) =            eleres.Fh^2
@functor with() gap( eleres,t) = SVector{1}(eleres.Fh^2)
@functor with() mode(       t) = SVector{1}(:equal     )
req                            = @request(Fh)

model_cost_and_Xconstraint = deepcopy(model)
model_Xconstraint          = deepcopy(model)
model_naked                = deepcopy(model)

OX,OU,IA        = 0,0,1
Nx,Nu,Na        = 3+1,0,2   
λ               =  SVector{Nx,𝕣}(1. for i=1:Nx)
x               = (SVector{Nx,𝕣}(1. for i=1:Nx),)
u               = (SVector{Nu,𝕣}(1. for i=1:Nu),)
a               =  SVector{Na,𝕣}(0. for i=1:Na)
_,Nz,(Λ,X,U,A)  = variate((λ,x,u,a))

Nxn,Nun,Nan     = 3,0,2   
λn              =  SVector{Nxn,𝕣}(1. for i=1:Nxn)
xn              = (SVector{Nxn,𝕣}(1. for i=1:Nxn),)
un              = (SVector{Nun,𝕣}(1. for i=1:Nun),)
an              =  SVector{Nan,𝕣}(0. for i=1:Nan)
_,Nzn,(Λn,Xn,Un,An)  = variate((λn,xn,un,an))


## ################ no constraint, no cost (naked)

∇L_naked   = [-438861.13074456755, 9278.602091074139, 1.871510789932893e6, -86753.51732058081, 8308.805033394932, 1.5339992031591167e6, -2.5203831430346757e7, 1.4419282612793995e6]
r_naked    = ∇L_naked[SVector(1,2,3)] 

el_naked   = ElementCostAndConstraint(model_naked.nod;TargetElement=AnchorLine, elementkwargs=(Δxₘtop=[5.,0,0], xₘbot=[250.,0], L=290., buoyancy=-5e3),req=@request(Fh) ) 
L,FB       = Muscade.lagrangian(el_naked, Λn,Xn,Un,An, 0.,nothing,(;naked=true))                 
@testset "lagrangian naked" begin
     @test ∂{1,Nzn}(L) ≈ ∇L_naked
end

## Test addin! DirectXUA with naked

e1                             = addelement!(model_naked,ElementCostAndConstraint,[n1,n2];  
                              TargetElement=AnchorLine,elementkwargs=(Δxₘtop=[5.,0,0], xₘbot=[250.,0], L=290., buoyancy=-5e3),req=@request(Fh))
initialstate                   = initialize!(model_naked)

wanted                         = Muscade.Wanted{1,OX+1,OU+1,IA}(:all,:all)
model_naked,dis                = initialstate.model, initialstate.dis
out,asm,dofgr                  = Muscade.prepare(Muscade.AssemblyDirect,model_naked,dis,wanted)          # mem and assembler for system at any given step
Muscade.zero!(out)
ieletyp,iele                   = 1,1
asmᵢ                           = asm[:,ieletyp]
scale                          = dis.dis[ieletyp].scale
t,Δt                           = 0.,1.
SP                             = nothing
Muscade.addin!{:matrices}(out,asmᵢ,iele,scale,model_naked.eleobj[1][1],Val(true),(λn,),xn,un,an,t,Δt,SP,(;naked=true))

@testset "addin! DirectXUA naked" begin
     @test out.L1[iλ][1] ≈ ∇L_naked[1:3]
end



## ################# With X-constraint, without cost

∇L_Xconstraint = [-1.926851845351649e11, 2.3221912378082836e10, -4.909682577367074e8, -9.902918379835745e10, -7.62996816062726e10, 3.3129419089208664e10, -1.0430356765960333e9, -2.028735837473837e11, 9.611778704183535e12, -5.379682903546137e11]
r_Xconstraint  = ∇L_Xconstraint[SVector{4}(1:4)]
K_Xconstraint  = [ -0.0            2.32224e10    -4.90978e8    -9.90311e10   0.0;  # the famous fifth column: corresponding to δr of Newmark-β adiffing
                   2.32224e10      26446.5       -521.472      -1.12679e5    0.0;
                  -4.90978e8      -521.472        1792.85       7037.43      0.0;
                  -9.90311e10     -1.12679e5      7037.43       1.63964e6    0.0]


## lagrangian with Xconstraint


el_Xconstraint            = ElementCostAndConstraint(model_Xconstraint.nod;TargetElement=AnchorLine, elementkwargs=(Δxₘtop=[5.,0,0], xₘbot=[250.,0], L=290., buoyancy=-5e3),
                                           req=@request(Fh), gap, λxinod=(1,), λxfield=(:λ,), mode) 
L,FB            = Muscade.lagrangian(el_Xconstraint, Λ,X,U,A, 0.,nothing,(Xconstraint=true,))                 

@testset "lagrangian Xconstraint" begin
     @test ∂{1,Nz}(L) ≈ ∇L_Xconstraint
end

## Test Xdiffedresidual for SweepX with Xconstraint

_,Nz,(X,δr)  = variate((x,0.))
R,FB = Muscade.Xdiffedresidual(el_Xconstraint,X,u,a,0.,nothing,(Xconstraint=true,))
r,K = value_∂{1,Nz}(R)



@testset "Xdiffedresidual Xconstraint" begin
   @test r ≈ r_Xconstraint
   @test K ≈ K_Xconstraint rtol = 1e-3
end

## Test addin! for directXUA with Xconstraint

e1                             = addelement!(model_Xconstraint,ElementCostAndConstraint,[n1,n2];  
                              TargetElement=AnchorLine,elementkwargs=(Δxₘtop=[5.,0,0], xₘbot=[250.,0], L=290., buoyancy=-5e3),
                              req, gap, λxinod=(1,), λxfield=(:λ,), mode)
initialstate                   = initialize!(model_Xconstraint)

wanted                         = Muscade.Wanted{1,OX+1,OU+1,IA}(:all,:all)
model_Xconstraint,dis          = initialstate.model, initialstate.dis
out,asm,dofgr                  = Muscade.prepare(Muscade.AssemblyDirect,model_Xconstraint,dis,wanted)          # mem and assembler for system at any given step
Muscade.zero!(out)
ieletyp,iele                   = 1,1
asmᵢ                           = asm[:,ieletyp]
scale                          = dis.dis[ieletyp].scale
t,Δt                           = 0.,1.
SP                             = nothing
Muscade.addin!{:matrices}(out,asmᵢ,iele,scale,model_Xconstraint.eleobj[1][1],Val(true),(λ,),x,u,a,t,Δt,SP,(;Xconstraint=true))

 @testset "addin! DirectXUA Xconstraint" begin
     @test out.L1[iλ][1] ≈ r_Xconstraint
     @test out.L1[ix][1] ≈ [0.0, 0.0, 0.0, 0.0]
     @test out.L1[iu][1] ≈ Float64[]
     @test out.L1[ia][1] ≈ [0.0, 0.0]

     @test out.L2[iλ,ix][1,1] ≈ K_Xconstraint[:,1:4] rtol = 1e-6
     @test out.L2[ix,iλ][1,1] ≈ K_Xconstraint[:,1:4]' rtol = 1e-6
     
     @test out.L2[iλ,iu][1,1] ≈ Matrix{Float64}(undef, 4, 0)
     @test out.L2[iu,iλ][1,1] ≈ Matrix{Float64}(undef, 0, 4)
     @test out.L2[iλ,ia][1,1] ≈ [6.73599e12    -3.8537e11; 7.67097e6 -4.38861e5; -162183.0 9278.6; -3.27126e7      1.87151e6] rtol = 1e-5
     @test out.L2[ia,iλ][1,1] ≈ out.L2[iλ,ia][1,1]'

     @test out.L2[ix,ix][1,1] ≈ [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0]

     @test out.L2[iu,ix][1,1] ≈ Matrix{Float64}(undef, 0, 4)
     @test out.L2[ix,iu][1,1] ≈ Matrix{Float64}(undef, 4, 0)

     @test out.L2[ia,ix][1,1] ≈ [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0]
     @test out.L2[ix,ia][1,1] ≈ [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]
end



## ############# With X-constraint and cost

## Test Lagrangian 

el_cost_and_Xconstraint  = ElementCostAndConstraint(model_cost_and_Xconstraint.nod;TargetElement=AnchorLine, elementkwargs=(Δxₘtop=[5.,0,0], xₘbot=[250.,0], L=290., buoyancy=-5e3),
                req, gap, λxinod=(1,), λxfield=(:λ,), mode, cost)
Nx,Nu,Na        = 3+1,0,2   
λ               =  SVector{Nx,𝕣}(1. for i=1:Nx)
x               = (SVector{Nx,𝕣}(1. for i=1:Nx),)
u               = (SVector{Nu,𝕣}(1. for i=1:Nu),)
a               =  SVector{Na,𝕣}(0. for i=1:Na)
_,Nz,(Λ,X,U,A)  = variate((λ,x,u,a))
dbg             = (alwaystesteverything=true,)
L,FB            = Muscade.lagrangian(el_cost_and_Xconstraint, Λ,X,U,A, 0.,nothing,(Cost_Xconstraint=true,))                 

@testset "lagrangian" begin
     @test Muscade.doflist(typeof(el_cost_and_Xconstraint)) == (inod = (1, 1, 1, 1, 2, 2), class = (:X, :X, :X, :X, :A, :A), field = (:λ, :tx1, :tx2, :rx3,  :ΔL, :Δbuoyancy))
     @test value{1}(L) ≈ -7.629823967801132e10
     @test ∂{1,Nz}(L) ≈ [-1.926851845351649e11, 2.3221912378082836e10, -4.909682577367074e8, -9.902918379835745e10, -7.62996816062726e10, 9.907067849995083e9, -5.520581402572349e8, -1.0384252843823631e11, 2.875791844697829e12, -1.5259792128428394e11]
end

req            = @request λ,eleres(cr,Fh)
L,FB,eleres    = Muscade.lagrangian(el_cost_and_Xconstraint, λ, x,u,a, 0.,nothing,(cost_Xconstraint=true,),req)   
@testset "lagrnagian results" begin
     @test eleres.λ         ≈ [1.]
     @test eleres.eleres.cr ≈ 87.79184120068672
     @test eleres.eleres.Fh ≈ 438959.2060034336
end


## Test draw

nel       = 1
nXder     = nUder = 1
EL        = [el_cost_and_Xconstraint]
Λm        = 𝕣2(undef,Nx,nel)
Xm        = ntuple(i->𝕣2(undef,Nx,nel),nXder)
Um        = ntuple(i->𝕣2(undef,Nu,nel),nUder)
Am        = 𝕣2(undef,Na,nel)
for i ∈ eachindex(EL)
    Λm[:,i]    .= λ
    Xm[1][:,i] .= x[1]
    Um[1][:,i] .= u[1]
    Am[:,i]    .= a
end

using Muscade: lines!,scatter!,mesh!,SpyAxis
axe = SpyAxis()
Muscade.draw_element!(axe,EL, Λm,Xm,Um,Am, 0.,nothing,(;))
@testset "drawing" begin
     @test  axe.call[1].fun == :lines!
     @test  axe.call[1].args[1][][:,1:11] ≈ [ 126.035    113.802     101.568    89.3348   77.1015  64.8682   52.6348   40.4015   28.1682   15.9348    3.70151  ;
                                                2.62093    2.87957     3.13821   3.39686   3.6555   3.91414   4.17278   4.43143   4.69007   4.94871   5.20735  ;
                                                0.0        0.854088    3.43297   7.78682  14.0004  22.1945   32.5286   45.2038   60.4668   78.6144   99.9998    ] atol=0.001
     @test  axe.call[1].kwargs[:color] == :blue
     @test  axe.call[1].kwargs[:linewidth] == 2
     @test  length(axe.call) == 5
end


############# Single dof oscillator

## Test addin!(AssemblyDirect,...) with AdjustableSdofOscillator

req                            = @request((damping,C))
@functor with() cost(eleres,t) =            eleres.damping^2
@functor with() gap( eleres,t) = SVector{2}(eleres.damping^2,eleres.C) # λX,λU
@functor with() mode(       t) = SVector{2}(:equal          ,:equal  )

model                 = Model()
n1                    = addnode!(model,[0.])
e1                    = addelement!(model,ElementCostAndConstraint,[n1];  
                              TargetElement=AdjustableSdofOscillator,elementkwargs=(K=1.,C=.1,M=1.),
                              req, gap, λxinod=(1,), λxfield=(:λx,), λuinod=(1,), λufield=(:λu,), mode, cost)
initialstate          = initialize!(model)

OX,OU,IA              = 2,0,1
Nx                    = 2
Nu                    = 2
Na                    = 1
λ                     =             SVector{Nx,𝕣}(1. for i=1:Nx)
x                     = ntuple(j -> SVector{Nx,𝕣}(1. for i=1:Nx), OX+1)
u                     = ntuple(j -> SVector{Nu,𝕣}(1. for i=1:Nu), OU+1)
a                     =             SVector{Na,𝕣}(0. for i=1:Na)

wanted                = Muscade.Wanted{1,OX+1,OU+1,IA}(:all,:all)
model,dis             = initialstate.model, initialstate.dis
out,asm,dofgr         = Muscade.prepare(Muscade.AssemblyDirect,model,dis,wanted)          # mem and assembler for system at any given step
Muscade.zero!(out)
ieletyp,iele          = 1,1
asmᵢ                  = asm[:,ieletyp]
scale                 = dis.dis[ieletyp].scale
t,Δt                  = 0.,1.
SP                    = nothing
Muscade.addin!{:matrices}(out,asmᵢ,iele,scale,model.eleobj[1][1],Val(true),(10λ,),x,u,a,t,Δt,SP,(;testall=true))

iλ,ix,iu,ia=1,2,3,4


@testset "addin! DirectXUA constraint on AdjustableSdofOscillator" begin
     @test out.L1[iλ][1] ≈  [-0.01, 1.1] rtol = 1e-5
     @test out.L1[ix][1] ≈  [0.0, 0.0]
     @test out.L1[ix][2] ≈  [0.0, 0.02] 
     @test out.L1[ix][3] ≈  [0,0]
     @test out.L1[iu][1] ≈  [0.0, 0.0]
     @test out.L1[ia][1] ≈  [0.046051701859880924] rtol = 1e-5

     @test out.L2[iλ,ix][1,1] ≈ [0 0;0 1]
     @test out.L2[iλ,ix][1,2] ≈ [0 -0.02;0 .1] rtol = 1e-5
     @test out.L2[iλ,ix][1,3] ≈ [0 0;0 1] rtol = 1e-5
     @test out.L2[ix,iλ][1,1] ≈ [0 0;0 1]' rtol = 1e-5
     @test out.L2[ix,iλ][2,1] ≈ [0 0;-0.02 .1] rtol = 1e-5
     @test out.L2[ix,iλ][3,1] ≈ [0 0;0 1]' rtol = 1e-5

     @test out.L2[iλ,iu][1,1] ≈ [0 0;0 -1] rtol = 1e-5
     @test out.L2[iu,iλ][1,1] ≈ [0 0;0 -1]' rtol = 1e-5
     
     @test out.L2[iλ,ia][1,1] ≈ [-0.046051701859880924,  0.2302585092994046] rtol = 1e-5
     @test out.L2[ia,iλ][1,1] ≈ [-0.046051701859880924   0.2302585092994046] rtol = 1e-5

     @test out.L2[ix,ix][1,1] ≈ [0 0;0 0]
     @test out.L2[ix,ix][1,2] ≈ [0 0.;0 0]
     @test out.L2[ix,ix][1,3] ≈ [0 0;0 0]
     @test out.L2[ix,ix][2,1] ≈ [0 0;0 0]
     @test out.L2[ix,ix][2,2] ≈ [0 0;0 .02] rtol = 1e-5
     @test out.L2[ix,ix][2,3] ≈ [0 0;0 0]
     @test out.L2[ix,ix][3,1] ≈ [0 0;0 0]
     @test out.L2[ix,ix][3,2] ≈ [0 0;0 0]
     @test out.L2[ix,ix][3,3] ≈ [0 0;0 0]
 
     @test out.L2[iu,ix][1,1] ≈ [0 0;0 0]
     @test out.L2[iu,ix][1,2] ≈ [0 0;0 0]
     @test out.L2[iu,ix][1,3] ≈ [0 0;0 0]
     @test out.L2[ix,iu][1,1] ≈ [0 0;0 0]
     @test out.L2[ix,iu][2,1] ≈ [0 0;0 0]
     @test out.L2[ix,iu][3,1] ≈ [0 0;0 0]

     @test out.L2[ia,ix][1,1] ≈ [0 0]
     @test out.L2[ia,ix][1,2] ≈ [0.0  0.0460517] rtol = 1e-5
     @test out.L2[ia,ix][1,3] ≈ [0 0]
     @test out.L2[ix,ia][1,1] ≈ [0,0]
     @test out.L2[ix,ia][2,1] ≈ [0.0,  0.0460517] rtol = 1e-5
     @test out.L2[ix,ia][3,1] ≈ [0,0]
end


end
;