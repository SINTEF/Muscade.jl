#module TestElementCostAndConstraint

using Test,StaticArrays
using Muscade
include("SomeElements.jl")

## Test Lagrangian

model           = Model(:TestModel)
n1              = addnode!(model,𝕣[0,0,+100]) # turbine
n2              = addnode!(model,𝕣[])         # Anod for anchor
@functor with() cost(eleres,t) =            eleres.Fh^2
@functor with() gap( eleres,t) = SVector{1}(eleres.Fh^2)
@functor with() mode(       t) = SVector{1}(:equal     )
elcoco  = ElementCostAndConstraint(model.nod;TargetElement=AnchorLine, elementkwargs=(Δxₘtop=[5.,0,0], xₘbot=[250.,0], L=290., buoyancy=-5e3),
                req=@request(Fh), gap, λxinod=(1,), λxfield=(:λ,), mode, cost)
Nx,Nu,Na        = 3+1,0,2   
λ               =  SVector{Nx,𝕣}(1. for i=1:Nx)
x               = (SVector{Nx,𝕣}(1. for i=1:Nx),)
u               = (SVector{Nu,𝕣}(1. for i=1:Nu),)
a               =  SVector{Na,𝕣}(0. for i=1:Na)
_,Nz,(Λ,X,U,A)  = variate((λ,x,u,a))
dbg             = (alwaystesteverything=true,)
L,FB            = Muscade.lagrangian(elcoco, Λ,X,U,A, 0.,nothing,dbg)                 

@testset "lagrangian" begin
     @test Muscade.doflist(typeof(elcoco)) == (inod = (1, 1, 1, 1, 2, 2), class = (:X, :X, :X, :X, :A, :A), field = (:λ, :tx1, :tx2, :rx3,  :ΔL, :Δbuoyancy))
     @test value{1}(L) ≈ -7.629823967801132e10
     @test ∂{1,Nz}(L) ≈ [-1.926851845351649e11, 2.3221912378082836e10, -4.909682577367074e8, -9.902918379835745e10, -7.62996816062726e10, 9.907067849995083e9, -5.520581402572349e8, -1.0384252843823631e11, 2.875791844697829e12, -1.5259792128428394e11]
end

req            = @request λ,eleres(cr,Fh)
L,FB,eleres    = Muscade.lagrangian(elcoco, λ, x,u,a, 0.,nothing,dbg,req)   
@testset "results" begin
     @test eleres.λ         ≈ [1.]
     @test eleres.eleres.cr ≈ 87.79184120068672
     @test eleres.eleres.Fh ≈ 438959.2060034336
end

## lagrangian with no cost (should compare with Xdiffedresidual)

elcoco  = ElementCostAndConstraint(model.nod;TargetElement=AnchorLine, elementkwargs=(Δxₘtop=[5.,0,0], xₘbot=[250.,0], L=290., buoyancy=-5e3),
                req=@request(Fh), gap, λxinod=(1,), λxfield=(:λ,), mode) #, cost)
L,FB            = Muscade.lagrangian(elcoco, Λ,X,U,A, 0.,nothing,dbg)                 

@testset "lagrangian no cost" begin
     @test ∂{1,Nz}(L) ≈ [-1.926851845351649e11, 2.3221912378082836e10, -4.909682577367074e8, -9.902918379835745e10, -7.62996816062726e10, 3.3129419089208664e10, -1.0430356765960333e9, -2.028735837473837e11, 9.611778704183535e12, -5.379682903546137e11]
end



## Test Xdiffedresidual for SweepX


_,Nz,(X,δr)  = variate((x,0.))
R,FB = Muscade.Xdiffedresidual(elcoco,X,u,a,0.,nothing,dbg)
r,K = value_∂{1,Nz}(R)

@testset "Xdiffedresidual" begin
   @test r ≈ [ -1.926851845351649e11,  2.3221912378082836e10, -4.909682577367074e8, -9.902918379835745e10]
   @test r ≈ ∂{1,10}(L)[1:4]  # from lagrangian no cost, so "lagrangian" and "Xdiffedresidual" agree on target-residual and X-constraint
   @test K ≈ [ -0.0             2.32224e10    -4.90978e8    -9.90311e10  -0.0
                2.32224e10      26446.5       -521.472      -1.12679e5    0.0
               -4.90978e8      -521.472        1792.85       7037.43      0.0
               -9.90311e10     -1.12679e5      7037.43       1.63964e6    0.0] rtol = 1e-3
end



# ## Test draw

# nel       = 1
# nXder     = nUder = 1
# EL        = [elcoco]
# Λm        = 𝕣2(undef,Nx,nel)
# Xm        = ntuple(i->𝕣2(undef,Nx,nel),nXder)
# Um        = ntuple(i->𝕣2(undef,Nu,nel),nUder)
# Am        = 𝕣2(undef,Na,nel)
# for i ∈ eachindex(EL)
#     Λm[:,i]    .= λ
#     Xm[1][:,i] .= x[1]
#     Um[1][:,i] .= u[1]
#     Am[:,i]    .= a
# end

# using Muscade: lines!,scatter!,mesh!,SpyAxis
# axe = SpyAxis()
# Muscade.draw_element!(axe,EL, Λm,Xm,Um,Am, 0.,nothing,(;))
# @testset "drawing" begin
#      @test  axe.call[1].fun == :lines!
#      @test  axe.call[1].args[1][][:,1:11] ≈ [ 126.035    113.802     101.568    89.3348   77.1015  64.8682   52.6348   40.4015   28.1682   15.9348    3.70151  ;
#                                                 2.62093    2.87957     3.13821   3.39686   3.6555   3.91414   4.17278   4.43143   4.69007   4.94871   5.20735  ;
#                                                 0.0        0.854088    3.43297   7.78682  14.0004  22.1945   32.5286   45.2038   60.4668   78.6144   99.9998    ] atol=0.001
#      @test  axe.call[1].kwargs[:color] == :blue
#      @test  axe.call[1].kwargs[:linewidth] == 2
#      @test  length(axe.call) == 5
# end

# ## Test addin!(AssemblyDirect,...)

# req                            = @request((damping,C))
# @functor with() cost(eleres,t) =            eleres.damping^2
# @functor with() gap( eleres,t) = SVector{2}(eleres.damping^2,eleres.C) # λX,λU
# @functor with() mode(       t) = SVector{2}(:equal          ,:equal  )

# model                 = Model()
# n1                    = addnode!(model,[0.])
# e1                    = addelement!(model,ElementCostAndConstraint,[n1];  
#                               TargetElement=AdjustableSdofOscillator,elementkwargs=(K=1.,C=.1,M=1.),
#                               req, gap, λxinod=(1,), λxfield=(:λx,), λuinod=(1,), λufield=(:λu,), mode, cost)
# initialstate          = initialize!(model)

# OX,OU,IA              = 2,0,1
# Nx                    = 2
# Nu                    = 2
# Na                    = 1
# λ                     =             SVector{Nx,𝕣}(1. for i=1:Nx)
# x                     = ntuple(j -> SVector{Nx,𝕣}(1. for i=1:Nx), OX+1)
# u                     = ntuple(j -> SVector{Nu,𝕣}(1. for i=1:Nu), OU+1)
# a                     =             SVector{Na,𝕣}(0. for i=1:Na)

# wanted                = Muscade.Wanted{1,OX+1,OU+1,IA}(:all,:all)
# model,dis             = initialstate.model, initialstate.dis
# out,asm,dofgr         = Muscade.prepare(Muscade.AssemblyDirect,model,dis,wanted)          # mem and assembler for system at any given step
# Muscade.zero!(out)
# ieletyp,iele          = 1,1
# asmᵢ                  = asm[:,ieletyp]
# scale                 = dis.dis[ieletyp].scale
# t,Δt                  = 0.,1.
# SP                    = nothing
# Muscade.addin!{:matrices}(out,asmᵢ,iele,scale,model.eleobj[1][1],Val(true),(10λ,),x,u,a,t,Δt,SP,(;testall=true))

# iλ,ix,iu,ia=1,2,3,4


# @testset "addin! DirectXUA constraint" begin
#      @test out.L1[iλ][1] ≈  [0.0, 1.1] 
#      @test out.L1[ix][1] ≈  [-0.01,0]
#      @test out.L1[ix][2] ≈  [0,0.00]
#      @test out.L1[ix][3] ≈  [0,0]
#      @test out.L1[iu][1] ≈  [-.1,0]
#      @test out.L1[ia][1] ≈  [-0.23025850929940464]

#      @test out.L2[iλ,ix][1,1] ≈ [0 0;0 1]
#      @test out.L2[iλ,ix][1,2] ≈ [0 0;0 .1]
#      @test out.L2[iλ,ix][1,3] ≈ [0 0;0 1]
#      @test out.L2[ix,iλ][1,1] ≈ [0 0;0 1]'
#      @test out.L2[ix,iλ][2,1] ≈ [0 0;0 .1]'
#      @test out.L2[ix,iλ][3,1] ≈ [0 0;0 1]'

#      @test out.L2[iλ,iu][1,1] ≈ [0 0;0 -1]
#      @test out.L2[iu,iλ][1,1] ≈ [0 0;0 -1]'
     
#      @test out.L2[iλ,ia][1,1] ≈ [0,  0.2302585092994046]
#      @test out.L2[ia,iλ][1,1] ≈ [0   0.2302585092994046]

#      @test out.L2[ix,ix][1,1] ≈ [0 0;0 0]
#      @test out.L2[ix,ix][1,2] ≈ [0 -.02;0 0]
#      @test out.L2[ix,ix][1,3] ≈ [0 0;0 0]
#      @test out.L2[ix,ix][2,1] ≈ [0 0;-.02 0]
#      @test out.L2[ix,ix][2,2] ≈ [0 0;0 0]
#      @test out.L2[ix,ix][2,3] ≈ [0 0;0 0]
#      @test out.L2[ix,ix][3,1] ≈ [0 0;0 0]
#      @test out.L2[ix,ix][3,2] ≈ [0 0;0 0]
#      @test out.L2[ix,ix][3,3] ≈ [0 0;0 0]

#      @test out.L2[iu,ix][1,1] ≈ [0 0;0 0]
#      @test out.L2[iu,ix][1,2] ≈ [0 0;0 0]
#      @test out.L2[iu,ix][1,3] ≈ [0 0;0 0]
#      @test out.L2[ix,iu][1,1] ≈ [0 0;0 0]
#      @test out.L2[ix,iu][2,1] ≈ [0 0;0 0]
#      @test out.L2[ix,iu][3,1] ≈ [0 0;0 0]

#      @test out.L2[ia,ix][1,1] ≈ [-0.046051701859880924 0]
#      @test out.L2[ia,ix][1,2] ≈ [0 0]
#      @test out.L2[ia,ix][1,3] ≈ [0 0]
#      @test out.L2[ix,ia][1,1] ≈ [-0.046051701859880924,0]
#      @test out.L2[ix,ia][2,1] ≈ [0;0]
#      @test out.L2[ix,ia][3,1] ≈ [0,0]
# end


#end
;