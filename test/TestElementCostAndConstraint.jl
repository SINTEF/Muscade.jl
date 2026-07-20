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
kwargs = (;ElementType=AnchorLine, elementkwargs=(Δxₘtop=[5.,0,0], xₘbot=[250.,0], L=290., buoyancy=-5e3),
                req=@request(Fh), gap, λxinod=(1,), λxfield=(:λ,), mode,
                cost)
elcoco  = ElementCostAndConstraint(model.nod;kwargs...)
Nx,Nu,Na        = 3+1,0,2   
λ               =  SVector{Nx,𝕣}(1. for i=1:Nx)
x               = (SVector{Nx,𝕣}(1. for i=1:Nx),)
u               = (SVector{Nu,𝕣}(1. for i=1:Nu),)
a               =  SVector{Na,𝕣}(0. for i=1:Na)
_,Nz,(Λ,X,U,A)  = variate((λ,x,u,a))
dbg             = (alwaystesteverything=true,)
L,FB            = Muscade.lagrangian(elcoco, Λ,X,U,A, 0.,nothing,dbg)                 

@testset "ElementCostAndConstrants" begin
     @test Muscade.doflist(typeof(elcoco)) == (inod = (1, 1, 1, 1, 2, 2), class = (:X, :X, :X, :X, :A, :A), field = (:λ, :tx1, :tx2, :rx3,  :ΔL, :Δbuoyancy))
     @test value{1}(L) ≈ 1.441928261291504e6
     @test ∂{1,Nz}(L) ≈ [       0.0, -438861.13074456755,    9278.602091074139,       1.871510789932893e6,      -1.926851845351649e11,  -86753.51731872559,    8308.805033385754,       1.5339992031555176e6,      -2.5203831430664062e7,       1.441928261291504e6]
end

req            = @request λ,eleres(cr)
L,FB,eleres    = Muscade.lagrangian(elcoco, λ, x,u,a, 0.,nothing,dbg,req)                 

@testset "ElementCostAndConstraintResult" begin
     @test eleres.λ         ≈ [1.]
     @test eleres.eleres.cr ≈ 87.79184120068672
     @test eleres.eleres.Fh ≈ 438959.2060034336
end

## Test draw

nel       = 1
nXder     = nUder = 1
EL        = [elcoco]
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

## Test addin!(AssemblyDirect,...)

addelement!(model,ElementCostAndConstraint,[n1,n2];kwargs...)
OX,OU,IA              = 0,0,0
wanted                = Muscade.Wanted{1,OX+1,OU+1,IA}(:all,:all)
initialstate          = initialize!(model)
model,dis             = initialstate.model, initialstate.dis
out,asm,dofgr         = Muscade.prepare(Muscade.AssemblyDirect,model,dis,wanted)          # mem and assembler for system at any given step
Muscade.zero!(out)
ieletyp,iele          = 1,1
asmᵢ                  = asm[:,ieletyp]
scale                 = dis.dis[ieletyp].scale
t,Δt                  = 0.,1.
SP                    = nothing
Muscade.addin!{:matrices}(out,asmᵢ,iele,scale,elcoco,Val(true),(λ,),x,u,a,t,Δt,SP,dbg)


;#end
