module TestElementCost

using Test,StaticArrays
using Muscade
include("SomeElements.jl")

model           = Model(:TestModel)
n1              = addnode!(model,𝕣[0,0,+100]) # turbine
n3              = addnode!(model,𝕣[])  # Anod for anchor

@once cost(eleres,X,U,A,t) = eleres.Fh^2
el = ElementCost(model.nod;req=@request(Fh),cost,ElementType=AnchorLine, 
                 elementkwargs=(Δxₘtop=[5.,0,0], xₘbot=[250.,0], L=290., buoyancy=-5e3))
d  = Muscade.doflist(typeof(el))
Nx,Nu,Na        = 3,0,2   
Nz              = 2Nx+Nu+Na     
iλ,ix,iu,ia     = Muscade.gradientpartition(Nx,Nx,Nu,Na) 
ΔZ              = δ{1,Nz,𝕣}()                 
ΔΛ,ΔX,ΔU,ΔA     = view(ΔZ,iλ),view(ΔZ,ix),view(ΔZ,iu),view(ΔZ,ia) 
# Λ =  SVector(0.,0.,0.)#zeros(Nx)
# X = (SVector(1.,1.,1.),)#(ones(Nx),)
# U = (SVector{0,𝕣}(),)#(zeros(Nu),)
# A =  SVector(0.,0.)# zeros(Na)
Λ =  SVector{Nx}(0. for i=1:Nx)
X = (SVector{Nx}(1. for i=1:Nx),)
U = (SVector{Nu,𝕣}(0. for i=1:Nu),)
A =  SVector{Na}(0. for i=1:Na)


L,FB  = Muscade.lagrangian(el, Λ+ΔΛ, (∂0(X)+ΔX,),(∂0(U)+ΔU,),A+ΔA, 0.,nothing,(testall=true,))                 

@testset "ElementCost" begin
     @test d == (inod = (1, 1, 1, 2, 2), class = (:X, :X, :X, :A, :A), field = (:tx1, :tx2, :rx3, :ΔL, :Δbuoyancy))
     @test value{1}(L) ≈ 1.926851845351649e11
     @test ∂{1,Nz}(L) ≈ [-438861.1307445675,9278.602091074139,1.8715107899328927e6,-2.322235123921358e10,4.9097753633879846e8,9.903105530914653e10,-6.735986859485705e12,3.853703690703298e11]
end

###

@once gap(eleres,X,U,A,t) = eleres.Fh^2
el = ElementConstraint(model.nod;req=@request(Fh),gap,ElementType=AnchorLine,λinod=1,λfield=:λ,mode=equal, 
                 elementkwargs=(Δxₘtop=[5.,0,0], xₘbot=[250.,0], L=290., buoyancy=-5e3))

d               = Muscade.doflist(typeof(el))
Nx,Nu,Na        = 3,0+1,2   
Nz              = 2Nx+Nu+Na     
iλ,ix,iu,ia     = Muscade.gradientpartition(Nx,Nx,Nu,Na) 
ΔZ              = δ{1,Nz,𝕣}()                 
ΔΛ,ΔX,ΔU,ΔA     = view(ΔZ,iλ),view(ΔZ,ix),view(ΔZ,iu),view(ΔZ,ia) 
Λ =  SVector{Nx}(0. for i=1:Nx)
X = (SVector{Nx}(1. for i=1:Nx),)
U = (SVector{Nu}(1. for i=1:Nu),)
A =  SVector{Na}(0. for i=1:Na)

L,FB  = Muscade.lagrangian(el, Λ+ΔΛ, (∂0(X)+ΔX,),(∂0(U)+ΔU,),A+ΔA, 0.,nothing,(testall=true,))                 

@testset "ElementConstraint" begin
     @test d == (inod = (1, 1, 1, 2, 2, 1), class = (:X, :X, :X, :A, :A, :U), field = (:tx1, :tx2, :rx3,  :ΔL, :Δbuoyancy, :λ))
     @test value{1}(L) ≈ -1.926851845351649e11
     @test ∂{1,Nz}(L) ≈ [-438861.1307445675, 9278.602091074139, 1.8715107899328927e6, 2.322235123921358e10, -4.9097753633879846e8, -9.903105530914653e10, -1.926851845351649e11, 6.735986859485705e12, -3.853703690703298e11]
end

req = @request λ,eleres(cr)

L,FB,eleres  = Muscade.lagrangian(el, Λ+ΔΛ, X,U,A, 0.,nothing,(testall=true,),req)                 
@testset "ElementConstraintResult" begin
     @test eleres.λ ≈ 1.
     @test eleres.eleres.cr ≈ 87.79184120068672
     @test eleres.eleres.Fh ≈ 438959.2060034336
end

end
