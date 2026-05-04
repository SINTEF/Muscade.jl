module TestRigidBodyKinematics

using Test, Muscade, Muscade.Toolbox, StaticArrays

 model          = Model(:TestModel)
 master_node    = addnode!(model, 𝕣[0., 0., 0.])
 emissary_node  = addnode!(model, 𝕣[1/√2, 1/√2, 0.])
 ERCelement        = ExcentricRigidConnection(model.nod)

@testset "ExcentricRigidConnection constructor" begin
dofs  = Muscade.doflist(typeof(ERCelement))
@test dofs.inod == (2, 2, 1, 1, 1, 2, 2) 
@test dofs.class == (:X, :X, :X, :X, :X, :X, :X) 
@test dofs.field == (:λt1, :λt2, :t1, :t2, :r3, :t1, :t2)
@test ERCelement.gargs[1] ≈ π/4
@test ERCelement.gargs[2] ≈ 1.0
end

Λ =  SVector(0.,0.)
U          = SVector{0,𝕣}()
A          = SVector{0,𝕣}()
t,dbg  = 0.,(status=:testing,)
δt1m, δt2m, δr3m  = 0., 1., π/4; 
δt1e, δt2e        = -1/√2, 2-1/√2; 
X =  SVector(0.,0.,δt1m, δt2m, δr3m, δt1e, δt2e)
R,FB  = Muscade.residual(ERCelement, (X,),(U,),A, t,nothing,dbg)
@testset "ExcentricRigidConnection residual" begin
@test R ≈ [0., 0., 0., 0., 0., 0., 0.] atol=10*eps()
end

end
# Xctc       = Muscade.variate{1,3}(SVector{3}(0.,0.,0.)) # contact
# L,FB = Muscade.lagrangian(ERCelement,Λ,(X,),(U,),A,0.0,nothing,(testall=true,))