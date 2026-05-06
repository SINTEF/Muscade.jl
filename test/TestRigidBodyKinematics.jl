module TestRigidBodyKinematics

using Test, Muscade, Muscade.Toolbox, StaticArrays, LinearAlgebra

 model          = Model(:TestModel)
 master_node    = addnode!(model, 𝕣[0., 0., 0.])
 emissary_node  = addnode!(model, 𝕣[1/√2, 1/√2, 0.])
 ERCelement        = ExcentricRigidConnection(model.nod)

@testset "constructor" begin
    dofs  = Muscade.doflist(typeof(ERCelement))
    @test dofs.inod == (2, 2, 1, 1, 1, 2, 2) 
    @test dofs.class == (:X, :X, :X, :X, :X, :X, :X) 
    @test dofs.field == (:λt1, :λt2, :t1, :t2, :r3, :t1, :t2)
    @test ERCelement.gargs[1] ≈ π/4
    @test ERCelement.gargs[2] ≈ 1.0
end

U          = SVector{0,𝕣}()
A          = SVector{0,𝕣}()
t,SP,dbg  = 0.,(;),(status=:testing,)
δt1m, δt2m, δr3m  = 0., 1., π/4; 
δt1e, δt2e        = -1/√2, 2-1/√2; 
x0 =  SVector(0.,0.,δt1m, δt2m, δr3m, δt1e, δt2e)
R,FB  = Muscade.residual(ERCelement, (x0,),(U,),A, t,SP,dbg)
@testset "residual" begin
    @test R ≈ [0., 0., 0., 0., 0., 0., 0.] atol=10*eps()
end


zsv7 = SVector(0.,0.,  0.,0.,0.,   0.,0.)
X = ntuple(i->zsv7, 7)
U = (SVector{0,𝕣}(),)
A = SVector{0,Float64}()
out  = Muscade.diffed_residual(ERCelement; X,U,A)
K = out.∇R[2][1]
# Muscade.print_element_array(ERCelement,:X,K) 
@testset "diffed residual" begin
    @test K[1,3] ≈ 1.
    @test K[1,4] ≈ 0.
    @test K[1,5] ≈ -√2/2
    @test K[1,6] ≈ -1.
    @test K[1,7] ≈ 0.

    @test K[2,3] ≈ 0.
    @test K[2,4] ≈ 1.
    @test K[2,5] ≈ √2/2
    @test K[2,6] ≈ 0.
    @test K[2,7] ≈ -1.

    @test norm(K[1:2,1:2])       ≈ 0.
    @test norm(K[3:7,3:7])       ≈ 0.
end


end