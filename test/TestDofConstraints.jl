module TestConstraints

using Test,StaticArrays
using Muscade,StaticArrays

Muscade.DofConstraint{:X,Nλ,Nx,Nu,Na,λinod,λfield,xinod,xfield,uinod,ufield,ainod,afield                       }(g,mode) where
                        {Nλ,Nx,Nu,Na,λinod,λfield,xinod,xfield,uinod,ufield,ainod,afield} =
Muscade.DofConstraint{:X,Nλ,Nx,Nu,Na,λinod,λfield,xinod,xfield,uinod,ufield,ainod,afield,typeof(g),typeof(()),typeof(mode)}(g,(),mode)

Muscade.DofConstraint{:U,Nλ,Nx,Nu,Na,λinod,λfield,xinod,xfield,uinod,ufield,ainod,afield                       }(g,mode) where
                        {Nλ,Nx,Nu,Na,λinod,λfield,xinod,xfield,uinod,ufield,ainod,afield} =
Muscade.DofConstraint{:U,Nλ,Nx,Nu,Na,λinod,λfield,xinod,xfield,uinod,ufield,ainod,afield,typeof(g),typeof(()),typeof(mode)}(g,(),mode)

t,dbg  = 0.,(status=:testing,)
SP1 = (γ=1.,)
SP0 = (γ=0.,)

#---------------------------------------------------------

@functor with() g(x,t)     = SVector(.3x[1] + .4x[2])
Xctc       = Muscade.variate{1,3}(SVector{3}(10.,4,-3)) # contact
Xgap       = Muscade.variate{1,3}(SVector{3}(10.,4,3))  # gap
U          = SVector{0,𝕣}()
A          = SVector{0,𝕣}()
C          = Muscade.DofConstraint{:X,    1 ,2 ,0 ,0 ,(1,) ,(:λ,) ,(1,1),(:t1,:t2),()   ,()    ,()   ,()    }
#                    DofConstraint{λclass,Nλ,Nx,Nu,Na,λinod,λfield,xinod,xfield   ,uinod,ufield,ainod,afield}(g,mode)

@testset "X equal contact" begin
    c     = C(g,equal)
    r,FB  = Muscade.residual(c, (Xctc,),(U,),A, t,SP1,dbg)
    R,R∂X = Muscade.value_∂{1,3}(r)
    @test Muscade.doflist(typeof(c)) == (inod = (1, 1, 1), class = (:X, :X, :X), field = (:λ, :t1, :t2))
    @test R   ≈ [0,-3,-4]
    @test R∂X ≈  [-0.0 -0.3 -0.4; -0.3 -0.0 -0.0; -0.4 -0.0 -0.0]
end

@testset "X equal gap" begin
    c     = C(g,equal)
    r,FB  = Muscade.residual(c, (Xgap,),(U,),A, t,SP1,dbg)
    R,R∂X = Muscade.value_∂{1,3}(r)
    @test R   ≈ [-2.4,-3, -4 ]
    @test R∂X ≈ [-0.0 -0.3 -0.4; -0.3 -0.0 -0.0; -0.4 -0.0 -0.0]
end

@testset "X off" begin
    c     = C(g,off)
    r,FB  = Muscade.residual(c, (Xctc,),(U,),A, t,SP1,dbg)
    R,R∂X = Muscade.value_∂{1,3}(r)
    @test R   ≈ [-10,0,0]
    @test R∂X ≈ [-1 0 0; 0 0 0; 0 0 0]
end

@testset "X positive contact" begin
    c     = C(g,positive)
    r,FB  = Muscade.residual(c, (Xctc,),(U,),A, t,SP1,dbg)
    R,R∂X = Muscade.value_∂{1,3}(r)
    @test R   ≈ [1.,-3.0, -4.0] 
    @test R∂X ≈ [0.0 -3.0 -4.0; -0.3 -0.0 -0.0; -0.4 -0.0 -0.0]
end

@testset "X positive gap" begin
    c     = C(g,positive)
    r,FB  = Muscade.residual(c, (Xgap,),(U,),A, t,SP1,dbg)
    R,R∂X = Muscade.value_∂{1,3}(r)
    @test R   ≈ [-23.,-3.0, -4.0]
    @test R∂X ≈ [-2.4 -3.0 -4.0; -0.3 -0.0 -0.0; -0.4 -0.0 -0.0] 
end

@testset "X positive contact γ==0" begin
    c     = C(g,positive)
    r,FB = Muscade.residual(c, (Xctc,),(U,),A, t,SP0,dbg)
    R,R∂X = Muscade.value_∂{1,3}(r)
    @test R   ≈ [0,-3,-4]
    @test R∂X ≈ [0.0 -3.0 -4.0; -0.3 -0.0 -0.0; -0.4 -0.0 -0.0]
end

@testset "X positive gap  γ==0" begin
    c     = C(g,positive)
    r,FB  = Muscade.residual(c, (Xgap,),(U,),A, t,SP0,dbg)
    R,R∂X = Muscade.value_∂{1,3}(r)
    @test R   ≈ [-24.,-3.0, -4.0]
    @test R∂X ≈  [-2.4 -3.0 -4.0; -0.3 -0.0 -0.0; -0.4 -0.0 -0.0]
end


#---------------------------------------------------------

@functor with() g(x,u,a,t) = SVector(.3u[1] + .4u[2])
Uctc       = Muscade.variate{2,3}(Muscade.variate{1,3}(SVector{3}(10.,4,-3))) # contact
Ugap       = Muscade.variate{2,3}(Muscade.variate{1,3}(SVector{3}(10.,4, 3))) # gap
Λ          = SVector{0,𝕣}()
X          = SVector{0,𝕣}()
A          = SVector{0,𝕣}()
C          = Muscade.DofConstraint{:U,    1, 0 ,2 ,0 ,(1,) ,(:λ,) ,()   ,()       ,(1,1),(:t1,:t2),()   ,()    }
#                    DofConstraint{λclass,Nλ,Nx,Nu,Na,λinod,λfield,xinod,xfield   ,uinod,ufield   ,ainod,afield}(g,mode)

@testset "U equal contact" begin
    c     = C(g,equal)
    r,FB  = Muscade.lagrangian(c, Λ,(X,),(Uctc,),A, t,SP1,dbg)
    R,R∂U = Muscade.value_∂{1,3}(Muscade.∂{2,3}(r))
    @test Muscade.doflist(typeof(c)) == (inod = (1, 1, 1), class = (:U, :U, :U), field = (:λ, :t1, :t2))
    @test R   ≈ [0,-3,-4]
    @test R∂U ≈ [0.0 -0.3 -0.4; -0.3 0.0 0.0; -0.4 0.0 0.0]
end

@testset "U equal gap" begin
    c     = C(g,equal)
    r,FB  = Muscade.lagrangian(c, Λ,(X,),(Ugap,),A, t,SP1,dbg)
    R,R∂U = Muscade.value_∂{1,3}(Muscade.∂{2,3}(r))
    @test R   ≈ [-2.4,-3, -4]
    @test R∂U ≈ [0.0 -0.3 -0.4; -0.3 0.0 0.0; -0.4 0.0 0.0] 
end

@testset "U off" begin
    c     = C(g,off)
    r,FB  = Muscade.lagrangian(c, Λ,(X,),(Uctc,),A, t,SP1,dbg)
    R,R∂U = Muscade.value_∂{1,3}(Muscade.∂{2,3}(r))
    @test R   ≈ [-10,0,0]
    @test R∂U ≈ [-1 0 0; 0 0 0; 0 0 0]
end

@testset "U positive contact" begin
    c     = C(g,positive)
    r,FB  = Muscade.lagrangian(c, Λ,(X,),(Uctc,),A, t,SP1,dbg)
    @show typeof(r)
    R,R∂U = Muscade.value_∂{1,3}(Muscade.∂{2,3}(r))
    @show size(R)
    @test R   ≈ [1.,-3.0, -4.0]
    @test R∂U ≈ [0.0 -3.0 -4.0; -0.3 -0.0 -0.0; -0.4 -0.0 -0.0]
end

@testset "U positive gap" begin
    c     = C(g,positive)
    r,FB  = Muscade.lagrangian(c, Λ,(X,),(Ugap,),A, t,SP1,dbg)
    R,R∂U = Muscade.value_∂{1,3}(Muscade.∂{2,3}(r))
    @test R   ≈ [-23.,-3.0, -4.0 ]
    @test R∂U ≈ [-2.4 -3.0 -4.0; -0.3 -0.0 -0.0; -0.4 -0.0 -0.0]
end

@testset "U positive contact γ==0" begin
    c     = C(g,positive)
    r,FB  = Muscade.lagrangian(c, Λ,(X,),(Uctc,),A, t,SP0,dbg)
    R,R∂U = Muscade.value_∂{1,3}(Muscade.∂{2,3}(r))
    @test R   ≈ [0,-3,-4]
    @test R∂U ≈ [0 -3.0 -4.0; -0.3 -0.0 -0.0; -0.4 -0.0 -0.0]
end

@testset "U positive gap  γ==0" begin
    c     = C(g,positive)
    r,FB  = Muscade.lagrangian(c, Λ,(X,),(Ugap,),A, t,SP0,dbg)
    R,R∂U = Muscade.value_∂{1,3}(Muscade.∂{2,3}(r))
    @test R   ≈ [-24.,-3.0, -4.0]
    @test R∂U ≈  [-2.4 -3.0 -4.0; -0.3 -0.0 -0.0; -0.4 -0.0 -0.0]
end

#---------------------------------------------------------


@functor with() g1(x,t)       = SVector(-0.1*sin(1.2*x[1])+.2x[1]+x[2]+1 )
@functor with() g2(x,t)       = SVector(-.4x[1] + .3x[2]+.1)
@functor with() f1(x)         = -(-0.1*sin(1.2*x)+.2x-1) 
@functor with() f2(x)         = (-1/.3)*(-.4x + .1)
@functor with() gravity(t)    = -2.

model           = Model(:TestModel)
n1              = addnode!(model,𝕣[0,0]) 
e1              = addelement!(model,DofConstraint,[n1],λclass=:X, λinod=(1,),λfield=(:λ1,), xinod=(1,1),xfield=(:t1,:t2), gap=g1,mode=positive)
e2              = addelement!(model,DofConstraint,[n1],λclass=:X, λinod=(1,),λfield=(:λ2,), xinod=(1,1),xfield=(:t1,:t2), gap=g2,mode=positive)
e3              = addelement!(model,DofLoad      ,[n1],field=:t2,value=gravity)
initialstate    = initialize!(model)
initialstate    = setdof!(initialstate,1.;field=:λ1)
initialstate    = setdof!(initialstate,1.;field=:λ2)
state           = solve(SweepX{0};initialstate,time=[0.],verbose=false,silenterror=false) # because there is zero physical stiffness in this model, setting γ0=0 gives singularity if one or more constraint is inactive

@testset "interior point" begin
    X = state[findlastassigned(state)].X[1][2:3]
    @test abs(g1(X,0)[1])   < 1e-6
    @test abs(g2(X,0)[1])   < 1e-5
end
end
