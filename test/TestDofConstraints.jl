module TestConstraints

using Test,StaticArrays
using Muscade

Muscade.DofConstraint{:X,    Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield                       }(g,mode,gₛ,λₛ) where
                         {Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield} =
    Muscade.DofConstraint{:X,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield,typeof(g),typeof(()),typeof(mode)}(g,(),mode,gₛ,λₛ)
Muscade.DofConstraint{:U,    Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield                       }(g,mode,gₛ,λₛ) where
                         {Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield} =
    Muscade.DofConstraint{:U,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield,typeof(g),typeof(()),typeof(mode)}(g,(),mode,gₛ,λₛ)

#t,γ,dbg    = 0.,1.,(status=:testing,)
t,χ,χcv,dbg  = 0.,nothing,identity,(status=:testing,)
SP1 = (γ=1.,)
SP0 = (γ=0.,)

#---------------------------------------------------------

g(x,t)     = .3x[1] + .4x[2]
Xctc       = Muscade.variate{1,3}(SVector{3}(4,-3,10.)) # contact
Xgap       = Muscade.variate{1,3}(SVector{3}(4,3,10.))  # gap
U          = SVector{0,𝕣}()
A          = SVector{0,𝕣}()
C          = Muscade.DofConstraint{:X,    2 ,0 ,0 ,(1,1),(:t1,:t2),()   ,()    ,()   ,()    ,1    ,:λ    }
#                    DofConstraint{λclass,Nx,Nu,Na,xinod,xfield   ,uinod,ufield,ainod,afield,λinod,λfield}(g,mode,gₛ,λₛ)

@testset "X equal contact" begin
    c     = C(g,equal,1,1)
    r,χn,FB = residual(c, (Xctc,),(U,),A, t,χ,χcv,SP1,dbg)
    R,R∂X = Muscade.value_∂{1,3}(r)
    @test doflist(typeof(c)) == (inod=(1,1,1),class=(:X,:X,:X),field=(:t1,:t2,:λ))
    @test R   ≈ [-3,-4,0]
    @test R∂X ≈ [0   0   -.3; 
                 0   0   -.4; 
                 -.3 -.4 0  ]
end

@testset "X equal gap" begin
    c     = C(g,equal,1,1)
    r,χn,FB = residual(c, (Xgap,),(U,),A, t,χ,χcv,SP1,dbg)
    R,R∂X = Muscade.value_∂{1,3}(r)
    @test R   ≈ [-3, -4, -2.4]
    @test R∂X ≈ [0   0   -.3; 
                 0   0   -.4; 
                 -.3 -.4 0  ]
end

@testset "X off" begin
    c     = C(g,off,1,1)
    r,χn,FB = residual(c, (Xctc,),(U,),A, t,χ,χcv,SP1,dbg)
    R,R∂X = Muscade.value_∂{1,3}(r)
    @test R   ≈ [0,0,-10]
    @test R∂X ≈ [0 0 0; 0 0 0; 0 0 -1]
end

@testset "X positive contact" begin
    c     = C(g,positive,1,1)
    r,χn,FB = residual(c, (Xctc,),(U,),A, t,χ,χcv,SP1,dbg)
    R,R∂X = Muscade.value_∂{1,3}(r)
    @test R   ≈ [-3.0, -4.0, 0.09901951359278449]
    @test R∂X ≈ [-0.0 -0.0 -0.3; -0.0 -0.0 -0.4; -0.29708710135363803 -0.39611613513818406 -0.009709662154539889]
end

@testset "X positive gap" begin
    c     = C(g,positive,1,1)
    r,χn,FB = residual(c, (Xgap,),(U,),A, t,χ,χcv,SP1,dbg)
    R,R∂X = Muscade.value_∂{1,3}(r)
    @test R   ≈ [-3.0, -4.0, -2.2706234591223002]
    @test R∂X ≈ [-0.0 -0.0 -0.3; -0.0 -0.0 -0.4; -0.29506118058939695 -0.39341490745252927 -0.01646273136867682]
end

@testset "X positive contact γ==0" begin
    c     = C(g,positive,1,1)
    r,χn,FB = residual(c, (Xctc,),(U,),A, t,χ,χcv,SP0,dbg)
    R,R∂X = Muscade.value_∂{1,3}(r)
    @test R   ≈ [-3,-4,0]
    @test R∂X ≈ [0   0   -.3; 
                 0   0   -.4; 
                 -.3 -.4 0  ]
end

@testset "X positive gap  γ==0" begin
    c     = C(g,positive,1,1)
    r,χn,FB = residual(c, (Xgap,),(U,),A, t,χ,χcv,SP0,dbg)
    R,R∂X = Muscade.value_∂{1,3}(r)
    @test R   ≈ [-3, -4, -2.4]
    @test R∂X ≈ [0   0   -.3; 
                 0   0   -.4; 
                 -.3 -.4 0  ]
end


#---------------------------------------------------------

g(x,u,a,t) = .3u[1] + .4u[2]
Uctc       = Muscade.variate{2,3}(Muscade.variate{1,3}(SVector{3}(4,-3,10.))) # contact
Ugap       = Muscade.variate{2,3}(Muscade.variate{1,3}(SVector{3}(4, 3,10.))) # gap
Λ          = SVector{0,𝕣}()
X          = SVector{0,𝕣}()
A          = SVector{0,𝕣}()
C          = Muscade.DofConstraint{:U,    0 ,2 ,0 ,()   ,()       ,(1,1),(:t1,:t2),()   ,()    ,1    ,:λ    }
#                    DofConstraint{λclass,Nx,Nu,Na,xinod,xfield   ,uinod,ufield   ,ainod,afield,λinod,λfield}(g,mode,gₛ,λₛ)

@testset "U equal contact" begin
    c     = C(g,equal,1,1)
    r,χn,FB = lagrangian(c, Λ,(X,),(Uctc,),A, t,χ,χcv,SP1,dbg)
    R,R∂U = Muscade.value_∂{1,3}(Muscade.∂{2,3}(r))
    @test doflist(typeof(c)) == (inod=(1,1,1),class=(:U,:U,:U),field=(:t1,:t2,:λ))
    @test R   ≈ [-3,-4,0]
    @test R∂U ≈ [0   0   -.3; 
                 0   0   -.4; 
                 -.3 -.4 0  ]
end

@testset "U equal gap" begin
    c     = C(g,equal,1,1)
    r,χn,FB = lagrangian(c, Λ,(X,),(Ugap,),A, t,χ,χcv,SP1,dbg)
    R,R∂U = Muscade.value_∂{1,3}(Muscade.∂{2,3}(r))
    @test R   ≈ [-3, -4, -2.4]
    @test R∂U ≈ [0   0   -.3; 
                 0   0   -.4; 
                 -.3 -.4 0  ]
end

@testset "U off" begin
    c     = C(g,off,1,1)
    r,χn,FB = lagrangian(c, Λ,(X,),(Uctc,),A, t,χ,χcv,SP1,dbg)
    R,R∂U = Muscade.value_∂{1,3}(Muscade.∂{2,3}(r))
    @test R   ≈ [0,0,-10]
    @test R∂U ≈ [0 0 0; 0 0 0; 0 0 -1]
end

@testset "U positive contact" begin
    c     = C(g,positive,1,1)
    r,χn,FB = lagrangian(c, Λ,(X,),(Uctc,),A, t,χ,χcv,SP1,dbg)
    R,R∂U = Muscade.value_∂{1,3}(Muscade.∂{2,3}(r))
    @test R   ≈ [-3.0, -4.0, 0.09901951359278449]
    @test R∂U ≈ [-0.0 -0.0 -0.3; -0.0 -0.0 -0.4; -0.29708710135363803 -0.39611613513818406 -0.009709662154539889]
end

@testset "U positive gap" begin
    c     = C(g,positive,1,1)
    r,χn,FB = lagrangian(c, Λ,(X,),(Ugap,),A, t,χ,χcv,SP1,dbg)
    R,R∂U = Muscade.value_∂{1,3}(Muscade.∂{2,3}(r))
    @test R   ≈ [-3.0, -4.0, -2.2706234591223002]
    @test R∂U ≈ [-0.0 -0.0 -0.3; -0.0 -0.0 -0.4; -0.29506118058939695 -0.39341490745252927 -0.01646273136867682]
end

@testset "U positive contact γ==0" begin
    c     = C(g,positive,1,1)
    r,χn,FB = lagrangian(c, Λ,(X,),(Uctc,),A, t,χ,χcv,SP0,dbg)
    R,R∂U = Muscade.value_∂{1,3}(Muscade.∂{2,3}(r))
    @test R   ≈ [-3,-4,0]
    @test R∂U ≈ [0   0   -.3; 
                 0   0   -.4; 
                 -.3 -.4 0  ]
end

@testset "U positive gap  γ==0" begin
    c     = C(g,positive,1,1)
    r,χn,FB = lagrangian(c, Λ,(X,),(Ugap,),A, t,χ,χcv,SP0,dbg)
    R,R∂U = Muscade.value_∂{1,3}(Muscade.∂{2,3}(r))
    @test R   ≈ [-3, -4, -2.4]
    @test R∂U ≈ [0   0   -.3; 
                 0   0   -.4; 
                 -.3 -.4 0  ]
end

#---------------------------------------------------------


g1(x,t)       = -0.1*sin(1.2*x[1])+.2x[1]+x[2]-1 # start outside feasible domain 😃
g2(x,t)       = -.4x[1] + .3x[2]+.1
f1(x)         = -(-0.1*sin(1.2*x)+.2x-1) 
f2(x)         = (-1/.3)*(-.4x + .1)

gravity(t)      = -2.
model           = Model(:TestModel)
n1              = addnode!(model,𝕣[0,0]) 
e1              = addelement!(model,DofConstraint,[n1],xinod=(1,1),xfield=(:t1,:t2),λinod=1, λclass=:X, λfield=:λ1,gap=g1,mode=positive)
e2              = addelement!(model,DofConstraint,[n1],xinod=(1,1),xfield=(:t1,:t2),λinod=1, λclass=:X, λfield=:λ2,gap=g2,mode=positive)
e3              = addelement!(model,DofLoad   ,[n1],field=:t2,value=gravity)
initialstate    = initialize!(model)
state           = solve(StaticX;initialstate,time=[0.],verbose=false) # because there is zero physical stiffness in this model, setting γ0=0 gives singularity if one or more constraint is inactive

@testset "interior point" begin
    X = state[findlastassigned(state)].X[1][1:2]
    @test g1(X,0)   ≈ 0.0
    @test abs(g2(X,0))   < 1e-16
end
end
