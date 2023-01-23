module TestConstraints

using Test,StaticArrays
using Muscade

Constraint{    λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield                       }(g,mode,gₛ,λₛ) where
              {λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield} =
    Constraint{λclass,Nx,Nu,Na,xinod,xfield,uinod,ufield,ainod,afield,λinod,λfield,typeof(g),typeof(mode)}(g,mode,gₛ,λₛ)

t,γ,dbg    = 0.,1.,(status=:testing,)

#---------------------------------------------------------

g(x,t)     = .3x[1] + .4x[2]
Xctc       = Muscade.variate{1,3}(SVector{3}(4,-3,10.)) # contact
Xgap       = Muscade.variate{1,3}(SVector{3}(4,3,10.))  # gap
U          = SVector{0,𝕣}()
A          = SVector{0,𝕣}()
C          = Constraint{Xclass,2 ,0 ,0 ,(1,1),(:t1,:t2),()   ,()    ,()   ,()    ,1    ,:λ    }
#            Constraint{λclass,Nx,Nu,Na,xinod,xfield   ,uinod,ufield,ainod,afield,λinod,λfield}(g,mode,gₛ,λₛ)

@testset "X equal contact" begin
    c     = C(g,equal,1,1)
    R,R∂X = Muscade.value_∂{1,3}(residual(c, (Xctc,),(U,),A, t,γ,dbg))
    @test doflist(typeof(c)) == (inod=(1,1,1),class=(:X,:X,:X),field=(:t1,:t2,:λ))
    @test R   ≈ [-3,-4,0]
    @test R∂X ≈ [0   0   -.3; 
                 0   0   -.4; 
                 -.3 -.4 0  ]
end

@testset "X equal gap" begin
    c     = C(g,equal,1,1)
    R,R∂X = Muscade.value_∂{1,3}(residual(c, (Xgap,),(U,),A, t,γ,dbg))
    @test R   ≈ [-3, -4, -2.4]
    @test R∂X ≈ [0   0   -.3; 
                 0   0   -.4; 
                 -.3 -.4 0  ]
end

@testset "X off" begin
    c     = C(g,off,1,1)
    R,R∂X = Muscade.value_∂{1,3}(residual(c, (Xctc,),(U,),A, t,γ,dbg))
    @test R   ≈ [0,0,-10]
    @test R∂X ≈ [0 0 0; 0 0 0; 0 0 -1]
end

@testset "X inequal contact" begin
    c     = C(g,inequal,1,1)
    R,R∂X = Muscade.value_∂{1,3}(residual(c, (Xctc,),(U,),A, t,γ,dbg))
    @test R   ≈ [-2.9708710135363803, -3.961161351381841, 0.09901951359278449]
    @test R∂X ≈ [0.001697158861772739 0.0022628784823636484 -0.3027442975595471; 
                 0.002262878482363652 0.003017171309818198 -0.4036590634127296; 
                -0.29708710135363803 -0.39611613513818406 -0.009709662154539889]
end

@testset "X inequal gap" begin
    c     = C(g,inequal,1,1)
    R,R∂X = Muscade.value_∂{1,3}(residual(c, (Xgap,),(U,),A, t,γ,dbg))
    @test R   ≈ [-2.9506118058939697, -3.9341490745252927, -2.2706234591223002]
    @test R∂X ≈ [0.0037086134933885934 0.004944817991184798 -0.30742322556735896; 
                 0.004944817991184791 0.006593090654913064 -0.40989763408981195; 
                -0.29506118058939695 -0.39341490745252927 -0.01646273136867682]
end

@testset "X inequal contact γ==0" begin
    c     = C(g,inequal,1,1)
    R,R∂X = Muscade.value_∂{1,3}(residual(c, (Xctc,),(U,),A, t,0,dbg))
    @test R   ≈ [-3,-4,0]
    @test R∂X ≈ [0   0   -.3; 
                 0   0   -.4; 
                 -.3 -.4 0  ]
end

@testset "X inequal gap  γ==0" begin
    c     = C(g,inequal,1,1)
    R,R∂X = Muscade.value_∂{1,3}(residual(c, (Xgap,),(U,),A, t,0,dbg))
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
C          = Constraint{Uclass,0 ,2 ,0 ,()   ,()       ,(1,1),(:t1,:t2),()   ,()    ,1    ,:λ    }
#            Constraint{λclass,Nx,Nu,Na,xinod,xfield   ,uinod,ufield   ,ainod,afield,λinod,λfield}(g,mode,gₛ,λₛ)

@testset "U equal contact" begin
    c     = C(g,equal,1,1)
    R,R∂U = Muscade.value_∂{1,3}(Muscade.∂{2,3}(lagrangian(c, Λ,(X,),(Uctc,),A, t,γ,dbg)))
    @test doflist(typeof(c)) == (inod=(1,1,1),class=(:U,:U,:U),field=(:t1,:t2,:λ))
    @test R   ≈ [-3,-4,0]
    @test R∂U ≈ [0   0   -.3; 
                 0   0   -.4; 
                 -.3 -.4 0  ]
end

@testset "U equal gap" begin
    c     = C(g,equal,1,1)
    R,R∂U = Muscade.value_∂{1,3}(Muscade.∂{2,3}(lagrangian(c, Λ,(X,),(Ugap,),A, t,γ,dbg)))
    @test R   ≈ [-3, -4, -2.4]
    @test R∂U ≈ [0   0   -.3; 
                 0   0   -.4; 
                 -.3 -.4 0  ]
end

@testset "U off" begin
    c     = C(g,off,1,1)
    R,R∂U = Muscade.value_∂{1,3}(Muscade.∂{2,3}(lagrangian(c, Λ,(X,),(Uctc,),A, t,γ,dbg)))
    @test R   ≈ [0,0,-10]
    @test R∂U ≈ [0 0 0; 0 0 0; 0 0 -1]
end

@testset "U inequal contact" begin
    c     = C(g,inequal,1,1)
    R,R∂U = Muscade.value_∂{1,3}(Muscade.∂{2,3}(lagrangian(c, Λ,(X,),(Uctc,),A, t,γ,dbg)))
    @test R   ≈ [-2.9708710135363803, -3.961161351381841, 0.09901951359278449]
    @test R∂U ≈ [0.001697158861772739 0.0022628784823636484 -0.3027442975595471; 
                 0.002262878482363652 0.003017171309818198 -0.4036590634127296; 
                -0.29708710135363803 -0.39611613513818406 -0.009709662154539889]
end

@testset "U inequal gap" begin
    c     = C(g,inequal,1,1)
    R,R∂U = Muscade.value_∂{1,3}(Muscade.∂{2,3}(lagrangian(c, Λ,(X,),(Ugap,),A, t,γ,dbg)))
    @test R   ≈ [-2.9506118058939697, -3.9341490745252927, -2.2706234591223002]
    @test R∂U ≈ [0.0037086134933885934 0.004944817991184798 -0.30742322556735896; 
                 0.004944817991184791 0.006593090654913064 -0.40989763408981195; 
                -0.29506118058939695 -0.39341490745252927 -0.01646273136867682]
end

@testset "U inequal contact γ==0" begin
    c     = C(g,inequal,1,1)
    R,R∂U = Muscade.value_∂{1,3}(Muscade.∂{2,3}(lagrangian(c, Λ,(X,),(Uctc,),A, t,0.,dbg)))
    @test R   ≈ [-3,-4,0]
    @test R∂U ≈ [0   0   -.3; 
                 0   0   -.4; 
                 -.3 -.4 0  ]
end

@testset "U inequal gap  γ==0" begin
    c     = C(g,inequal,1,1)
    R,R∂U = Muscade.value_∂{1,3}(Muscade.∂{2,3}(lagrangian(c, Λ,(X,),(Ugap,),A, t,0.,dbg)))
    @test R   ≈ [-3, -4, -2.4]
    @test R∂U ≈ [0   0   -.3; 
                 0   0   -.4; 
                 -.3 -.4 0  ]
end

#---------------------------------------------------------


g1(x,t)         = -0.1*sin(1.2*x[1])+.2x[1]+x[2]-1 # start outside feasible domain 😃
g2(x,t)         = -.4x[1] + .3x[2]+.1
f1(x)         = -(-0.1*sin(1.2*x)+.2x-1) 
f2(x)         = (-1/.3)*(-.4x + .1)

gravity(t)      = -2.
model           = Model(:TestModel)
n1              = addnode!(model,𝕣[0,0]) 
e1              = addelement!(model,Constraint,[n1],xinod=(1,1),xfield=(:t1,:t2),λinod=1, λclass=:X, λfield=:λ1,g=g1,mode=inequal)
e2              = addelement!(model,Constraint,[n1],xinod=(1,1),xfield=(:t1,:t2),λinod=1, λclass=:X, λfield=:λ2,g=g2,mode=inequal)
e3              = addelement!(model,DofLoad   ,[n1],field=:t2,value=gravity)
state           = solve(staticX;model,time=[0.],verbose=false) # because there is zero physical stiffness in this model, setting γ0=0 gives singularity if one or more constraint is inactive

@testset "interior point" begin
    X = state[findlastassigned(state)].X[1][1:2]
    @test g1(X,0)   ≈ 5.294270244426968e-7
    @test g2(X,0)   ≈ 1.4402215001430019e-6
end



end
