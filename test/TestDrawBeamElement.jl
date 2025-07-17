module TestDrawBeamElement
using Muscade, StaticArrays, Test
using Printf

include("../examples/BeamElement.jl")

# Beam simply supported at both ends  
L   = 1;  # Beam length [m]
q   = 0.0;  # Uniform lateral load [N/m]
EI₂ = 1;  # Bending stiffness [Nm²]
EI₃ = 1;  # Bending stiffness [Nm²]
EA  = 1e6;  # Axial stiffness [N]
GJ  = 1e6;  # Torsional stiffness [Nm²]
μ   = 1;
ι₁  = 1;

# Create model
nel         = 3
nodeCoord   = hcat( vcat(0:L/nel:L,zeros(nel)),zeros(Float64,nel+1 + nel,2))
mat         = BeamCrossSection(EA=EA,EI₂=EI₂,EI₃=EI₃,GJ=GJ,μ=μ,ι₁=ι₁)
model       = Model(:TestModel)
nodid       = addnode!(model,nodeCoord)
mesh        = hcat(nodid[1:nel],nodid[2:nel+1],nodid[nel+2:2nel+1])
eleid       = addelement!(model,EulerBeam3D{true},mesh;mat=mat,orient2=SVector(0.,1.,0.))

# Static analysis
state    = initialize!(model);
setdof!(state,  [0.0150628,  0.00279812,  0.0323863,   0.00527616 ],nodID=nodid[1:nel+1     ]         ,field=:t1)
setdof!(state,  [0.0118553, -0.0150178 ,  0.0118042,   0.000708444],nodID=nodid[1:nel+1     ]         ,field=:t2)
setdof!(state,  [-0.000945, -0.00591755,  0.0133323,  -0.0545434  ],nodID=nodid[1:nel+1     ]         ,field=:t3)
setdof!(state,3*[-0.149568, -0.32468   , -0.385368 ,  -0.369911   ],nodID=nodid[1:nel+1     ]         ,field=:r1)
setdof!(state,  [0.167858 , 0.310638   ,  0.15357  ,   0.162442   ],nodID=nodid[1:nel+1     ]         ,field=:r2)
setdof!(state,  [0.0480129,  -0.0636363, -0.596994 ,  -0.149229   ],nodID=nodid[1:nel+1     ]         ,field=:r3)
setdof!(state,  [0.       ,  0.2        , 0.4                     ],nodID=nodid[nel+2:2nel+1],class=:U,field=:t1)
setdof!(state,  [0.       ,  0.         , 0.                      ],nodID=nodid[nel+2:2nel+1],class=:U,field=:t2)
setdof!(state,  [1.       ,  1.3        , 1.6                     ],nodID=nodid[nel+2:2nel+1],class=:U,field=:t3)

state2 = copy(state)
state2.X[1] .+= 1.

# using GLMakie
# fig      = Figure(size = (500,500))
# display(fig) # open interactive window (gets closed down by "save")
# axis      = Axis3(fig[1,1],title="Test",xlabel="X",ylabel="Y",zlabel="Z",aspect=:data,viewmode=:fit,perspectiveness=.5)

α = 2π*(0:19)/20
circle = 0.05*[cos.(α) sin.(α)]'
square = 0.1*[1 -1 -1 1;1 1 -1 -1]


# f(i) = axis.call[i].fun, axis.call[i].args[1][][:,1:2]
axis     = Muscade.SpyAxis()
graphic = draw!(axis,state;EulerBeam3D=(;style=:shape,nseg=10,frame=true,Uscale=0.1))
@testset "drawing shape" begin
    @test axis.call[1].fun == :scatter!
    @test axis.call[1].args[1][][:,1:2] ≈ [  0.0150628   0.336131;
                                          0.0118553  -0.0150178;
                                         -0.000945   -0.00591755 ] rtol=1e-4
    @test axis.call[2].fun == :lines!
    @test axis.call[2].args[1][][:,1:2] ≈ [  0.0150628   0.0466323;
                                          0.0118553   0.0121641;
                                         -0.000945   -0.0050873 ] rtol=1e-4
    @test axis.call[3].fun == :lines!
    @test axis.call[3].args[1][][:,1:2] ≈ [0.0150628 0.037066306833087495; 
                                            0.0118553 0.07642591077859431; 
                                            -0.000945 0.07217464107355528] rtol=1e-4
end
draw!(graphic,state2;EulerBeam3D=(;style=:shape))
@testset "drawupdate shape" begin
    @test axis.call[1].args[1][][:,1:2].-1 ≈ [  0.0150628   0.336131;
                                              0.0118553  -0.0150178;
                                             -0.000945   -0.00591755 ] rtol=1e-4
    @test axis.call[2].args[1][][:,1:2] ≈ [1.01506   1.02511;
                                          1.01186   1.0142;
                                          0.999055  1.00503 ] rtol=1e-4
    @test axis.call[3].args[1][][:,1:2] ≈ [1.0150628 1.1033134483690596; 
                                            1.0118553 1.0437310556703177; 
                                            0.999055 1.0336347522097369] rtol=1e-4
end



axis     = Muscade.SpyAxis()
graphic = draw!(axis,state;EulerBeam3D=(;style=:solid,nseg=10,section = circle,marking=true,Uscale=0.1))
@testset "drawing solid" begin
    @test axis.call[1].fun == :scatter!
    @test axis.call[1].args[1][][:,1:2] ≈ [  0.0150628   0.336131;
                                          0.0118553  -0.0150178;
                                         -0.000945   -0.00591755 ] rtol=1e-4
    @test axis.call[2].fun == :mesh!
    @test axis.call[2].args[1][][:,1:10] ≈ [ 0.010907    0.0134416   0.016135    0.0187234  0.0209534  0.0226069  0.0235219  0.0236088    0.0228593   0.0213465;
                                          0.0568619   0.0613883   0.061066    0.0559267  0.0464733  0.0336313  0.0186577  0.00301825  -0.0117562  -0.0242193;
                                         -0.0223258  -0.00756755  0.00783893  0.0223856  0.0346485  0.0434272  0.0478625  0.0475202    0.0424337   0.0331011] rtol = 1e-4
    @test axis.call[2].args[2][][1:10,:]' == [  1   1   2   2   3   3   4   4   5   5;  2  22   3  23   4  24   5  25   6  26; 22  21  23  22  24  23  25  24  26  25]
end
draw!(graphic,state2;EulerBeam3D=(;style=:solid))
@testset "drawing solid" begin
    @test axis.call[1].args[1][][:,1:2] ≈ [ 1.01506   1.33613 ;
                                           1.01186   0.984982;
                                           0.999055  0.994082] rtol=1e-4
    @test axis.call[2].args[1][][:,1:10] ≈ [  0.996402  1.01164  1.02722  1.04161  1.05339  1.06143  1.06493  1.06354   1.05741   1.04714 ;
                                             1.03417   1.03545  1.03441  1.03117  1.02604  1.01952  1.01225  1.00494   0.998305  0.993   ;
                                             1.03972   1.04301  1.04199  1.03677  1.02786  1.01613  1.00272  0.988963  0.976189  0.965654] rtol = 1e-4
end


end# module

