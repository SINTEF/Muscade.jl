#module TestDrawBeamElement
using Muscade, StaticArrays, Test
using Printf
using Muscade: lines!,scatter!,mesh!

include("../examples/BeamElements.jl")
include("../examples/BeamElementsDraw.jl")

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

# using GLMakie
# fig      = Figure(size = (500,500))
# display(fig) # open interactive window (gets closed down by "save")
# axe      = Axis3(fig[1,1],title="Test",xlabel="X",ylabel="Y",zlabel="Z",aspect=:data,viewmode=:fit,perspectiveness=.5)

α = 2π*(0:19)/20
circle = 0.05*[cos.(α) sin.(α)]'
square = 0.1*[1 -1 -1 1;1 1 -1 -1]

axe = Muscade.SpyAxis()
draw!(axe,state;EulerBeam3D=(;style=:simple))
draw!(axe,state;EulerBeam3D=(;style=:shape,nseg=10,frame=true,Uscale=0.1))
draw!(axe,state;EulerBeam3D=(;style=:solid,nseg=10,section = circle,marking=true,Uscale=0.1))
@testset "drawing" begin
    @test axe.call[1].fun == :scatter!
    @test axe.call[1].args[1][][:,1:2] ≈ [  0.0150628   0.336131;
                                          0.0118553  -0.0150178;
                                         -0.000945   -0.00591755 ] rtol=1e-4
    @test axe.call[2].fun == :lines!
    @test axe.call[2].args[1][][:,1:2] ≈ [  0.0150628   0.336131;
                                          0.0118553  -0.0150178;
                                         -0.000945   -0.00591755 ] rtol=1e-4
    @test axe.call[4].fun == :lines!
    @test axe.call[4].args[1][][:,1:5] ≈ [  0.0150628   0.0466323   0.0787699    0.111307     0.144076;
                                          0.0118553   0.0121641   0.0118317    0.0108489    0.00920646;
                                         -0.000945   -0.0050873  -0.00640313  -0.00564867  -0.00358008 ] rtol=1e-4
    @test axe.call[1].fun == :scatter!
    @test axe.call[1].args[1][][:,1:2] ≈ [  0.0150628   0.336131;
                                          0.0118553  -0.0150178;
                                         -0.000945   -0.00591755 ] rtol=1e-4
    @test axe.call[5].fun == :lines!
    @test axe.call[5].args[1][][:,1:2] ≈ [  0.175597     0.283664 ;
                                         -0.00158125  -0.0113104 ;
                                         -0.00343127  -0.0273596 ] rtol=1e-4
    @test  axe.call[6].fun == :lines!
    @test  axe.call[6].args[1][][:,1:5] ≈ [ 0.0134968  0.0355003  0.359701    0.337697   0.0134968;
                                          0.0130125  0.0775831  0.0483956  -0.016175   0.0130125;
                                          0.0324613  0.105581   0.0337958  -0.0393238  0.0324613] rtol=1e-4
    @test axe.call[7].fun == :scatter!
    @test axe.call[7].args[1][][:,1:2] ≈ [  0.0150628   0.336131;
                                          0.0118553  -0.0150178;
                                         -0.000945   -0.00591755 ] rtol=1e-4
    @test axe.call[8].fun == :mesh!
    @test axe.call[8].args[1][][:,1:10] ≈ [ 0.010907    0.0134416   0.016135    0.0187234  0.0209534  0.0226069  0.0235219  0.0236088    0.0228593   0.0213465;
                                          0.0568619   0.0613883   0.061066    0.0559267  0.0464733  0.0336313  0.0186577  0.00301825  -0.0117562  -0.0242193;
                                         -0.0223258  -0.00756755  0.00783893  0.0223856  0.0346485  0.0434272  0.0478625  0.0475202    0.0424337   0.0331011] rtol = 1e-4
    @test axe.call[8].args[2][][1:10,:]' == [  1   1   2   2   3   3   4   4   5   5;  2  22   3  23   4  24   5  25   6  26; 22  21  23  22  24  23  25  24  26  25]
end

#end# module

