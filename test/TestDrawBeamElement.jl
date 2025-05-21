module TestDrawBeamElement
using Muscade, StaticArrays, Test
using Printf

include("../examples/BeamElements.jl")

# Beam simply supported at both ends  
L   = 1;  # Beam length [m]
q   = 0.0;  # Uniform lateral load [N/m]
EI₂ = 1;  # Bending stiffness [Nm²]
EI₃ = 1;  # Bending stiffness [Nm²]
EA  = 1e6;  # Axial stiffness [N]
GJ  = 1e6;  # Torsional stiffness [Nm²]
μ   = 1;
ι₁  = 1;

# Create model.03*randn(Nnod)
nel         = 3
Nnod        = nel+1   
nodeCoord   = hcat((0:L/nel:L),zeros(Float64,Nnod,2))
mat         = BeamCrossSection(EA=EA,EI₂=EI₂,EI₃=EI₃,GJ=GJ,μ=μ,ι₁=ι₁)
model       = Model(:TestModel)
nodid       = addnode!(model,nodeCoord)
mesh        = hcat(nodid[1:Nnod-1],nodid[2:Nnod])
eleid       = addelement!(model,EulerBeam3D,mesh;mat=mat,orient2=SVector(0.,1.,0.))
[addelement!(model,Hold,[nodid[1]]  ;field) for field∈[:t1,:t2,:t3,:r1]]                                # Simply supported end 1
[addelement!(model,Hold,[nodid[end]];field) for field∈[:t1,:t2,:t3,:r1]]                                # Simply supported end 2
[addelement!(model,Hold,[nodid[nodeidx]];field=:t3) for nodeidx∈2:Nnod-1]                               # Enforce beam motions in one dimension to obtain planar modeshapes
[addelement!(model,DofLoad,[nodid[nodeidx]];field=:t2,value=t->sin(t)*q*L/Nnod) for nodeidx=1:Nnod]; # Distributed vertical load q

# Static analysis
state    = initialize!(model);
setdof!(state,[0.0150628,  0.00279812,  0.0323863,   0.00527616],nodID=nodid,field=:t1)
setdof!(state,[0.0118553, -0.0150178 ,  0.0118042,   0.000708444],nodID=nodid,field=:t2)
setdof!(state,[-0.000945, -0.00591755,  0.0133323,  -0.0545434],nodID=nodid,field=:t3)
setdof!(state,3*[-0.149568, -0.32468   , -0.385368 ,  -0.369911],nodID=nodid,field=:r1)
setdof!(state,[0.167858 , 0.310638   ,  0.15357  ,   0.162442],nodID=nodid,field=:r2)
setdof!(state,[0.0480129,  -0.0636363, -0.596994 ,  -0.149229],nodID=nodid,field=:r3)

# using GLMakie
# fig      = Figure(size = (500,500))
# display(fig) # open interactive window (gets closed down by "save")
# axe      = Axis3(fig[1,1],title="Test",xlabel="X",ylabel="Y",zlabel="Z",aspect=:data,viewmode=:fit,perspectiveness=.5)

α = 2π*(0:19)/20
circle = 0.05*[cos.(α) sin.(α)]'
square = 0.1*[1 -1 -1 1;1 1 -1 -1]


include("GLMakieTester.jl")
axe = SpyAxe()
draw(axe,state;EulerBeam3D=(;style=:simple))
draw(axe,state;EulerBeam3D=(;style=:shape,nseg=10,frame=true))
draw(axe,state;EulerBeam3D=(;style=:solid,nseg=10,section = circle,marking=true))
@testset "drawing" begin
    @test  axe.call[1].fun == :lines!
    @test  isequal(axe.call[1].args[1][1,:], [   0.0150628,   0.33613145333333333, NaN,   0.3361314533333334,   0.6990529666666666, NaN,   0.6990529666666665,   1.00527616, NaN])
    @test axe.call[2].fun == :scatter!
    @test  isequal(axe.call[2].args[1],axe.call[1].args[1])  
    @test axe.call[3].fun == :lines!
    @test axe.call[3].args[1][:,1:10] ≈ [  0.0150628   0.0464197    0.078392     0.110811     0.14351      0.176318     0.209069     0.241595     0.273725     0.305294;
                                           0.0118553   0.0112004    0.0101185    0.0086003    0.00663665   0.00421828   0.00133595  -0.00201958  -0.00585756  -0.0101872;
                                           -0.000945   -0.00565553  -0.00741332  -0.00697455  -0.00509537  -0.00253197  -4.05128e-5   0.00162283   0.00170188  -0.000559519] rtol=1e-4
    @test axe.call[6].fun == :mesh!
    @test axe.call[6].args[1][:,1:10] ≈ [ 0.010907    0.0134416   0.016135    0.0187234  0.0209534  0.0226069  0.0235219  0.0236088    0.0228593   0.0213465;
                                          0.0568619   0.0613883   0.061066    0.0559267  0.0464733  0.0336313  0.0186577  0.00301825  -0.0117562  -0.0242193;
                                         -0.0223258  -0.00756755  0.00783893  0.0223856  0.0346485  0.0434272  0.0478625  0.0475202    0.0424337   0.0331011] rtol = 1e-4
    @test axe.call[6].args[2][1:10,:]' == [  1   1   2   2   3   3   4   4   5   5;  2  22   3  23   4  24   5  25   6  26; 22  21  23  22  24  23  25  24  26  25]
end

end# module

