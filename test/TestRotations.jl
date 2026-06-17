module TestRotations

using Test, StaticArrays, LinearAlgebra
using Muscade, Muscade.Toolbox

# using GLMakie
# fig      = Figure(size = (1200,1000))
# display(fig) # open interactive window (gets closed down by "save")
# axe1      = Axis(fig[1,1],title="sinc1(x)",xlabel="x",ylabel="Y")
# axe2      = Axis(fig[2,1],title="sinc1′(x)",xlabel="x",ylabel="Y")
# axe3      = Axis(fig[3,1],title="sinc1″(x)",xlabel="x",ylabel="Y")
# axe4      = Axis(fig[4,1],title="sinc1‴(x)",xlabel="x",ylabel="Y")
# axe5      = Axis(fig[5,1],title="sinc1⁗(x)",xlabel="x",ylabel="Y")
# x = range(0,π,1001)
# with_theme(Theme(fontsize = 30,font="Arial")) do
#     lines!(  axe1,x,sinc1.(x) ,color = :black)
#     lines!(  axe2,x,sinc1′.(x) ,color = :black)
#     lines!(  axe3,x,sinc1″.(x) ,color = :black)
#     lines!(  axe4,x,sinc1‴.(x) ,color = :black)
#     lines!(  axe5,x,sinc1⁗.(x) ,color = :black)
# end

function scac1(x)
    X = variate{1}(x)
    Y = Toolbox.scac(X)
    return ∂{1}(Y)
end
function scac2(x)
    X = variate{2}(variate{1}(x))
    Y = Toolbox.scac(X)
    return ∂{1}(∂{2}(Y))
end
function scac3(x)
    X = variate{3}(variate{2}(variate{1}(x)))
    Y = Toolbox.scac(X)
    return ∂{1}(∂{2}(∂{3}(Y)))
end
function scac4(x)
    X = variate{4}(variate{3}(variate{2}(variate{1}(x))))
    Y = Toolbox.scac(X)
    return ∂{1}(∂{2}(∂{3}(∂{4}(Y))))
end
function scac5(x)
    X = variate{5}(variate{4}(variate{3}(variate{2}(variate{1}(x)))))
    Y = Toolbox.scac(X)
    return ∂{1}(∂{2}(∂{3}(∂{4}(∂{5}(Y)))))
end

x = [-1+1e-11,-1+1e-9,0,1-1e-9,1-1e-11,1]
@testset "Toolbox.scac" begin
    @test Toolbox.scac.(x)≈[1.4235271721825872e-6,1.4235453308738641e-5,0.6366197723675814,0.9999999996666666,0.9999999999966667,1.0]
    @test scac1.(x)        ≈[71176.45403923128, 7117.8281743984935, 0.40528473456935105, 0.33333333337777776, 0.33333333333377774, 0.3333333333333333]
    @test scac2.(x)        ≈[  -3.55881227537807e15,    -3.5588128658975117e12,    -0.12059522143638957,    -0.04444444447619157,    -0.04444444444476192,    -0.044444444444444446]
    @test scac3.(x)        ≈[ 5.338217971349224e26,    5.338219446646726e21,    0.17496482731099405,    0.03174712798947714,    0.031747127950916276,    0.031747127950526775] 
    @test scac4.(x)        ≈[-1.334554382408856e38, -1.33455489871291e31, 0., -0.038950361148861766, -0.03895036108203678, -0.038950361081361774]
end
# fig      = Figure(size = (1200,1000))
# display(fig) # open interactive window (gets closed down by "save")
# axe1      = Axis(fig[1,1],title="Toolbox.scac(x)",xlabel="x",ylabel="Y")
# axe2      = Axis(fig[2,1],title="scac1(x)",xlabel="x",ylabel="Y")
# axe3      = Axis(fig[3,1],title="scac2(x)",xlabel="x",ylabel="Y")
# axe4      = Axis(fig[4,1],title="scac3(x)",xlabel="x",ylabel="Y")
# axe5      = Axis(fig[5,1],title="scac4(x)",xlabel="x",ylabel="Y")
# axe6      = Axis(fig[6,1],title="scac5(x)",xlabel="x",ylabel="Y")
# x = range(.98,1,1001)
# with_theme(Theme(fontsize = 30,font="Arial")) do
#     lines!(  axe1,x,Toolbox.scac.(x) ,color = :black)
#     lines!(  axe2,x,scac1.(x) ,color = :black)
#     lines!(  axe3,x,scac2.(x) ,color = :black)
#     lines!(  axe4,x,scac3.(x) ,color = :black)
#     lines!(  axe5,x,scac4.(x) ,color = :black)
#     lines!(  axe6,x,scac5.(x) ,color = :black)
# end

a = SA[1,0,0]
b = SA[0,1,1]
r = Toolbox.adjust(a,b)
R = Toolbox.Rodrigues(r)
u = R*a

v1      = variate{1,3}(SA[.1,.2,.3])
M1      = Toolbox.Rodrigues(v1)
w1,w∂v1 = value_∂{1,3}(Toolbox.Rodrigues⁻¹(M1))

v2      = variate{1,3}(SA[1e-7,2e-7,1e-8])
M2      = Toolbox.Rodrigues(v2)
w2,w∂v2 = value_∂{1,3}(Toolbox.Rodrigues⁻¹(M2))

v3      = variate{1,3}(SA[0.,0.,0.])
M3      = Toolbox.Rodrigues(v3)
w3,w∂v3 = value_∂{1,3}(Toolbox.Rodrigues⁻¹(M3))


@testset "rotations" begin
    @test r ≈ [0.0, -1.1107207345395913, 1.1107207345395913]
    @test u ≈ [2.220446049250313e-16, 0.7071067811865476, 0.7071067811865476]
    @test v1 ≈ w1
    @test w∂v1 ≈ I
    @test v2 ≈ w2
    @test w∂v2 ≈ I
    @test v3 ≈ w3
    @test w∂v3 ≈ I
end

end


