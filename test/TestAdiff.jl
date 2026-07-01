module TestAdiff
using Muscade
using Test,StaticArrays

const ∂𝕣11= ∂ℝ{1,1,𝕣}
const ∂𝕣12= ∂ℝ{1,2,𝕣}
const ∂𝕣22= ∂ℝ{2,2,𝕣}

## Constructors and promotion
dx1= ∂ℝ{1,2}(3.,SVector(.1,.2))
dx2= ∂ℝ{1,2}(2.,SVector(.3,.4))
dx3= ∂ℝ{1,2}(4.,SVector(.5,.6))
x  = ∂ℝ{2,2}(dx1,SVector(dx2,dx3))
a  = ∂ℝ{2,2}(dx1,SVector(1,2))
b  = ∂ℝ{2,4}(dx1)
c  = ∂ℝ{2,4}(dx1,3)
t1 = promote_rule(typeof(dx1),typeof(3))
t2 = promote_rule(typeof(dx1),typeof(x))
(d,e)=promote(dx1,7)
(g,h)=promote(dx1,x)

## Extraction
dscaled=Muscade.δ_{1,3,𝕣}(SVector(1.,2.,3.))
vscaled=Muscade.variate_{1,3}(SVector(4.,4.,4.),SVector(1.,2.,3.))
Δ   = Muscade.δ_{1,2,𝕣}()
C1  = Muscade.variate_{precedence(Δ)+1,2}(SVector(3.,4.))
C   = Muscade.variate_{precedence(C1)+1,2}(C1)
PC  = precedence(C)
PC1 = precedence(C1)
vC  = value{PC}(C)
∂C  = ∂{PC,2}(  C)
vvC = value{PC1}(vC)
v∂C = value{PC1}(∂C)
∂vC = ∂{PC1,2}(    vC)
#∂∂C = ∂{PC1,2}(    ∂C)
dX1 = Muscade.toggle(false,dx1,3.)

## Operations
oa = Muscade.variate_{1}(1.)
ob = Muscade.variate_{precedence(oa   )+1}(2.)
oc = Muscade.variate_{precedence(oa,ob)+1}(3.)
od = oa+oc
oe = od+ob
oj = od^2
og = 2^od
oh = oa^oc
ok = oc*oa/oc
ox = SVector(1.,2.,3.)
oX = Muscade.variate_{1,3}(ox)


## norm
using LinearAlgebra
nrm = norm(oX)

## atan(s,c)

x3 = Muscade.variate_{1}(1.6455)
s3,c3=sin(x3),cos(x3)
y3 = atan(s3,c3)

##
@testset "Adiff" begin
    @testset "Adiff construct and promote" begin
        @test dscaled[2] === ∂ℝ{1, 3, Float64}(0.0, SVector(0.0, 2.0, 0.0))
        @test vscaled[2] === ∂ℝ{1, 3, Float64}(4.0, SVector(0.0, 2.0, 0.0))
        @test dx1        === ∂𝕣12(3.0, SVector(0.1, 0.2))
        @test x          === ∂ℝ{2,2,∂𝕣12}(∂𝕣12(3.0, SVector(0.1, 0.2)), SVector(∂𝕣12(2.0, SVector(0.3, 0.4)), ∂𝕣12(4.0, SVector(0.5, 0.6))))
        @test a          === ∂ℝ{2,2,∂𝕣12}(∂𝕣12(3.0, SVector(0.1, 0.2)), SVector(∂𝕣12(1.0, SVector(0.0, 0.0)), ∂𝕣12(2.0, SVector(0.0, 0.0))))
        @test b          === ∂ℝ{2,4,∂𝕣12}(∂𝕣12(3.0, SVector(0.1, 0.2)), SVector(∂𝕣12(0.0, SVector(0.0, 0.0)), ∂𝕣12(0.0, SVector(0.0, 0.0)), ∂𝕣12(0.0, SVector(0.0, 0.0)), ∂𝕣12(0.0, SVector(0.0, 0.0))))
        @test c          === ∂ℝ{2,4,∂𝕣12}(∂𝕣12(3.0, SVector(0.1, 0.2)), SVector(∂𝕣12(0.0, SVector(0.0, 0.0)), ∂𝕣12(0.0, SVector(0.0, 0.0)), ∂𝕣12(1.0, SVector(0.0, 0.0)), ∂𝕣12(0.0, SVector(0.0, 0.0))))
        @test t1         == ∂𝕣12
        @test t2         ==∂ℝ{2,2,∂𝕣12}
        @test d          === dx1
        @test e          === ∂𝕣12(7.0, SVector(0.0, 0.0))
        @test typeof(g) == ∂ℝ{2,2,∂𝕣12}
        @test g         === ∂ℝ{2,2,∂𝕣12}(∂𝕣12(3.0, SVector(0.0, 0.0)), SVector(∂𝕣12(0.1, SVector(0.0, 0.0)), ∂𝕣12(0.2, SVector(0.0, 0.0))))
        @test h          === x
    end

    @testset "Adiff extraction" begin
        @test Δ            === SVector(∂𝕣12(0.0, SVector(1.0, 0.0)),∂𝕣12(0.0, SVector(0.0, 1.0)))
        @test precedence(Δ) == 1
        @test typeof(C)    == SVector{2, ∂ℝ{3, 2, ∂ℝ{2, 2, Float64}}}#Array{∂ℝ{3,2,∂𝕣22},1}
        @test C[1]         === ∂ℝ{3,2,∂𝕣22}(∂𝕣22(3.0, SVector(1.0, 0.0)), SVector(∂𝕣22(1.0, SVector(0.0, 0.0)), ∂𝕣22(0.0, SVector(0.0, 0.0))))
        @test vC[1]        === ∂𝕣22(3.0, SVector(1.0, 0.0))
        @test ∂C[1]        === ∂𝕣22(1.0, SVector(0.0, 0.0))
        @test vvC          === SVector(3.0,4.0)
        @test v∂C          === SMatrix{2,2}(1.0,0.0,0.0,1.0)
        @test ∂vC          === ∂vC
#        @test ∂∂C          === zeros(2,2,2)   broken=true # extract partial derivatives of an SMatrix or higher order
        @test typeof(dX1)  == ∂𝕣12
    end

    @testset "Adiff operations" begin
        @test oa === ∂𝕣11(1.0, SVector(1.0))
        @test ob === ∂ℝ{2,1,𝕣}(2, SVector(1))
        @test oc === ∂ℝ{3,1,𝕣}(3.0, SVector(1.0))
        @test od === ∂ℝ{3,1,∂𝕣11}(∂𝕣11(4.0, SVector(1.0)), SVector(∂𝕣11(1.0, SVector(0.0))))
        @test oe === ∂ℝ{3,1,∂ℝ{2,1,∂𝕣11}}(∂ℝ{2,1,∂𝕣11}(∂𝕣11(6.0, SVector(1.0)), SVector(∂𝕣11(1.0, SVector(0.0)))), SVector(     ∂ℝ{2,1,∂𝕣11}(∂𝕣11(1.0, SVector(0.0)), SVector(∂𝕣11(0.0, SVector(0.0))))))
        @test od === ∂ℝ{3,1,∂𝕣11}(∂𝕣11(4.0, SVector(1.0)), SVector(∂𝕣11(1.0, SVector(0.0))))
        @test og === ∂ℝ{3,1,∂𝕣11}(∂𝕣11(15.999999999999998, SVector(11.090354888959123)), SVector(∂𝕣11(11.090354888959123, SVector(7.687248222691221)))) 
        @test oj === od*od
        @test ok === ∂ℝ{3,1,∂𝕣11}(∂𝕣11(1.0, SVector(1.0)), SVector(∂𝕣11(0.0, SVector(0.0))))
        @test value{1}(2*oX)==2*ox
        @test Muscade.variate_{1}(1.)^0 === ∂ℝ{1,1,𝕣}(1. , SVector(0.)) 
        @test Muscade.variate_{1}(0.)^0 === ∂ℝ{1,1,𝕣}(1. , SVector(0.)) # this can of course be discussed...
    end


    @testset " norm" begin
        @test nrm === sqrt(sum(oX.^2))
    end

    @testset "atan" begin
        y3 === ∂ℝ{1, 1, 𝕣}(1.6455, [1.0])
    end

    θ = [0,1e-8,1e-6,1e-2,.1,1,π]
    @testset "sinc1" begin
        @test Muscade.sinc1.( θ,sin.(θ),cos.(θ),inv.(θ),θ.*θ)≈[1.0,1.0,0.9999999999998334,0.9999833334166666,0.9983341664682817,0.8414709848078965,0.0]
        @test Muscade.sinc1′.(θ,sin.(θ),cos.(θ),inv.(θ),θ.*θ)≈[ 0.0,-3.3333333333333334e-9,-3.3333333333329995e-7,-0.003333300000107897,-0.03330001190255594,-0.30116867893975674,-0.3183098861837907]
        @test Muscade.sinc1″.(θ,sin.(θ),cos.(θ),inv.(θ),θ.*θ)≈[-0.3333333333333333, -0.3333333333333333, -0.33333333333323334, -0.33332333339285697, -0.33233392841710757, -0.23913362692838303, 0.20264236728467552]
        @test Muscade.sinc1‴.(θ,sin.(θ),cos.(θ),inv.(θ),θ.*θ)≈[0.0, 2.0e-9, 1.999999999999762e-7, 0.001999976190568783, 0.019976199733646196, 0.1770985749170091, 0.12480067958459377]
        @test Muscade.sinc1⁗.(θ,sin.(θ),cos.(θ),inv.(θ),θ.*θ)≈[0.2,0.2,0.1999999999999286,0.19999285718915333,0.19928617712243368,0.13307670326700266,-0.14309499132952147]
    end



end # testset Adiff

end
