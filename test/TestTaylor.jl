module TestTaylor
using Muscade
using Test,StaticArrays

const X1 = (SVector{3,𝕣}(1,2,3),)
const X2 = (SVector{3,𝕣}(1,2,3),SVector{3,𝕣}(4,5,6))
const X3 = (SVector{3,𝕣}(1,2,3),SVector{3,𝕣}(4,5,6),SVector{3,𝕣}(7,8,9))
const X4 = (Muscade.variate_{1,3}(SVector{3,𝕣}(1,2,3)),)
const X5 = (Muscade.variate_{1,3}(SVector{3,𝕣}(1,2,3)),Muscade.variate_{1,3}(SVector{3,𝕣}(4,5,6)))
const X6 = (Muscade.variate_{1,3}(SVector{3,𝕣}(1,2,3)),Muscade.variate_{1,3}(SVector{3,𝕣}(4,5,6)),Muscade.variate_{1,3}(SVector{3,𝕣}(7,8,9)))

Y1=motion{2}(X1)
Y2=motion{2}(X2)
Y3=motion{2}(X3)
Y4=motion{3}(X4)
Y5=motion{3}(X5)
Y6=motion{3}(X6)

@testset "motion" begin
    @test typeof(X1) == Tuple{SVector{3, 𝕣}}
    @test typeof(Y1) == SVector{3, 𝕣}

    @test typeof(X2) == Tuple{SVector{3, 𝕣}, SVector{3, 𝕣}}
    @test typeof(Y2) == SVector{3, ∂ℝ{2, 1, 𝕣}} 

    @test typeof(X3) == Tuple{SVector{3, 𝕣}, SVector{3, 𝕣}, SVector{3, 𝕣}}
    @test typeof(Y3) == SVector{3, ∂ℝ{3, 1, ∂ℝ{2, 1, 𝕣}}}

    @test typeof(X4) == Tuple{SVector{3, ∂ℝ{1, 3, 𝕣}}}
    @test typeof(Y4) == SVector{3, ∂ℝ{1, 3, 𝕣}}

    @test typeof(X5) == Tuple{SVector{3, ∂ℝ{1, 3, 𝕣}}, SVector{3, ∂ℝ{1, 3, 𝕣}}}
    @test typeof(Y5) == SVector{3, ∂ℝ{3, 1, ∂ℝ{1, 3, 𝕣}}}

    @test typeof(X6) == Tuple{SVector{3, ∂ℝ{1, 3, 𝕣}}, SVector{3, ∂ℝ{1, 3, 𝕣}}, SVector{3, ∂ℝ{1, 3, 𝕣}}}
    @test typeof(Y6) == SVector{3, ∂ℝ{4, 1, ∂ℝ{3, 1, ∂ℝ{1, 3, 𝕣}}}}
end

@testset "motion⁻¹" begin
    @test motion⁻¹{2,1,0}(Y1) === SVector{3,𝕣}(1,2,3)
    @test motion⁻¹{2,1,1}(Y1) === SVector{3,𝕣}(0,0,0)
    @test motion⁻¹{2,1,2}(Y1) === SVector{3,𝕣}(0,0,0)

    @test motion⁻¹{2,2,0}(Y2) === SVector{3,𝕣}(1,2,3)
    @test motion⁻¹{2,2,1}(Y2) === SVector{3,𝕣}(4,5,6)
    @test motion⁻¹{2,2,2}(Y2) === SVector{3,𝕣}(0,0,0)

    @test motion⁻¹{2,3,0}(Y3) === SVector{3,𝕣}(1,2,3)
    @test motion⁻¹{2,3,1}(Y3) === SVector{3,𝕣}(4,5,6)
    @test motion⁻¹{2,3,2}(Y3) === SVector{3,𝕣}(7,8,9)

    @test motion⁻¹{3,1,0}(Y4) === Muscade.variate_{1,3}(SVector{3,𝕣}(1,2,3))
    @test motion⁻¹{3,1,1}(Y4) === SVector{3,𝕣}(0,0,0)
    @test motion⁻¹{3,1,2}(Y4) === SVector{3,𝕣}(0,0,0)

    @test motion⁻¹{3,2,0}(Y5) === Muscade.variate_{1,3}(SVector{3,𝕣}(1,2,3))
    @test motion⁻¹{3,2,1}(Y5) === Muscade.variate_{1,3}(SVector{3,𝕣}(4,5,6))
    @test motion⁻¹{3,2,2}(Y5) === SVector{3,𝕣}(0,0,0)

    @test motion⁻¹{3,3,0}(Y6) === Muscade.variate_{1,3}(SVector{3,𝕣}(1,2,3))
    @test motion⁻¹{3,3,1}(Y6) === Muscade.variate_{1,3}(SVector{3,𝕣}(4,5,6))
    @test motion⁻¹{3,3,2}(Y6) === Muscade.variate_{1,3}(SVector{3,𝕣}(7,8,9))

    @test motion⁻¹{2,1}(Y1) === (SVector{3,𝕣}(1,2,3),)
    @test motion⁻¹{2,2}(Y2) === (SVector{3,𝕣}(1,2,3),SVector{3,𝕣}(4,5,6))
    @test motion⁻¹{2,3}(Y3) === (SVector{3,𝕣}(1,2,3),SVector{3,𝕣}(4,5,6),SVector{3,𝕣}(7,8,9))
end

a = SVector(3.,4.)
@testset "revariate" begin
    @test Muscade.revariate{0}(a)    === a
    @test Muscade.revariate{1}(a)[1] ===                                                                   ∂ℝ{1, 2, 𝕣}(3.0, [1.0, 0.0])
    @test Muscade.revariate{2}(a)[1] ===                                       ∂ℝ{2, 2, ∂ℝ{1, 2, 𝕣}}(∂ℝ{1, 2, 𝕣}(3.0, [1.0, 0.0]), ∂ℝ{1, 2, 𝕣}[∂ℝ{1, 2, 𝕣}(1.0, [0.0, 0.0]), ∂ℝ{1, 2, 𝕣}(0.0, [0.0, 0.0])])
    @test Muscade.revariate{3}(a)[1] === ∂ℝ{3, 2, ∂ℝ{2, 2, ∂ℝ{1, 2, 𝕣}}}(∂ℝ{2, 2, ∂ℝ{1, 2, 𝕣}}(∂ℝ{1, 2, 𝕣}(3.0, [1.0, 0.0]), ∂ℝ{1, 2, 𝕣}[∂ℝ{1, 2, 𝕣}(1.0, [0.0, 0.0]), ∂ℝ{1, 2, 𝕣}(0.0, [0.0, 0.0])]), ∂ℝ{2, 2, ∂ℝ{1, 2, 𝕣}}[∂ℝ{2, 2, ∂ℝ{1, 2, 𝕣}}(∂ℝ{1, 2, 𝕣}(1.0, [0.0, 0.0]), ∂ℝ{1, 2, 𝕣}[∂ℝ{1, 2, 𝕣}(0.0, [0.0, 0.0]), ∂ℝ{1, 2, 𝕣}(0.0, [0.0, 0.0])]), ∂ℝ{2, 2, ∂ℝ{1, 2, 𝕣}}(∂ℝ{1, 2, 𝕣}(0.0, [0.0, 0.0]), ∂ℝ{1, 2, 𝕣}[∂ℝ{1, 2, 𝕣}(0.0, [0.0, 0.0]), ∂ℝ{1, 2, 𝕣}(0.0, [0.0, 0.0])])])
end

f(X) = SVector(cos(X[1])*X[3],sin(X[2])*X[3])
g(X) = 4.5*X[1]^2/2 + X[1] + 2.
h(X) = SVector(3.,4.)
k(X) = 5.
w(X) = (f(X),g(X),h(X),k(X))
X₀   = SVector(0.,0.,1.)
vX₀  = Muscade.variate_{1,3}(X₀)
wX₀  = Muscade.variate_{2,3}(vX₀)
@testset "Muscade.fast" begin
    @test Muscade.apply{:chainrule}(w, X₀) === w( X₀)
    @test Muscade.apply{:chainrule}(w,vX₀) === w(vX₀)
    @test Muscade.apply{:chainrule}(w,wX₀) === w(wX₀)
end

yy    = Muscade.to_order{1,3}((revariate{1}(SVector(4.,5.,6.)),revariate{2}(SVector(7.,8.,9.)) ))  
@testset "firstorderonly" begin
    @test yy[1] === revariate{1}(SVector(4.,5.,6.)) 
    @test yy[2] === revariate{1}(SVector(7.,8.,9.)) 
    @test Muscade.to_order{0,0}(3.) == 3.0
    @test Muscade.to_order{1,0}(3.) == ∂ℝ{1, 0, Float64}(3.0, Float64[])
    @test Muscade.to_order{2,0}(3.) == ∂ℝ{2, 0, ∂ℝ{1, 0, Float64}}(∂ℝ{1, 0, Float64}(3.0, Float64[]), ∂ℝ{1, 0, Float64}[])
end



#### chainrule with NamedTuple

x      = SVector(1.,2.,2.5,3.)
X      = Muscade.variate_{1,4}(x)
ε      = SMatrix{2,2}((X.^2)...)
eleres = (part=(ε = ε, x = X), y = 2x[2],z = 3.)

cost(eleres) = sum(eleres.part.ε)+eleres.y

Neleres = Muscade.flat_length(eleres)
Teleres = Muscade.flat_eltype(eleres)
Peleres = Muscade.precedence(eleres)
Feleres = Muscade.flatten(eleres)

@testset "flatten" begin
    @test Neleres               == 10
    @test Teleres               == ∂ℝ{1,4,𝕣}
    @test Peleres               == 1
    @test length(Feleres)       == 10
    @test eltype(Feleres)       == ∂ℝ{1,4,𝕣}
    @test Muscade.precedence(Feleres)   == 1
end

# eleres, P=1 comes from the element
# 4: Ndof
# 10: Neleres
Releres  = Muscade.revariate{2}(eleres)
Rq       = cost(Releres)
q        = Muscade.chainrule(Rq,Muscade.to_order{2,4}(eleres))   
q2       = cost(Muscade.to_order{2,4}(eleres))

@testset "chainrule NamedTuple" begin
    @test Muscade.flat_eltype(Muscade.revariate{2}(eleres))             == ∂ℝ{2, 10, ∂ℝ{1, 10, 𝕣}}
    @test Muscade.flat_eltype(Rq)                                       == ∂ℝ{2, 10, ∂ℝ{1, 10, 𝕣}}
    @test Muscade.flat_eltype(q)                                        == ∂ℝ{2, 4 , ∂ℝ{1, 4 , 𝕣}} 
    @test Muscade.flat_eltype(Muscade.to_order{2,4}(eleres))            == ∂ℝ{2, 4 , ∂ℝ{1, 4 , 𝕣}} 
    @test Muscade.flat_eltype(q2)                                       == ∂ℝ{2, 4 , ∂ℝ{1, 4 , 𝕣}} 
    @test q == q2
end

@testset "inferred" begin
    @inferred revariate{2}(eleres)
    @inferred Muscade.to_order{2,4}(Muscade.flatten(eleres))
    @inferred Muscade.chainrule(Rq,Muscade.to_order{2,4}(Muscade.flatten(eleres)))
end



@testset "multiple derivative" begin
    a       = SVector(1.,2.,3.)
    b       = SVector(3.,4.)
    Pa,Na,A = variate{2}(a)
    Pb,Nb,B = variate(b) # A and B not compatible
    _ ,_ ,C = variate{2}(a,scale=AllElements(2.))
    _ ,_ ,D = variate{2}(a,scale=SVector(2.,2.,2.))
    A0      = Muscade.zerovalue(A[1])

    @test Muscade.npartials((a=A,b=b)) == (3,3) 
    @test_throws AssertionError Muscade.npartials((a=A,b=B))   

    @test Pa==2
    @test Na==3
    @test A0  ===∂ℝ{2,3,∂ℝ{1,3,𝕣}}(∂ℝ{1,3,𝕣}(0,[1,0,0]),∂ℝ{1,3,𝕣}[∂ℝ{1,3,𝕣}(1,[0,0,0]),∂ℝ{1,3,𝕣}(0,[0,0,0]),∂ℝ{1,3,𝕣}(0,[0,0,0])])
    @test A[1]===∂ℝ{2,3,∂ℝ{1,3,𝕣}}(∂ℝ{1,3,𝕣}(1,[1,0,0]),∂ℝ{1,3,𝕣}[∂ℝ{1,3,𝕣}(1,[0,0,0]),∂ℝ{1,3,𝕣}(0,[0,0,0]),∂ℝ{1,3,𝕣}(0,[0,0,0])])
    @test A[2]===∂ℝ{2,3,∂ℝ{1,3,𝕣}}(∂ℝ{1,3,𝕣}(2,[0,1,0]),∂ℝ{1,3,𝕣}[∂ℝ{1,3,𝕣}(0,[0,0,0]),∂ℝ{1,3,𝕣}(1,[0,0,0]),∂ℝ{1,3,𝕣}(0,[0,0,0])])
    @test A[3]===∂ℝ{2,3,∂ℝ{1,3,𝕣}}(∂ℝ{1,3,𝕣}(3,[0,0,1]),∂ℝ{1,3,𝕣}[∂ℝ{1,3,𝕣}(0,[0,0,0]),∂ℝ{1,3,𝕣}(0,[0,0,0]),∂ℝ{1,3,𝕣}(1,[0,0,0])])
    @test Pb==1
    @test Nb==2
    @test B[1]===∂ℝ{1,2,𝕣}(3,[1,0])
    @test B[2]===∂ℝ{1,2,𝕣}(4,[0,1])
    @test C[1]===∂ℝ{2,3,∂ℝ{1,3,𝕣}}(∂ℝ{1,3,𝕣}(1,[2,0,0]),∂ℝ{1,3,𝕣}[∂ℝ{1,3,𝕣}(2,[0,0,0]),∂ℝ{1,3,𝕣}(0,[0,0,0]),∂ℝ{1,3,𝕣}(0,[0,0,0])])
    @test C[2]===∂ℝ{2,3,∂ℝ{1,3,𝕣}}(∂ℝ{1,3,𝕣}(2,[0,2,0]),∂ℝ{1,3,𝕣}[∂ℝ{1,3,𝕣}(0,[0,0,0]),∂ℝ{1,3,𝕣}(2,[0,0,0]),∂ℝ{1,3,𝕣}(0,[0,0,0])])
    @test C[3]===∂ℝ{2,3,∂ℝ{1,3,𝕣}}(∂ℝ{1,3,𝕣}(3,[0,0,2]),∂ℝ{1,3,𝕣}[∂ℝ{1,3,𝕣}(0,[0,0,0]),∂ℝ{1,3,𝕣}(0,[0,0,0]),∂ℝ{1,3,𝕣}(2,[0,0,0])])
    @test D[1]===∂ℝ{2,3,∂ℝ{1,3,𝕣}}(∂ℝ{1,3,𝕣}(1,[2,0,0]),∂ℝ{1,3,𝕣}[∂ℝ{1,3,𝕣}(2,[0,0,0]),∂ℝ{1,3,𝕣}(0,[0,0,0]),∂ℝ{1,3,𝕣}(0,[0,0,0])])
    @test D[2]===∂ℝ{2,3,∂ℝ{1,3,𝕣}}(∂ℝ{1,3,𝕣}(2,[0,2,0]),∂ℝ{1,3,𝕣}[∂ℝ{1,3,𝕣}(0,[0,0,0]),∂ℝ{1,3,𝕣}(2,[0,0,0]),∂ℝ{1,3,𝕣}(0,[0,0,0])])
    @test D[3]===∂ℝ{2,3,∂ℝ{1,3,𝕣}}(∂ℝ{1,3,𝕣}(3,[0,0,2]),∂ℝ{1,3,𝕣}[∂ℝ{1,3,𝕣}(0,[0,0,0]),∂ℝ{1,3,𝕣}(0,[0,0,0]),∂ℝ{1,3,𝕣}(2,[0,0,0])])
end


@testset "joint derivative" begin
    a       = SVector(1.,2.,3.)
    b       = SVector(3.,4.)
    P,N,(A,B) = variate((a,b)) 

    @test P==1
    @test N==5
    @test A[1]===∂ℝ{1,5,𝕣}(1,[1,0,0,0,0])
    @test A[2]===∂ℝ{1,5,𝕣}(2,[0,1,0,0,0])
    @test A[3]===∂ℝ{1,5,𝕣}(3,[0,0,1,0,0])
    @test B[1]===∂ℝ{1,5,𝕣}(3,[0,0,0,1,0])
    @test B[2]===∂ℝ{1,5,𝕣}(4,[0,0,0,0,1])
end

@testset "nested derivative" begin
    a       = SVector(1.,2.,3.)
    b       = SVector(3.,4.)
    Pa,Na,A  = variate{1}(a)
    Pb,Nb,B  = variate{1}(b,constants=(A,0.)) 
    _ ,_ ,Bs = variate{1}(b,constants=A,scale=AllElements(2.)) 

    Pb0,Nb0,B0 = variate0{1}(b,constants=A) 
    @test Pb==2
    @test Nb==2
    @test B[1]===∂ℝ{2,2,∂ℝ{1,3,𝕣}}(∂ℝ{1,3,𝕣}(3,[0,0,0]),∂ℝ{1,3,𝕣}[∂ℝ{1,3,𝕣}(1,[0,0,0]),∂ℝ{1,3,𝕣}(0,[0,0,0])])
    @test B[2]===∂ℝ{2,2,∂ℝ{1,3,𝕣}}(∂ℝ{1,3,𝕣}(4,[0,0,0]),∂ℝ{1,3,𝕣}[∂ℝ{1,3,𝕣}(0,[0,0,0]),∂ℝ{1,3,𝕣}(1,[0,0,0])])
    @test Muscade.npartials(B)==(2,3)

    @test Pb0==2
    @test Nb0==2
    @test B0[1]===∂ℝ{2,2,∂ℝ{1,3,𝕣}}(∂ℝ{1,3,𝕣}(0,[0,0,0]),∂ℝ{1,3,𝕣}[∂ℝ{1,3,𝕣}(1,[0,0,0]),∂ℝ{1,3,𝕣}(0,[0,0,0])])
    @test B0[2]===∂ℝ{2,2,∂ℝ{1,3,𝕣}}(∂ℝ{1,3,𝕣}(0,[0,0,0]),∂ℝ{1,3,𝕣}[∂ℝ{1,3,𝕣}(0,[0,0,0]),∂ℝ{1,3,𝕣}(1,[0,0,0])])

    @test Bs[1]===∂ℝ{2,2,∂ℝ{1,3,𝕣}}(∂ℝ{1,3,𝕣}(3,[0,0,0]),∂ℝ{1,3,𝕣}[∂ℝ{1,3,𝕣}(2,[0,0,0]),∂ℝ{1,3,𝕣}(0,[0,0,0])])
    @test Bs[2]===∂ℝ{2,2,∂ℝ{1,3,𝕣}}(∂ℝ{1,3,𝕣}(4,[0,0,0]),∂ℝ{1,3,𝕣}[∂ℝ{1,3,𝕣}(0,[0,0,0]),∂ℝ{1,3,𝕣}(2,[0,0,0])])
end

@testset "revariate" begin
    a         = SVector(1.,2.,3.)
    b         = SVector(3.,4.)
    P,N,(A,B) = variate{2}((a,b),scale=AllElements(AllElements(2.))) 
    B_        = revariate(B)
    Bm1       = revariate{1}(B)
    A2,B2     = revariate((A,B))
    @test B_[1]===∂ℝ{2,2,∂ℝ{1,2,𝕣}}(∂ℝ{1,2,𝕣}(3,[1,0]),∂ℝ{1,2,𝕣}[∂ℝ{1,2,𝕣}(1,[0,0]),∂ℝ{1,2,𝕣}(0,[0,0])])
    @test B_[2]===∂ℝ{2,2,∂ℝ{1,2,𝕣}}(∂ℝ{1,2,𝕣}(4,[0,1]),∂ℝ{1,2,𝕣}[∂ℝ{1,2,𝕣}(0,[0,0]),∂ℝ{1,2,𝕣}(1,[0,0])])
    @test B2[1]==∂ℝ{2,5,∂ℝ{1,5,𝕣}}(∂ℝ{1,5,𝕣}(3,[0,0,0,1,0]),∂ℝ{1,5,𝕣}[∂ℝ{1,5,𝕣}(0,[0,0,0,0,0]),∂ℝ{1,5,𝕣}(0,[0,0,0,0,0]),∂ℝ{1,5,𝕣}(0,[0,0,0,0,0]),∂ℝ{1,5,𝕣}(1,[0,0,0,0,0]),∂ℝ{1,5,𝕣}(0,[0,0,0,0,0])])
    @test Bm1[1]===∂ℝ{1,2,𝕣}(3,[1,0])
end


end # module
