#module TestFFT
using Test
using Muscade

g(t)  = 1/√(2π)*exp(-1/2*t^2)   
G(ω)  = 1/√(2π)*exp(-1/2*ω^2)  

# G is obtained from the unitary, angular frequency Fourrier transform
# 𝔉, operates only on real ts, which have a Hermitian transform.  So 𝔉 returns only for ω>0, multiplied by 2 for unitariness

p1      = 6
n1      = 2^p1
δt1     = .4
δω1     = Muscade.getδω(n1,δt1)
t1      = δt1*(-n1/2:n1/2-1)
ω1      = δω1*(0:n1/2-1)
X1      = Muscade.𝔉(g.(t1),δt1) 
x1′     = Muscade.𝔉⁻¹(X1,δω1)

p2      = 6
n2      = 2^p2
δt2     = .2
δω2     = Muscade.getδω(n2,δt2)
t2      = δt2*(-n2/2:n2/2-1)
ω2      = δω2*(0:n2/2-1)
X2      = Muscade.𝔉(g.(t2),δt2) 
x2′     = Muscade.𝔉⁻¹(X2,δω2)

p3      = 7
n3      = 2^p3
δt3     = .4
δω3     = Muscade.getδω(n3,δt3)
t3      = δt3*(-n3/2:n3/2-1)
ω3      = δω3*(0:n3/2-1)
X3      = Muscade.𝔉(g.(t3),δt3) 
x3′     = Muscade.𝔉⁻¹(X3,δω3)

# using GLMakie
# fig      = Figure(size = (1000,600))
# display(fig) # open interactive window (gets closed down by "save")

# axe_time      = Axis(fig[1,1],xlabel="t",ylabel="g")
# axe_freq      = Axis(fig[2,1],xlabel="ω",ylabel="|G|")

# with_theme(Theme(fontsize = 30,font="Arial")) do
#    oedge    = lines!(  axe_time,t1,g.(t1)          ,color = :black,linewidth = 1)
#    overtex  = scatter!(axe_time,t1,     x1′        ,markersize=5,color=:green)
#    overtex  = scatter!(axe_time,t2,     x2′        ,markersize=5,color=:red)
#    overtex  = scatter!(axe_time,t3,     x3′        ,markersize=5,color=:orange)
# end
# with_theme(Theme(fontsize = 30,font="Arial")) do
#    oedge    = lines!(  axe_freq,ω1,  G.(ω1)      ,color = :black,linewidth = 1)
#    overtex  = scatter!(axe_freq,ω1,abs.(X1)      ,markersize=10,color=:green)
#    overtex  = scatter!(axe_freq,ω2,abs.(X2)      ,markersize=10,color=:red)
#    overtex  = scatter!(axe_freq,ω3,abs.(X3)      ,markersize=10,color=:orange)
# end


@testset "𝔉" begin
   @test X1[1:4:end] ≈ [0.39894228040143276 - 3.2102472328760894e-14im, 0.24638675308419544 - 1.3287447155944078e-17im, 0.05804157918388552 - 4.9827926834790294e-18im, 0.005215256379624876 - 8.784767309154435e-18im, 0.00017874202879189385 + 0.0im, 2.3366425464292235e-6 - 5.3149788623776313e-17im, 1.1651248216739184e-8 - 1.3287447155944078e-17im, 2.215987424538262e-11 - 7.086638483170176e-17im]
   @test X2[1:4:end] ≈ [0.3989422803309979 + 2.3091315432716196e-11im, 0.05804157911823811 + 1.10728726299534e-18im, 0.00017874197398465492 - 7.569346524382207e-21im, 1.1607540302830413e-8 - 4.262085169816502e-18im, -3.5122106703035916e-11 + 0.0im, -2.93170920089475e-11 - 1.7549532028762994e-17im, -2.5684643488406082e-11 - 1.7301363484302186e-20im, -2.3714292963473583e-11 - 1.771659620792544e-17im]
   @test X3[1:8:end] ≈ [0.39894228040143276 - 3.2102472328760894e-14im, 0.24638675308419553 - 2.6574894311888157e-17im, 0.05804157918388555 + 2.21457452599068e-18im, 0.00521525637962485 - 1.4100611239706282e-17im, 0.00017874202879189385 + 0.0im, 2.3366425464373555e-6 - 5.121203591353447e-18im, 1.16512481979153e-8 + 0.0im, 2.2159812237295893e-11 - 1.771659620792544e-17im]
end
@testset "𝔉⁻¹" begin
   @test x1′[1:8:end]  ≈ [0.0, 5.944950858368643e-18, 5.088140135787997e-10, 0.002384088201464852, 0.39894228040143276, 0.0023840882014648478, 5.08813991837266e-10, 4.416249209073849e-17]
   @test x2′[1:8:end]  ≈ [5.08814046191102e-10, 3.961299091042141e-6, 0.0023840882014648313, 0.1109208346794555, 0.39894228040143276, 0.11092083467945556, 0.0023840882014648755, 3.961299091150848e-6]
   @test x3′[1:16:end] ≈ [0.0, -1.3315999428623215e-17, 0.0, 5.088140083372391e-10, 0.39894228040143276, 5.088140377655672e-10, 0.0, -2.7370741107777605e-17]
end

@testset "𝔉 is unitary" begin
   @test  sum(abs2.(g.(t1)))           .*δt1   ≈ 0.28209479177387814 # energy of signal
   @test (sum(abs2.(X1))-abs2(X1[1])/2).*δω1   ≈ 0.14104739588693924 # energy of half spectre
   @test (sum(abs2.(X2))-abs2(X2[1])/2).*δω2   ≈ 0.14104739588693924 # energy of half spectre
   @test (sum(abs2.(X3))-abs2(X3[1])/2).*δω3   ≈ 0.14104739588693924 # energy of half spectre
end

#end #module
