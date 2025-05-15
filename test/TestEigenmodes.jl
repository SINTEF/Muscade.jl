#module TestEigenmodes
using Muscade,Test,SparseArrays,Random

"""
    m = MAC(u,v)

Compute the Modal Assurance Criterion for quantitative comparison between two eigenvectors/modeshapes.
"""   
function MAC(v1::AbstractVector, v2::AbstractVector)
    return  dot(v1,v2)^2 ./ (dot(v1,v1)*dot(v2,v2))
end

N      = 1000
K      = spdiagm(N,N,-1=>range(-N,-N,N-1),0=>range(2N,2N,N), 1=>range(-N,-N,N-1))
M      = spdiagm(N,N,0=>range(1/N,1/N,N))

Random.seed!(1234)
seed = randn(N)




neig   = 5
ω²,v,ncv = Muscade.geneig{:SDP}(K,M,neig;seed)

@testset "eigenmodes SDP" begin
    @test ncv ≥ neig
    @test sqrt.(ω²[1:neig]) ≈ [3.1384529113316875, 6.276898094305056, 9.41532782058609, 12.55373436187465, 15.692109989929111]
    @test MAC(v[1][1:200:end],      [  0.00014028558300245532,    0.026364134154732224,    0.042537225116976426,    0.04249391630968802,    0.02625071828168894]) ≈ 1.000 rtol=1e-6
    @test MAC(v[2][1:200:end],      [  0.00028056978420787495,    0.042580114937435326,    0.02613704384203099,   -0.026364134154727957,   -0.04249391630968351]) ≈ 1.000 rtol=1e-6
    @test maximum(abs.((K-ω²[1]*M)*v[1])./abs.(K*v[1])) < 1e-10
    @test maximum(abs.((K-ω²[2]*M)*v[2])./abs.(K*v[2])) < 1e-10
    @test maximum(abs.((K-ω²[3]*M)*v[3])./abs.(K*v[3])) < 1e-10
    @test maximum(abs.((K-ω²[4]*M)*v[4])./abs.(K*v[4])) < 1e-10
    @test maximum(abs.((K-ω²[5]*M)*v[5])./abs.(K*v[5])) < 2e-10
end

ω²,v,ncv = Muscade.geneig{:Hermitian}(K,M,neig;seed)
@testset "eigenmodes Hermitian" begin
    @test ncv ≥ neig
    @test sqrt.(ω²[1:neig]) ≈ [3.138452911330836, 6.276898094304366, 9.415327820585714, 12.55373436187439, 15.692109989929069]
    @test MAC(v[1][1:200:end],      [  0.0002811008347510544,    0.026372151256201664,    0.04254237814552391,    0.042489472779309846,    0.02624087646034281])        ≈ 1.000 rtol=1e-5
    @test MAC(v[2][1:200:end],      [ -0.0005640094745986227,    -0.04274190474190939,    -0.02628003327036983,     0.026212397939501447,     0.042340762998788375 ])   ≈ 1.000 rtol=1e-4
    @test maximum(abs.((K-ω²[1]*M)*v[1])./abs.(K*v[1])) < 1e-10
    @test maximum(abs.((K-ω²[2]*M)*v[2])./abs.(K*v[2])) < 1e-10
    @test maximum(abs.((K-ω²[3]*M)*v[3])./abs.(K*v[3])) < 1e-10
    @test maximum(abs.((K-ω²[4]*M)*v[4])./abs.(K*v[4])) < 1e-10
    @test maximum(abs.((K-ω²[5]*M)*v[5])./abs.(K*v[5])) < 2e-10
end

ω²,v,ncv = Muscade.geneig{:Complex}(K,M,neig;seed)
@testset "eigenmodes Unsymmetric" begin
    @test ncv ≥ neig
    @test sqrt.(real(ω²[1:neig])) ≈ [3.138452911330836, 6.276898094304366, 9.415327820585714, 12.55373436187439, 15.692109989929069]
    @test MAC(real.(v[1][1:200:end]),[ -0.0002805684107607904,    -0.026364133746437556,    -0.042537225451125094,    -0.04249391664172418,    -0.026250717868500813 ]) ≈ 1.000 rtol=1e-4
    @test MAC(real.(v[2][1:200:end]),[ -0.0005611175282438147,    -0.042580114600632754,    -0.026137046653016693,     0.026364136984402934,     0.042493915923130676 ])≈ 1.000 rtol=1e-4
    @test maximum(abs.((K-ω²[1]*M)*v[1])./abs.(K*v[1])) < 1e-10
    @test maximum(abs.((K-ω²[2]*M)*v[2])./abs.(K*v[2])) < 1e-10
    @test maximum(abs.((K-ω²[3]*M)*v[3])./abs.(K*v[3])) < 1e-10
    @test maximum(abs.((K-ω²[4]*M)*v[4])./abs.(K*v[4])) < 1e-10
    @test maximum(abs.((K-ω²[5]*M)*v[5])./abs.(K*v[5])) < 2e-10
end




#end

# using GLMakie
# fig = Figure()
# display(fig) # open interactive window (gets closed down by "save")
# axe      = Axis(fig[1,1])

# with_theme(Theme(fontsize = 30,font="Arial")) do
#     for vᵢ∈v
#     oedge    = lines!(  axe,Real.(vᵢ) )
#     end
# end

