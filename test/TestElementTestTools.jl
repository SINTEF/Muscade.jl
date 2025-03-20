module TestElementTestTools

using Test, Muscade, StaticArrays

model           = Model(:TestModel)
node1           = addnode!(model,𝕣[0,0,0])
#elnod           = model.nod

uload = SingleUdof(model.nod;Xfield=:Xdof,Ufield=:Udof,cost=(u,t)->.5*u^2)
rload = DofLoad(model.nod;field=:Xdof,value=t->cos(t))

Λ = SVector(0.)
X = (SVector(0.1),)
U = (SVector(0.3),)
A = SVector{0,Float64}()
out1 = diffed_lagrangian(uload; Λ,X,U,A)
# @printf "\nU\n"
# print_element_array(uload,:U,U[1])
# @printf "∂L/∂U₀\n"
# print_element_array(uload,:U,out1.∇L[3][1])
# @printf "\n∂²L/∂U₀²\n"
# print_element_array(uload,:U,out1.HL[3,3][1,1])

@testset "diffed_lagrangian" begin
    @test out1.∇L[3][1]   ≈ [.3]
    @test out1.HL[3,3][1,1] ≈ [1. ;;]
end


X = (SVector(0.1),)
U = (SVector{0,Float64}(),)
A = SVector{0,Float64}()
out2 = diffed_residual(rload; X,U,A)
# @printf "\nR\n"
# print_element_array(uload,:X,out2.R)
# @printf "\n∂R/∂X₀\n"
# print_element_array(uload,:X,out2.∇R[2][1])

@testset "diffed_residual" begin
    @test out2.R   ≈ [-1.]
    @test out2.∇R[2][1] ≈ [0.;;]
end


end

