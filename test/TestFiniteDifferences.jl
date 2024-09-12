module TestFiniteDifferences



using Test,Muscade


n = 10
Δt = 0.1
t = (0:(n-1))*Δt 
x = 1 .+t .+0.5*t.^2 .+1/3*t.^3

x1  = [zeros(n),zeros(n),zeros(n)]
x2  = [zeros(n),zeros(n),zeros(n)]
for order = 0:2
    for s = 1:n
        for (Δs,w) ∈ Muscade.finitediff(order,n,s)
            x1[order+1][s   ] += x[s+Δs] * w/Δt^order  
        end
        for (Δs,w) ∈ Muscade.finitediff(order,n,s,transposed=true)
            x2[order+1][s+Δs] += x[s   ] * w/Δt^order  
        end
    end
end

@testset "FinDiff" begin
    @test x1[1]≈x
    @test 300x1[2]≈ [316,334,373,418,469,526,589,658,733,772]
    @test 10x1[3]≈ [12,12,14,16,18,20,22,24,26,26]
    @test x1[1]≈x2[1]
    @test x1[2]≈x2[2]
    @test x1[3]≈x2[3]
end

# using GLMakie
# fig      = Figure(size = (500,500))
# display(fig) # open interactive window (gets closed down by "save")
# axe      = Axis(fig[1,1],title="FDsparsity",xlabel="i",ylabel="j")
# (i,j) = Muscade.FDsparsity(2,7)
# scatter!(axe,i,-j)



@testset "FDsparsity" begin
    @test Muscade.FDsparsity(0,5) == ([1,2,3,4,5],[1,2,3,4,5])
    @test Muscade.FDsparsity(1,5) == ([1,2,3,4,5,1,2,3,3,4,5],[1,2,3,4,5,3,4,5,1,2,3])
    @test Muscade.FDsparsity(2,5) == ([1,2,3,4,5,1,2,3,3,4,5,1,2,3,4,2,3,4,5],[1,2,3,4,5,3,4,5,1,2,3,2,3,4,5,1,2,3,4])
end


end