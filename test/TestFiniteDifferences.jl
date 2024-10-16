module TestFiniteDifferences



using Test,Muscade


n = 10
Δt = 0.1
t = (0:(n-1))*Δt 
x = 1 .+t .+0.5*t.^2 .+1/3*t.^3

x1  = [zeros(n),zeros(n),zeros(n)]
for order = 0:2
    for s = 1:n
        for (Δs,w) ∈ Muscade.finitediff(order,n,s)
            x1[order+1][s   ] += x[s+Δs] * w/Δt^order  
        end
    end
end

@testset "FinDiff" begin
    @test x1[1]≈x
    @test 300x1[2]≈ [316,334,373,418,469,526,589,658,733,772]
    @test 10x1[3]≈ [12,12,14,16,18,20,22,24,26,26]
end
end