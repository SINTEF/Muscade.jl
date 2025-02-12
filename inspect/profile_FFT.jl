# not really profiling - rather analysing type stability
using Profile,ProfileView,BenchmarkTools
using Muscade

g(t)  = 1/√(2π)*exp(-1/2*t^2)   
G(ω)  = 1/√(2π)*exp(-1/2*ω^2)  

# G is obtained from the unitary, angular frequency Fourrier transform
# 𝔉, operates only on real ts, which have a Hermitian transform.  So 𝔉 returns only for ω>0, multiplied by 2 for unitariness

nstep = 2^16
ndof  = 1000

X       = Matrix{Complex{Float64}}(undef,div(nstep,2),ndof)
x       = reinterpret(Float64,X)
x       = randn(nstep,ndof)
Muscade.𝔉!(X,0.1)


#@code_warntype Muscade.assemble!(out,asm,dis,model1,state0,(;)) # works - but only on top function

mission = :profile
if  mission == :time
    Muscade.𝔉!(X,0.1)
    @btime Muscade.𝔉!(X,0.1) # Matlab 261 ms

elseif mission == :profile
    Profile.clear()
    Profile.@profile for i=1:10
        Muscade.𝔉!(X,0.1)
    end
    ProfileView.view(fontsize=30);
end
# After clicking on a bar in the flame diagram, you can type warntype_last() and see the result of 
# code_warntype for the call represented by that bar.



;