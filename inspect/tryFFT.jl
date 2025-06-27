using Muscade

pref   = 6
nref   = 2^pref
δωref  = 0.1

Δp     = 0  # changing this term to positive values leaves ampltudes and energies unchanged
p      = pref+Δp
n      = 2^p
δω     = δωref/2^Δp 

nc     = 2^(p-1)
δt     = Muscade.getδt(n,δω)
t      = range(start=0.,length=n ,step=δt)
ω      = range(start=0.,length=nc,step=δω)

X      = zeros(Complex{Float64},nc,1)
X[1+2^Δp,1] = 10. * 2^Δp


Z      = copy(X)
Muscade.𝔉⁻¹!(Z,δω)
Muscade.𝔉!(Z,δt)
@show maximum(abs.(X.-Z))


Y      = copy(X)
x      = reinterpret(Float64,Y)
Muscade.𝔉⁻¹!(Y,δω)
    

# signal power - Check Parseval and invariance with δt
@show sum(abs2.(x))*δt/2 / (n*δt)
@show sum(abs2.(X))*δω   / (n*δt)

@show (n+1)*δt /2^Δp
using GLMakie
fig      = Figure(size = (1000,600))
display(fig) # open interactive window (gets closed down by "save")

axe_time      = Axis(fig[1,1],xlabel="t",ylabel="x")
axe_freq      = Axis(fig[2,1],xlabel="ω",ylabel="X")
i = 1
with_theme(Theme(fontsize = 30,font="Arial")) do
   oedge    = lines!(  axe_time,t,x[:,i]          ,color = :black,linewidth = 1)
end
with_theme(Theme(fontsize = 30,font="Arial")) do
    oedge    = lines!(  axe_freq,ω,  real(X[:,i])      ,color = :black,linewidth = 1)
    oedge    = lines!(  axe_freq,ω,  imag(X[:,i])      ,color = :red,linewidth = 1)
end

