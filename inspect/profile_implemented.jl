using Profile,ProfileView,BenchmarkTools
using Muscade
using Muscade.ElTest
using StaticArrays,Printf

include("..\\test\\SomeElements.jl")

### Turbine
sea(t,x) = SVector(1.,0.)
sky(t,x) = SVector(0.,1.)
turbine  = Turbine(SVector(0.,0.),-10., 2.,sea, 3.,sky)
Λ        = @SVector [1.,1.]
X        = @SVector [1.,2.]
U        = @SVector 𝕣[]
A        = @SVector [0.,0.]  # [Δseadrag,Δskydrag]

#                            eleobj, Λ, X,  U,  A, t, χ,      χcv,     SP,     dbg
L,Lλ,Lx,Lu,La,χn  = gradient(turbine,Λ ,[X],[U],A, 0.,nothing,identity,nothing,(;))

Profile.clear()
Profile.@profile for i=1:1000000
    local L,Lλ,Lx,Lu,La,χn  = gradient(turbine,Λ ,[X],[U],A, 0.,nothing,identity,nothing,(;))
end
ProfileView.view(fontsize=30);


;