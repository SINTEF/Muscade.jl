module TestEspy2
using Muscade,Test#,StaticArrays

#=
Code for a material model and a finite element
The following code is part of a model for the elongation
of a rod hanging from an attachement, under its own weight.
Parts of the element code that are not relevant to EspyInsideFubction 
are left out.

The function "force" takes in nodal displacements, and returns nodal forces,
this being what the finite element solution requires.  But we are interested
in reporting intermediate results from inisde this function (and the 
function "material" that it calls).
=#
struct Material
    E   :: Float64 # Young's modulus
end
@espy function material(m::Material,ε)
    :σ = m.E*ε
    return σ
end
requestable(m::Material) = (σ=scalar,)


struct Element  
    L₀  :: Float64  # underformed length
    A   :: Float64  # cross section area
    ρg  :: Float64  # density times gravity
    mat :: Material 
end
const ngp = 1 # 1 Gauss quadrature point
@espy function force(e::Element,ΔX)
    :w      = e.ρg*e.A*e.L₀
    :R      = [w/2,w/2]
    for igp = 1:ngp # loop over 1 point, but still a loop
        :ε  = (ΔX[2]-ΔX[1])/e.L₀
        σ   = :material(e.mat,ε)  
        :T  = e.A*σ
        R   = +[T,-T]
    end
    return R
end
requestable(e::Element) = (w=scalar, R=(2,), gp=forloop(ngp, (ε=scalar, T=scalar, material=requestable(e.mat)) ) )


#=
Now we take these elements into use to solve ten elements hanging in a 
straight vertical line. 
(We skip the actual solution process, which is not relevant to this demo: 
the value of ΔX is hard-coded).
=#
nel  = 10
topo = [[i,i+1] for i = 1:nel]
m    = Material(2.1e11)
e    = Element(1.,1e-4,8000*9.81,m)
ΔX   = [1/2*e.ρg/m.E*((iel*e.L₀)^2-(nel*e.L₀)^2) for iel = 0:nel]  # say we found a solution,

#=
The analysis has been completed, and ΔX has been stored: with it we
can compute any intermediate result.

The user defines which results are wanted.
=# 

request   = @request w,gp[].(ε,material.(σ))

#= 
The finite element software provides a function
which includes the following code, and returns
`key` and `out`
=#    
key,nkey  = makekey(request,requestable(e))

out = Matrix{Float64}(undef,nkey,nel)
for (iel,t) ∈ enumerate(topo)
    _ = force(@view(out[:,iel]),key, e,ΔX[t])
end

#=
The user can now access the intermediate results
=#
iel = 4
igp = 1
σ = out[key.gp[igp].material.σ,iel]
ε = out[key.gp[igp].ε         ,iel]

#
@testset "Results" begin
    @test request == :((w, (gp[]).(ε, material.(σ))))
    @test σ ≈ 274680.
    @test ε ≈ 1.308e-6
    @test key == (w = 1, gp = [(ε = 2, material = (σ = 3,))])
    @test out ≈ [7.848     7.848       7.848       7.848     7.848       7.848       7.848       7.848       7.848       7.848;
                1.86857e-7 5.60571e-7  9.34286e-7  1.308e-6  1.68171e-6  2.05543e-6  2.42914e-6  2.80286e-6  3.17657e-6  3.55029e-6;
                39240.0    117720.0    196200.0    274680.0  353160.0    431640.0    510120.0    588600.0    667080.0    745560.0]
end    
end