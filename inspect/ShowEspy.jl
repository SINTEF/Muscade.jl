using Test,StaticArrays,MacroTools
include("../src/Dialect.jl")
include("../src/Espy.jl")
include("../src/Exceptions.jl")
#using Muscade

# @espydbg function residual(x::Vector{R},y) where{R<:Real}
#     ngp = 2
#     α,β = 1,2
#     va = 1,2,3
#     tu = (a=1.,b=2)
#     e,f = foo()
#     accum = ntuple(ngp) do igp
#         ☼z = x[igp]+y[igp]
#         ☼s,☼t  = ☼material(z)
#         @named(s) 
#     end
#     r = sum(i->accum[i].s,ngp)
#     return r,nothing,nothing
# end
# @espydbg function material(z)
#     ☼a = z+1
#     ☼b = a*z
#     return b,3.
# end

# @espydbg lagrangian(o::Vector{R}, δX,X,U,A, t,γ,dbg) where{R} = o.cost(∂n(X,R)[1],t),nothing,nothing

exbar= @espydbg function bar(x,y,z)
    ngp = 4
    vec = SVector{2}
    ☼p = 3
    for i ∈ 1:2
        j = i^2
        k = i+j
    end
    t = ntuple(ngp) do igp
        a = vec(igp,igp^2)
        b = vec(1.,1.)
        ☼c = b*b'
        ♢square = c^2
        χ = randn()
        r = vcat(a,reshape(c,4))
        @named(χ,r)
    end
    χ = ntuple(i->t[i].χ,ngp)
    r = sum(   i->t[i].r,ngp)
    return r,χ,nothing
end


;
