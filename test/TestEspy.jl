#module TestEspy
using Test,StaticArrays,MacroTools
include("../src/Dialect.jl")
include("../src/Espy.jl")
include("../src/Exceptions.jl")

macro expression(ex) return QuoteNode(ex) end
pretty(ex)  = println(prettify(ex))
####



exresidual = @macroexpand @espy function residual(x::Vector{R},y) where{R<:Real}
    ngp=2
    accum = ntuple(ngp) do igp
        ☼z = x[igp]+y[igp]
        ☼s,☼t  = ☼material(z)
        @named(s) 
    end
    r = sum(i->accum[i].s,ngp)
    return r,nothing,nothing
end
exmaterial = @macroexpand @espy function material(z)
    ☼a = z+1
    ☼b = a*z
    return b,3.
end

exlagrangian = @macroexpand @espy lagrangian(o::Vector{R}, δX,X,U,A, t,γ,dbg) where{R} = o.cost(∂n(X,R)[1],t),nothing,nothing

exbar = @macroexpand @espy function bar(x,y,z)
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

############# @espy outputs

exresidual_ = quote
    function residual(x::Vector{R}, y) where R <: Real
        ngp = 2
        accum = ntuple(ngp) do igp
                z = x[igp] + y[igp]
                (s, t) = material(z)
                (s = s,)
            end
        r = sum((i->(accum[i]).s), ngp)
        return (r, nothing, nothing)
    end
    function residual(x::Vector{R}, y, req_002; ) where R <: Real
        out_001 = (;)
        ngp = 2
        req_002_gp = if haskey(req_002, :gp)
                req_002.gp
            else
                nothing
            end
        accum = ntuple(ngp) do igp
                out_001_gp = (;)
                z = x[igp] + y[igp]
                out_001_gp_004 = if haskey(req_002_gp, :z)
                        (out_001_gp..., z = z)
                    else
                        out_001_gp
                    end
                z
                if haskey(req_002_gp, :material)
                    (s, t, out_001_gp_006) = material(z, req_002_gp.material)
                    out_001_gp_005 = (out_001_gp_004..., material = out_001_gp_006)
                else
                    (s, t) = material(z)
                    out_001_gp_005 = out_001_gp_004
                end
                (s, t)
                (s = s, out = out_001_gp_005)
            end
        out_003 = if haskey(req_002, :gp)
                (out_001..., gp = NTuple{ngp}(((accum[igp]).out for igp = 1:ngp)))
            else
                out_001
            end
        r = sum((i->(accum[i]).s), ngp)
        return (r, nothing, nothing, out_003)
    end
end
exmaterial_ = quote
    function material(z)
        a = z + 1
        b = a * z
        return (b, 3.0)
    end
    function material(z, req_008; )
        out_007 = (;)
        a = z + 1
        out_009 = if haskey(req_008, :a)
                (out_007..., a = a)
            else
                out_007
            end
        a
        b = a * z
        out_010 = if haskey(req_008, :b)
                (out_009..., b = b)
            else
                out_009
            end
        b
        return (b, 3.0, out_010)
    end
end

exlagrangian_ = quote
    (lagrangian(o::Vector{R}, δX, X, U, A, t, γ, dbg) where R) = (o.cost((∂n(X, R))[1], t), nothing, nothing)
    (lagrangian(o::Vector{R}, δX, X, U, A, t, γ, dbg, req) where R) = (o.cost((∂n(X, R))[1], t), nothing, nothing, nothing)
end

exbar_ = quote
    function bar(x, y, z)
        ngp = 4
        vec = SVector{2}
        p = 3
        for i = 1:2
            j = i ^ 2
            k = i + j
        end
        t = ntuple(ngp) do igp
                a = vec(igp, igp ^ 2)
                b = vec(1.0, 1.0)
                c = b * b'
                square = c ^ 2
                χ = randn()
                r = vcat(a, reshape(c, 4))
                (χ = χ, r = r)
            end
        χ = ntuple((i->(t[i]).χ), ngp)
        r = sum((i->(t[i]).r), ngp)
        return (r, χ, nothing)
    end
    function bar(x, y, z, req_014; )
        out_013 = (;)
        ngp = 4
        vec = SVector{2}
        p = 3
        out_015 = if haskey(req_014, :p)
                (out_013..., p = p)
            else
                out_013
            end
        p
        for i = 1:2
            j = i ^ 2
            k = i + j
        end
        req_014_gp = if haskey(req_014, :gp)
                req_014.gp
            else
                nothing
            end
        t = ntuple(ngp) do igp
                out_015_gp = (;)
                a = vec(igp, igp ^ 2)
                b = vec(1.0, 1.0)
                c = b * b'
                out_015_gp_017 = if haskey(req_014_gp, :c)
                        (out_015_gp..., c = c)
                    else
                        out_015_gp
                    end
                c
                if haskey(req_014_gp, :square)
                    square = c ^ 2
                    out_015_gp_018 = (out_015_gp_017..., square = square)
                else
                    out_015_gp_018 = out_015_gp_017
                end
                χ = randn()
                r = vcat(a, reshape(c, 4))
                (χ = χ, r = r, out = out_015_gp_018)
            end
        out_016 = if haskey(req_014, :gp)
                (out_015..., gp = NTuple{ngp}(((t[igp]).out for igp = 1:ngp)))
            else
                out_015
            end
        χ = ntuple((i->(t[i]).χ), ngp)
        r = sum((i->(t[i]).r), ngp)
        return (r, χ, nothing, out_016)
    end
end

@testset "Transformed codes" begin
    @test prettify(exresidual)   == prettify(exresidual_)
    @test prettify(exmaterial)   == prettify(exmaterial_)
    @test prettify(exlagrangian) == prettify(exlagrangian_)
    @test prettify(exbar)        == prettify(exbar_)
end

eval(exresidual)
eval(exmaterial)
eval(exlagrangian)
eval(exbar)

req         = @request gp(z,s,t,material(a,b))
r,χ,SFB,out = residual([1.,2.],[3.,4.],req)

@testset "Get output" begin
    @test out == (gp = ((z = 4.0, material = (a = 5.0, b = 20.0)), (z = 6.0, material = (a = 7.0, b = 42.0))),)
end


#end # Module
