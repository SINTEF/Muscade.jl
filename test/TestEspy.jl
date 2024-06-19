module TestEspy
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

exmerge = @macroexpand @espy function lagrangian(o::ElementConstraint{Teleobj,λinod,λfield,Nu}, Λ,X,U,A,t,χ,SP,dbg) where{Teleobj,λinod,λfield,Nu} 
    req        = merge(o.req)
    γ          = default{:γ}(SP,0.)
    u          = getsomedofs(U,SVector{Nu}(1:Nu)) 
    ☼λ         = ∂0(U)[Nu+1]
    L,χn,FB,☼eleres  = getlagrangian(implemented(o.eleobj)...,o.eleobj,Λ,X,u,A,t,χ,identity,SP,(dbg...,via=ElementConstraint),req)
    ☼gap       = o.gap(eleres,X,u,A,t,o.gargs...)
    if         o.mode(t)==:equal;    return L-gap*λ                  ,noχ,(α=∞                        ,)
    elseif     o.mode(t)==:positive; return L-KKT(λ,gap,γ,o.λₛ,o.gₛ)  ,noχ,(α=decided(λ/o.λₛ,gap/o.gₛ,γ),)
    elseif     o.mode(t)==:off;      return L-o.gₛ/(2o.λₛ)*λ^2        ,noχ,(α=∞                        ,)  
    end
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
    function residual(x::Vector{R}, y, req_001; ) where R <: Real
        out_001 = (;)
        ngp = 2
        req_001_gp = if haskey(req_001, :gp)
                req_001.gp
            else
                nothing
            end
        accum = ntuple(ngp) do igp
                out_001_gp = (;)
                z = x[igp] + y[igp]
                out_001_gp_001 = if haskey(req_001_gp, :z)
                        (out_001_gp..., z = z)
                    else
                        out_001_gp
                    end
                z
                if haskey(req_001_gp, :material)
                    (s, t, out_001_gp_002) = material(z, req_001_gp.material)
                    out_001_gp_002 = (out_001_gp_001..., material = out_001_gp_002)
                else
                    (s, t) = material(z)
                    out_001_gp_002 = out_001_gp_001
                end
                out_001_gp_003 = if haskey(req_001_gp, :s)
                        (out_001_gp_002..., s = s)
                    else
                        out_001_gp_002
                    end
                out_001_gp_004 = if haskey(req_001_gp, :t)
                        (out_001_gp_003..., t = t)
                    else
                        out_001_gp_003
                    end
                (s, t)
                (s = s, out = out_001_gp_004)
            end
        out_002 = if haskey(req_001, :gp)
                (out_001..., gp = NTuple{ngp}(((accum[igp]).out for igp = 1:ngp)))
            else
                out_001
            end
        r = sum((i->(accum[i]).s), ngp)
        return (r, nothing, nothing, out_002)
    end
end

exmaterial_ = quote
    function material(z)
        a = z + 1
        b = a * z
        return (b, 3.0)
    end
    function material(z, req_001; )
        out_001 = (;)
        a = z + 1
        out_002 = if haskey(req_001, :a)
                (out_001..., a = a)
            else
                out_001
            end
        a
        b = a * z
        out_003 = if haskey(req_001, :b)
                (out_002..., b = b)
            else
                out_002
            end
        b
        return (b, 3.0, out_003)
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
                nothing
                χ = randn()
                r = vcat(a, reshape(c, 4))
                (χ = χ, r = r)
            end
        χ = ntuple((i->(t[i]).χ), ngp)
        r = sum((i->(t[i]).r), ngp)
        return (r, χ, nothing)
    end
    function bar(x, y, z, req_001; )
        out_001 = (;)
        ngp = 4
        vec = SVector{2}
        p = 3
        out_002 = if haskey(req_001, :p)
                (out_001..., p = p)
            else
                out_001
            end
        p
        for i = 1:2
            j = i ^ 2
            k = i + j
        end
        req_001_gp = if haskey(req_001, :gp)
                req_001.gp
            else
                nothing
            end
        t = ntuple(ngp) do igp
                out_002_gp = (;)
                a = vec(igp, igp ^ 2)
                b = vec(1.0, 1.0)
                c = b * b'
                out_002_gp_001 = if haskey(req_001_gp, :c)
                        (out_002_gp..., c = c)
                    else
                        out_002_gp
                    end
                c
                if haskey(req_001_gp, :square)
                    square = c ^ 2
                    out_002_gp_002 = (out_002_gp_001..., square = square)
                else
                    out_002_gp_002 = out_002_gp_001
                end
                χ = randn()
                r = vcat(a, reshape(c, 4))
                (χ = χ, r = r, out = out_002_gp_002)
            end
        out_003 = if haskey(req_001, :gp)
                (out_002..., gp = NTuple{ngp}(((t[igp]).out for igp = 1:ngp)))
            else
                out_002
            end
        χ = ntuple((i->(t[i]).χ), ngp)
        r = sum((i->(t[i]).r), ngp)
        return (r, χ, nothing, out_003)
    end
end

exmerge_ = quote
    function lagrangian(o::ElementConstraint{Teleobj, λinod, λfield, Nu}, Λ, X, U, A, t, χ,  SP, dbg) where {Teleobj, λinod, λfield, Nu}
        req = merge(o.req)
        γ = default{:γ}(SP, 0.0)
        u = getsomedofs(U, SVector{Nu}(1:Nu))
        λ = (∂0(U))[Nu + 1]
        (L, χn, FB, eleres) = getlagrangian(implemented(o.eleobj)..., o.eleobj, Λ, X, u, A, t, χ,identity,  SP, (dbg..., via = ElementConstraint), req)
        gap = o.gap(eleres, X, u, A, t, o.gargs...)
        if o.mode(t) == :equal
            return (L - gap * λ, noχ, (α = ∞,))
        elseif o.mode(t) == :positive
            return (L - KKT(λ, gap, γ, o.λₛ, o.gₛ), noχ, (α = decided(λ / o.λₛ, gap / o.gₛ, γ),))
        elseif o.mode(t) == :off
            return (L - (o.gₛ / (2 * o.λₛ)) * λ ^ 2, noχ, (α = ∞,))
        end
    end
    function lagrangian(o::ElementConstraint{Teleobj, λinod, λfield, Nu}, Λ, X, U, A, t, χ,  SP, dbg, req_001; ) where {Teleobj, λinod, λfield, Nu}
        out_001 = (;)
        req = merge(req_001, o.req)
        γ = default{:γ}(SP, 0.0)
        u = getsomedofs(U, SVector{Nu}(1:Nu))
        λ = (∂0(U))[Nu + 1]
        out_002 = if haskey(req_001, :λ)
                (out_001..., λ = λ)
            else
                out_001
            end
        λ
        (L, χn, FB, eleres) = getlagrangian(implemented(o.eleobj)..., o.eleobj, Λ, X, u, A, t, χ, identity,  SP, (dbg..., via = ElementConstraint), req)
        out_003 = if haskey(req_001, :eleres)
                (out_002..., eleres = eleres)
            else
                out_002
            end
        (L, χn, FB, eleres)
        gap = o.gap(eleres, X, u, A, t, o.gargs...)
        out_004 = if haskey(req_001, :gap)
                (out_003..., gap = gap)
            else
                out_003
            end
        gap
        if o.mode(t) == :equal
            return (L - gap * λ, noχ, (α = ∞,), out_004)
        elseif o.mode(t) == :positive
            return (L - KKT(λ, gap, γ, o.λₛ, o.gₛ), noχ, (α = decided(λ / o.λₛ, gap / o.gₛ, γ),), out_004)
        elseif o.mode(t) == :off
            return (L - (o.gₛ / (2 * o.λₛ)) * λ ^ 2, noχ, (α = ∞,), out_004)
        end
    end
end

@testset "Transformed codes" begin
    @test prettify(exresidual)   == prettify(exresidual_)
    @test prettify(exmaterial)   == prettify(exmaterial_)
    @test prettify(exlagrangian) == prettify(exlagrangian_)
    @test prettify(exbar)        == prettify(exbar_)
    @test prettify(exmerge)        == prettify(exmerge_)
end

r1 = @request (a,b,c)
r2 = @request (a,b,c(x,y))  
@testset "Request" begin
    @test r1 == (a = nothing, b = nothing, c = nothing)
    @test r2 == (a = nothing, b = nothing, c = (x = nothing, y = nothing))
end

r1 = @request a(x(α),y),b,e
r2 = @request a(x(α),z),c(x,y),d(f),e
r = merge(r1,r2)
@testset "MergeRequest" begin
    @test r == (a=(x=(α=nothing,),y=nothing,z=nothing),b=nothing,e=nothing,c=(x=nothing,y=nothing),d=(f=nothing,))
end


eval(exresidual)
eval(exmaterial)
eval(exlagrangian)
eval(exbar)

req         = @request gp(z,s,t,material(a,b))
r,χ,SFB,out = residual([1.,2.],[3.,4.],req)

@testset "Get output" begin
    @test out ==  (gp = ((z = 4.0, material = (a = 5.0, b = 20.0), s = 20.0, t = 3.0), (z = 6.0, material = (a = 7.0, b = 42.0), s = 42.0, t = 3.0)),)
end


end # Module
