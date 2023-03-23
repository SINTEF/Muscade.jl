module TestEspy2
using Test,StaticArrays
using Muscade,MacroTools

macro expression(ex) return QuoteNode(ex) end
clean(ex)   = prettify(Muscade.code_clean_function(ex))
extract(ex) = prettify(Muscade.code_extractor_function(ex,:out,:key,false))
pretty(ex)  = println(prettify(ex))
####

ex1 = @expression function residual(x::R,y) where{R<:Real}
    ngp=2
    r = 0
    for igp=1:ngp
        ☼z = x[igp]+y[igp]
        ☼s,☼t  = ☼material(z)
        r += s
    end
    return r
end
ex1c = @expression function residual(x::R,y) where{R<:Real}
    ngp=2
    r = 0
    for igp=1:ngp
        z = x[igp]+y[igp]
        s,t  = material(z)
        r += s
    end
    return r
end
ex1e = @expression function residual(out, key, x::R, y; ) where R <: Real
    ngp = 2
    r = 0
    for igp = 1:ngp
        if haskey(key, :gp)
            key_gp = key.gp[igp]
        else
            key_gp = Nothing
        end
        z = x[igp] + y[igp]
        if haskey(key_gp, :z)
            if typeof(key_gp.z) == Int64
                out[key_gp.z] = z
            else
                out[key_gp.z] .= z
            end
        end
        z
        if haskey(key_gp, :material)
            (s, t) = material(out, key_gp.material, z)
        else
            (s, t) = material(z)
        end
        if haskey(key_gp, :s)
            if typeof(key_gp.s) == Int64
                out[key_gp.s] = s
            else
                out[key_gp.s] .= s
            end
        end
        if haskey(key_gp, :t)
            if typeof(key_gp.t) == Int64
                out[key_gp.t] = t
            else
                out[key_gp.t] .= t
            end
        end
        (s, t)
        r += s
    end
    return r
end
#######
ex2 = @expression function material(z)
    ☼a = z+1
    ☼b = a*z
    return b,3.
end
ex2c = @expression function material(z)
    a = z + 1
    b = a * z
    return (b, 3.0)
end
ex2e = @expression function material(out, key, z; )
    a = z + 1
    if haskey(key, :a)
        if typeof(key.a) == Int64
            out[key.a] = a
        else
            out[key.a] .= a
        end
    end
    a
    b = a * z
    if haskey(key, :b)
        if typeof(key.b) == Int64
            out[key.b] = b
        else
            out[key.b] .= b
        end
    end
    b
    return (b, 3.0)
end
###
ex3 = @expression lagrangian(o::Vector{R}, δX,X,U,A, t,γ,dbg) where{R} = ☼J = o.cost(∂n(X,R)[1],t)
ex3c = @expression (lagrangian(o::Vector{R}, δX, X, U, A, t, γ, dbg) where R) = (J = o.cost((∂n(X, R))[1], t))
ex3e = @expression function lagrangian(out, key, o::Vector{R}, δX, X, U, A, t, γ, dbg; ) where R
    J = o.cost((∂n(X, R))[1], t)
    if haskey(key, :J)
        if typeof(key.J) == Int64
            out[key.J] = J
        else
            out[key.J] .= J
        end
    end
    J
end
#####
ex4 = @expression function foo()
    ngp = 4
    SV2 = SVector{2} 
    ☼r   = sum(1:ngp) do igp
        a = SV2(igp,igp^2)
        ☼b = SV2(1.,1.)
        c = b*b'
        vcat(a,reshape(c,4))
    end
end
ex4c = @expression function foo()
    ngp = 4
    SV2 = SVector{2}
    r = sum(1:ngp) do igp        
            a = SV2(igp, igp ^ 2)
            b = SV2(1.0, 1.0)    
            c = b * b'
            vcat(a, reshape(c, 4))
        end
end
ex4e = @expression function foo(out, key; )
    ngp = 4
    SV2 = SVector{2}    
    r = sum(1:ngp) do igp
            if haskey(key, :gp)
                key_gp = key.gp[igp]
            else
                key_gp = Nothing
            end
            a = SV2(igp, igp ^ 2)
            b = SV2(1.0, 1.0)
            if haskey(key_gp, :b)
                if typeof(key_gp.b) == Int64
                    out[key_gp.b] = b
                else
                    out[key_gp.b] .= b
                end
            end
            b
            c = b * b'
            vcat(a, reshape(c, 4))
        end
    if haskey(key, :r)
        if typeof(key.r) == Int64
            out[key.r] = r
        else
            out[key.r] .= r
        end
    end
    r
end
#####
ex5 = @expression function bar()
    ngp = 4
    vec = SVector{2}
    t = ntuple(ngp) do igp
        a = vec(igp,igp^2)
        b = vec(1.,1.)
        ☼c = b*b'
        ♢square = c^2
        χ = randn()
        r = vcat(a,reshape(c,4))
        (χ,r)
    end
    χ = ntuple(i->t[i][1],ngp)
    r = sum(   i->t[i][2],ngp)
    return r,χ
end
ex5c = @expression function bar()
    ngp = 4
    vec = SVector{2}
    t = ntuple(ngp) do igp
            a = vec(igp, igp ^ 2)
            b = vec(1.0, 1.0)
            c = b * b'
            square = c ^ 2
            χ = randn()
            r = vcat(a, reshape(c, 4))
            (χ, r)
        end
    χ = ntuple((i->(t[i])[1]), ngp)
    r = sum((i->(t[i])[2]), ngp)
    return (r, χ)
end
ex5e = @expression function bar(out, key; )
    ngp = 4
    vec = SVector{2}
    t = ntuple(ngp) do igp
            if haskey(key, :gp)
                key_gp = key.gp[igp]
            else
                key_gp = Nothing
            end
            a = vec(igp, igp ^ 2)
            b = vec(1.0, 1.0)
            c = b * b'
            if haskey(key_gp, :c)
                if typeof(key_gp.c) == Int64
                    out[key_gp.c] = c
                else
                    out[key_gp.c] .= c
                end
            end
            c
            if haskey(key_gp, :square)
                square = c ^ 2
                if typeof(key_gp.square) == Int64
                    out[key_gp.square] = square
                else
                    out[key_gp.square] .= square
                end
            end
            χ = randn()
            r = vcat(a, reshape(c, 4))
            (χ, r)
        end
    χ = ntuple((i->(t[i])[1]), ngp)
    r = sum((i->(t[i])[2]), ngp)
    return (r, χ)
end

@testset "Transformed codes" begin
    @test clean(ex1)   == prettify(ex1c)
    @test extract(ex1) == prettify(ex1e)
    @test clean(ex2)   == prettify(ex2c)
    @test extract(ex2) == prettify(ex2e)
    @test clean(ex3)   == prettify(ex3c)
    @test extract(ex3) == prettify(ex3e)
    @test clean(ex4)   == prettify(ex4c)
    @test extract(ex4) == prettify(ex4e)
    @test clean(ex5)   == prettify(ex5c)
    @test extract(ex5) == prettify(ex5e)
end
end


;
