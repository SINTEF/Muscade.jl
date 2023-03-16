using Test,StaticArrays
using Muscade
## Test request

struct Element  end
requestable(el::Element)= (X=(3,2),  gp=forloop(2, (F=(2,2), material=(ε=(2,2),Σ=(2,2)))) )
el = Element()
r8 = @request X,gp[].(F,material.(Σ,ε))

ra    = requestable(el)

## annotated code (written by user)

struct El  end
requestable(el::El)= (gp=forloop(2, (z=scalar,s=scalar, material=(a=scalar,b=scalar))),)
el = El()

@espydbg function residual(x::R,y) where{R<:Real}
    ngp=2
    r = 0
    for igp=1:ngp
        ☼z = x[igp]+y[igp]
        ☼s,☼t  = ☼material(z)
        r += s
    end
    return r
end
@espydbg function material(z)
    ☼a = z+1
    ☼b = a*z
    return b,3.
end
@espydbg lagrangian(o::Vector{R}, δX,X,U,A, t,γ,dbg) where{R} = ☼J = o.cost(∂n(X,R)[1],t)

@espydbg function foo()
    ngp = 4
    SV2 = SVector{2} 
    ☼r   = sum(1:ngp) do igp
        a = SV2(igp,igp^2)
        ☼b = SV2(1.,1.)
        c = b*b'
        vcat(a,reshape(c,4))
    end
end

@espydbg function bar()
    ngp = 4
    vec = SVector{2}
    t = ntuple(ngp) do igp
        a = vec(igp,igp^2)
        b = vec(1.,1.)
        ☼c = b*b'
        ♢square =
        χ = randn()
        r = vcat(a,reshape(c,4))
        (χ,r)
    end
    χ = ntuple(i->t[i][1],ngp)
    r = sum(   i->t[i][2],ngp)
    return r,χ
end


;
