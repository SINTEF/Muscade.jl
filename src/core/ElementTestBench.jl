using  Printf,Zygote,StaticArrays
#using Muscade.Tools.Dialect


@Zygote.adjoint (T::Type{<:SVector})(x::Number...     ) = T(x...), dv -> (nothing, dv...)
@Zygote.adjoint (T::Type{<:SVector})(x::AbstractVector) = T(x   ), dv -> (nothing, dv   )


"""
    nodes = nodesforelementtest(coord)

Create a vector of `Node` objects at specified coordinates, where
- `coord[inod,ix]` is the `ix`-th coordinate of the `inod`-th node
- `nodes[inod]` is a `Node`
"""
function nodesforelementtest(coord)
    model = Model(field=(), fieldscale=(), residualscale=())
    inods = [model(Node,coord[inod,:]) for inod ∈ axes(coord,1)]
    return model.nod[inods]
end

function showout(key::NamedTuple,out::𝕣1,str)
    for k ∈ keys(key)
        showout(key[k],out,str*"."*string(k))
    end
end
function showout(key::Array,out::𝕣1,str)
    for k ∈ eachindex(key)
        showout(key[k],out,str*"["*string(k)*"]")
    end
end
function showout(key::𝕫1,out::𝕣1,str)
    @printf "\n%s = \n" str
    for k∈key
        @printf "%10.3g\n" out[k]
    end
end
function showout(key::𝕫,out::𝕣1,str)
    @printf "\n%s = %10.3g\n" str out[key]
end
function showout(key::𝕫2,out::𝕣1,str)
    @printf "\n%s = \n" str
    for i∈axes(key,1)
        for j∈axes(key,2)
            @printf "%10.3g  " out[key[i,j]]
        end
        @printf "\n"
    end
end
function showχ(χ::NTuple,str)
    for x∈χ
        showχ(x)
    end
end
function showχ(χ::NamedTuple,str)
    for k ∈ keys(χ)
        showχ(χ[k],str*"."*string(k))
    end
end
function showχ(χ::Array,str)
    for k ∈ eachindex(χ)
        showχ(χ[k],str*"["*string(k)*"]")
    end
end
function showχ(χ::𝕣,str)
    @printf "\n%s = %10.3g \n" str χ
end
function showχ(χ::AbstractVector{𝕣},str)
    @printf "\n%s = \n" str
    for i∈χ
        @printf "%10.3g\n" χ[i]
    end
end
function showχ(χ::AbstractMatrix{𝕣},str)
    @printf "\n%s = \n" str
    for i∈axes(χ,1)
        for j∈axes(χ,2)
            @printf "%10.3g  " χ[i,j]
        end
        @printf "\n"
    end
end

gr(x,∇x::Nothing)= x .*0
gr(x,∇x         )=∇x


function testStaticElement(el; δX,X,U,A, χo=initstate(el), t::𝕣=0.,ε::𝕣=0., req=nothing,verbose::𝕓=true)
    id       = dofid(el)
    n        = neldof(el) 
    χcv      = identity
    dbg      = ()
    L,χn = lagrangian(el, [δX],[X],[U],A, χo,χcv, t,ε,dbg)
    function closure(δX,X,U,A)
        L,χn = lagrangian(el, [δX],[X],[U],A, χo,χcv, t,ε,dbg)
        return L
    end
    Lδx,Lx,Lu,La = gradient(closure,δX,X,U,A)
    Lδx,Lx,Lu,La = gr(δX,Lδx),gr(X,Lx),gr(U,Lu),gr(A,La)
    if verbose
        @printf "\nElement type: %s\n" typeof(el)
        if n.X > 0
            @printf "\n    idof               doftyp   inod          δX           X         Lδx          Lx \n"
            for idof = 1:n.X
                @printf "    %4d     %16s  %5d  %10.3g  %10.3g  %10.3g  %10.3g\n" idof id.X.typ[idof] id.X.nod[idof] δX[idof] X[idof] Lδx[idof] Lx[idof]
            end
        end
        if n.U > 0
            @printf "\n    idof               doftyp   inod           U          Lu \n"
            for idof = 1:n.U
                @printf "    %4d     %16s  %5d  %10.3g  %10.3g\n" idof id.U.typ[idof] id.U.nod[idof] U[idof] Lu[idof]
            end
        end
        if n.A > 0
            @printf "\n    idof               doftyp   inod           A          La \n"
            for idof = 1:n.A
                @printf "    %4d     %16s  %5d  %10.3g  %10.3g\n" idof id.A.typ[idof] id.A.nod[idof] A[idof] La[idof]
            end
        end
        showχ(χo,"χo")
        showχ(χn,"χn")
    end

    return Lδx,Lx,Lu,La,χo,χn
    # if isnothing(req)
    #     L,χn,Lδx,Lx,Lu,La = lagrangian(el, [δX],[X],[U],A, χo,χcv, t,ε,dbg)
    #     return L,χo,χn
    # else
    #     key,nkey = makekey(req,requestable(el))
    #     out      = 𝕣1(undef,nkey)
    #     L,χn     = lagrangian(out,key,el, [δX],[X],[U],A, χo,χcv, t,ε,dbg)
    #     return L,χo,χn,out,key
    # end
end


