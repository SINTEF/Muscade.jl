"""
    inftyp,rettyp = @typeof(foo(args...[;kwargs...]))

    Determine the inferred type and the returned type of the output[s] returned by the relevant method-instance of foo.
    
"""
macro typeof(ex)
    _inferred_type(ex, __module__)
end
function _inferred_type(ex, mod)
    if Meta.isexpr(ex, :ref)
        ex = Expr(:call, :getindex, ex.args...)
    end
    Meta.isexpr(ex, :call)|| error("@inferred requires a call expression")
    farg = ex.args[1]
    if isa(farg, Symbol) && farg !== :.. && first(string(farg)) == '.'
        farg = Symbol(string(farg)[2:end])
        ex = Expr(:call, GlobalRef(Test, :_materialize_broadcasted),
            farg, ex.args[2:end]...)
    end
    result = let ex = ex
        quote
            $(if any(@nospecialize(a)->(Meta.isexpr(a, :kw) || Meta.isexpr(a, :parameters)), ex.args)
                # Has keywords
                args   = gensym()
                kwargs = gensym()
                quote
                    $(esc(args)), $(esc(kwargs)), result = $(esc(Expr(:call, _args_and_call, ex.args[2:end]..., ex.args[1])))
                    inftype = $(gen_call_with_extracted_types(mod, Base.infer_return_type, :($(ex.args[1])($(args)...; $(kwargs)...)); is_source_reflection = false))
                end
            else
                # No keywords
                quote
                    args    = ($([esc(ex.args[i]) for i = 2:length(ex.args)]...),)
                    result  = $(esc(ex.args[1]))(args...)
                    inftype = Base.infer_return_type($(esc(ex.args[1])), Base.typesof(args...))
                end
            end)
            rettype = result isa Type ? Type{result} : typeof(result)
            (inftype,rettype)
        end
    end
    return result
end


"""
    print_element_array(eleobj,class,V)

Show a vector (or a matrix) `V`, the rows of `V` being described as corresponding to `eleobj` dof of class `class`.
This can be used to print degrees of freedom, residuals, their derivatives, or gradients and Hessian of the Lagrangian.

See also: [`diffed_residual`](@ref), [`diffed_lagrangian`](@ref)
"""    
print_element_array(ele::AbstractElement,class::Symbol,V::AbstractVector) = print_element_array(ele,class,reshape(V,(length(V),1)))
function print_element_array(ele::Eletyp,class::Symbol,V::AbstractMatrix) where{Eletyp<:AbstractElement}
    inod,~,field     = Muscade.getdoflist(Eletyp)
    iVdof            = Muscade.getidof(Eletyp,class)
    (nV,ncol)      = size(V)
    @assert nV==Muscade.getndof(Eletyp,class)
    @printf "    i  ieldof               doftyp   inod |"
    for icol = 1:ncol
        @printf "  %10i" icol 
    end
    @printf "\n__________________________________________|"
    for icol = 1:ncol
        @printf "____________" 
    end
    @printf "\n"

    for iV ∈ 1:nV
        idof = iVdof[iV]
        @printf " %4d    %4d     %16s  %5d |" iV idof field[idof] inod[idof] 
        for icol = 1:ncol
            @printf "  %10.3g" V[iV,icol] 
        end
        @printf "\n"
    end
end


"""
    diffed_lagrangian(eleobj;Λ,X,U,A,t=0.,SP=nothing)

Compute the Lagrangian, its gradients and Hessian, and the memory of an element.
For element debugging and testing. 

The output is a `NamedTuple` with fields `Λ`, `X`, `U`, `A`, `t`, `SP` echoing the inputs and fields
- `∇L` of format `∇L[iclass][ider]`so that for example `∇L[2][3]` contains the gradient of the Lagrangian wrt to the acceleration.
   `iclass` is 1,2,3 and 4 for `Λ`, `X`, `U` and `A` respectively.
- `HL` of format `HL[iclass,jclass][ider,jder]`so that for example `HL[1,2][1,3]` contains the mass matrix.
- `FB` as returned by `lagrangian`

See also: [`diffed_residual`](@ref), [`print_element_array`](@ref)
"""     
function diffed_lagrangian(ele::Eletyp; Λ,X,U,A, t::𝕣=0.,SP=nothing) where{Eletyp<:AbstractElement}
    OX,OU,IA         = length(X)-1,length(U)-1,1
    Nλ               = length(   Λ ) 
    Nx               = length(∂0(X)) 
    Nu               = length(∂0(U)) 
    Na               = length(   A ) 

    if (Nλ,Nx,Nu,Na) ≠ getndof(Eletyp,(:X,:X,:U,:A))
        display(Eletyp)
        @printf("diffed_lagrangian received %i Λ, %i X, %i U and %i A dofs\n",Nλ,Nx,Nu,Na)
        @printf("element requires           %i Λ, %i X, %i U and %i A dofs\n",getndof(Eletyp,:X),getndof(Eletyp,:X),getndof(Eletyp,:U),getndof(Eletyp,:A))
        @assert false
    end
    λxua      = ( 1,    2,    3,  4)
    ndof      = (Nx,   Nx,   Nu, Na)
    nder      = ( 1, OX+1, OU+1, IA)
    Np        = Nx + Nx*(OX+1) + Nu*(OU+1) + Na*IA # number of partials 
    T         = ∂ℝ{2, Np, ∂ℝ{1, Np, Float64}}
    Λ∂        =              SVector{Nx,T}(∂²ℝ{1,Np}(Λ[      idof],                           idof)   for idof=1:Nx)
    X∂        = ntuple(ider->SVector{Nx,T}(∂²ℝ{1,Np}(X[ider][idof],Nx+Nx*(ider-1)            +idof)   for idof=1:Nx),Val(OX+1))
    U∂        = ntuple(ider->SVector{Nu,T}(∂²ℝ{1,Np}(U[ider][idof],Nx+Nx*(OX+1)  +Nu*(ider-1)+idof)   for idof=1:Nu),Val(OU+1))
    A∂        =              SVector{Na,T}(∂²ℝ{1,Np}(A[      idof],Nx+Nx*(OX+1)  +Nu*(OU+1)  +idof)   for idof=1:Na)

    L,FB      = lagrangian(ele, Λ∂,X∂,U∂,A∂,t,SP,(;calledby=:test_element))
    

    ∇Lz,HLz   = value_∂{1,Np}(∂{2,Np}(L))

    ∇L        = Vector{Vector{Any}}(undef,4  )
    HL        = Matrix{Matrix{Any}}(undef,4,4)
    pα        = 0   # points into the partials, 1 entry before the start of relevant partial derivative in α,ider-loop
    for α∈λxua 
        ∇L[α] = Vector{Any}(undef,nder[α])
        for i=1:nder[α]   # we must loop over all time derivatives to correctly point into the adiff-partials...
            iα       = pα.+(1:ndof[α])
            pα      += ndof[α]
            ∇L[α][i] = ∇Lz[iα]
            pβ       = 0
            for β∈λxua 
                HL[α,β] = Matrix{Any}(undef,nder[α],nder[β])
                for j=1:nder[β]
                    iβ   = pβ.+(1:ndof[β])
                    pβ  += ndof[β]
                    HL[α,β][i,j] = HLz[iα,iβ]
                end
            end
        end
    end
    return (Λ=Λ,X=X,U=U,A=A,t=t,SP=SP,∇L=∇L,HL=HL,FB=FB)#,inftyp=inftyp,rettyp=rettyp)
end


"""
    diffed_residual(eleobj;X,U,A,t=0.,SP=nothing)

Compute the residual, its gradients, and the memory of an element.
For element debugging and testing. 

The output is a `NamedTuple` with fields `X`, `U`, `A`, `t`, `SP` echoing the inputs and fields
- `R` containing the residual
- `∇R` of format `∇R[iclass][ider]`so that for example `∇R[2][3]` contains the mass matrix.
   `iclass` is 2,3 and 4 for `X`, `U` and `A` respectively.
- `FB` as returned by `residual`

See also: [`diffed_lagrangian`](@ref), [`print_element_array`](@ref)
"""     
function diffed_residual(ele::Eletyp; X,U,A, t::𝕣=0.,SP=nothing) where{Eletyp<:AbstractElement}
    OX,OU,IA         = length(X)-1,length(U)-1,1
    Nx               = length(∂0(X)) 
    Nu               = length(∂0(U)) 
    Na               = length(   A ) 
    if (Nx,Nu,Na) ≠ getndof(Eletyp,(:X,:U,:A))
        display(Eletyp)
        @printf("diffed_residual received %i X, %i U and %i A dofs\n",Nx,Nu,Na)
        @printf("element requires         %i X, %i U and %i A dofs\n",getndof(Eletyp,:X),getndof(Eletyp,:U),getndof(Eletyp,:A))
        @assert false
    end

    xua       = (    2,    3,  4)
    ndof      = (0, Nx,   Nu, Na)
    nder      = (0 ,OX+1, OU+1, IA)
    Np        = Nx*(OX+1) + Nu*(OU+1) + Na*IA # number of partials 
    X∂        = ntuple(ider->SVector{Nx,∂ℝ{1,Np,𝕣}}(∂ℝ{1,Np}(X[ider][idof],Nx*(ider-1)            +idof)   for idof=1:Nx),Val(OX+1))
    U∂        = ntuple(ider->SVector{Nu,∂ℝ{1,Np,𝕣}}(∂ℝ{1,Np}(U[ider][idof],Nx*(OX+1)  +Nu*(ider-1)+idof)   for idof=1:Nu),Val(OU+1))
    A∂        =              SVector{Na,∂ℝ{1,Np,𝕣}}(∂ℝ{1,Np}(A[      idof],Nx*(OX+1)  +Nu*(OU+1)  +idof)   for idof=1:Na)

    r_,FB     = residual(ele, X∂,U∂,A∂,t,SP,(;calledby=:test_element))
    #inftyp,rettyp = @typeof(residual(ele, X∂,U∂,A∂,t,SP,(;calledby=:test_element)))
    R,∇r      = value_∂{1,Np}(r_)

    ∇R        = Vector{Vector{Any}}(undef,4  )
    pα        = 0   # points into the partials, 1 entry before the start of relevant partial derivative in α,ider-loop
    for α∈xua 
        ∇R[α] = Vector{Any}(undef,nder[α])
        for i=1:nder[α]   # we must loop over all time derivatives to correctly point into the adiff-partials...
            iα       = pα.+(1:ndof[α])
            pα      += ndof[α]
            ∇R[α][i] = ∇r[:,iα]
        end
    end
    return (X=X,U=U,A=A,t=t,SP=SP,R=R,∇R=∇R,FB=FB)#,inftyp=inftyp,rettyp=rettyp)
end

