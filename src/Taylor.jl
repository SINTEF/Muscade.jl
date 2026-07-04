struct motion_{P,Q}       end 
struct motion⁻¹{P,ND,OD}  end 
"""
    P,ND,X_ = motion(X [;context])

Input variables `X` and `U` provided to`residual` or `Lagrangian` are `NTuples` which can
have a form like `(X₀,X₁,X₂)`, with the zeroth, first and second time derivative. `motion` 
transform such a `NTuple` into a `SVector` of `∂ℝ`, containing derivatives with respect to time.

Using a scalar instead of `SVector` to provide an idea, (`motion` actually only operates on `NTuple` of `SVector`s)
this transforms `(x₀,x₁,x₂)` into `x₀+∂₁⟨x₁⟩+∂₂⟨x₁+∂₁⟨x₂⟩⟩`.

This enables an element to use automatic differentiation to use time derivatives of various results, such
as Coriolis and centrifugal accelerations, or strain rates. 

If the third output from `motion` (or variables computed from this 3rd ouput), are combined in computation
with variables obtained by `variate` or `revariate` (including the inputs `Λ`, `X`, `U`, `A` and `t`) to
`lagrangian` and `residual` (or variables computed from these), then the variables obtained by `variate` or `revariate`
must be included in `context`. `context` can be a variable or a ``tuple` of variables.  To prevent bugs, mark variables "touched"
directly or indirectly by `motion`'s 3rd output with e.g. an underscore to keep track of the separation.

    P,ND,X_ = motion(X ,context=(A,t)])
    b       = f(A)         # b is "touched by A", which may be an addifed input to `residual`
    Y_      = g(X_,b,t)    # Y_ is touched by t, and indirectly by A.  This is OK, they where 
                           # declared in the context of X_ 

See [`toolbox.EulerBeam3D`] for an example of usage.

!!! Warning
    Use `∂0`, `∂1` and `∂2` to extract time derivatives from `X` and `U`, do not use
    indexing `X[1]` as the `NTuple` may be truncated in e.g. static analyses.

See [`motion⁻¹`](@ref), [`variate`](@ref)
"""
function motion(a::NTuple{ND,SV{N,R}};context=0) where{ND,N,R} 
    P,_,a_ = variate{0}(a;context) 
    A  = SV{N}(motion_{P,P+ND-1}(let a=a 
                  ntuple(j->a_[j][i],Val(ND))
              end) for i=1:N)
    return P,ND,A          
end
motion_{P,Q  }(a::NTuple{D,      R }) where{D ,P  ,R<:Real,Q} = ∂ℝ{Q,1}(motion_{P,Q-1}(a),SV(motion_{P,Q-1}(a[2:D]))) 
motion_{P,P  }(a::NTuple{D,      R }) where{D ,P  ,R<:Real  } = a[1]

"""
    P,ND,X_  = motion(X)
    Y_       = f(Y_)    
    Y        = motion⁻¹{P,ND}(Y_) 

`motion⁻¹` transforms `Y_` into a `NTuple` `Y`.  The values `P` and ``ND`
must be as returned by `motion`. Time derivatives can be accessed as

    Y₀ = ∂0(Y)
    Y₁ = ∂1(Y)
    Y₂ = ∂2(Y)

One can also obtain the i-th time derivative directly as

    Yᵢ = motion⁻¹{P,ND,i}(Y_)

See also [`motion`](@ref)
"""
motion⁻¹{P,1,0}(a::ℝ) where{P} =                      a
motion⁻¹{P,2,0}(a::ℝ) where{P} =            value{P+1 }(a)
motion⁻¹{P,3,0}(a::ℝ) where{P} = value{P+1}(value{P+2 }(a))
# velocities
motion⁻¹{P,1,1}(a::ℝ) where{P} = 0. 
motion⁻¹{P,2,1}(a::ℝ) where{P} =            ∂{    P+1,1}(a)[1]  # [1]: only partial is wrt time
motion⁻¹{P,3,1}(a::ℝ) where{P} = value{P+1}(∂{    P+2,1}(a)[1])
# accelerations
motion⁻¹{P,1,2}(a::ℝ) where{P} = 0. 
motion⁻¹{P,2,2}(a::ℝ) where{P} = 0.
motion⁻¹{P,3,2}(a::ℝ) where{P} = ∂{   P+1,1}(∂{   P+2,1}(a)[1])[1]
motion⁻¹{P,ND,OD}(a::SArray{S,R}) where{S,P,ND,OD,R<:ℝ}   = SArray{S}(motion⁻¹{P,ND,OD}(aᵢ) for aᵢ∈a)
 
motion⁻¹{P,1    }(a::Union{ℝ,SArray}) where{P   } = (motion⁻¹{P,1,0}(a),)
motion⁻¹{P,2    }(a::Union{ℝ,SArray}) where{P   } = (motion⁻¹{P,2,0}(a),motion⁻¹{P,2,1}(a))
motion⁻¹{P,3    }(a::Union{ℝ,SArray}) where{P   } = (motion⁻¹{P,3,0}(a),motion⁻¹{P,3,1}(a),motion⁻¹{P,3,2}(a))
motion⁻¹{P,ND   }(a::Tuple          ) where{P,ND} = (motion⁻¹{P,ND}(first(a)),motion⁻¹{P,ND}(Base.tail(a))...)   
motion⁻¹{P,ND   }(a::Tuple{}        ) where{P,ND} = ()   
motion⁻¹{P,ND   }(a::NamedTuple     ) where{P,ND} = NamedTuple{keys(a)}(motion⁻¹{P,ND}(values(a)))  
motion⁻¹{P,ND   }(a...              ) where{P,ND} = motion⁻¹{P,ND}(a) 

#############

flat_length(a::NamedTuple   )                = flat_length(values(a))
flat_length(a::Tuple        )                = flat_length(first(a))+flat_length(Base.tail(a))
flat_length(a::Tuple{}      )                = 0
flat_length(a::SArray       )                = length(a) 
flat_length(a::ℝ            )                = 1

flat_eltype(a::NamedTuple   )                = flat_eltype(values(a))
flat_eltype(a::Tuple        )                = promote_type(flat_eltype(first(a)),flat_eltype(Base.tail(a)))
flat_eltype(a::Tuple{X}     ) where{X}       = flat_eltype(a[1]) 
flat_eltype(a::SArray{S,T}  ) where{S,T}     = T  
flat_eltype(a::R            ) where{R<:ℝ}    = R

precedence(a::NamedTuple   )                 = precedence(values(a))
precedence(a::Tuple{X}     ) where{X}        = precedence(a[1]) 
function precedence(a::Tuple        ) 
    pa1 = precedence(first(a))
    pat = precedence(Base.tail(a))
    pa1==0   && return pat
    pat==0   && return pa1
    pa1==pat && return pat
    error("precedence of a tuple with components of various precedence")
end


npartial(  a::NamedTuple   )                 = npartial(values(a))
npartial(  a::Tuple{X}     ) where{X}        = npartial(a[1]) 
function npartial(a::Tuple        ) 
    na1 = npartial(first(a))
    nat = npartial(Base.tail(a))
    na1==0   && return nat
    nat==0   && return na1
    na1==nat && return nat
    error("npartial of a tuple with components of various npartial")
end


struct flatten{T} end
flatten(a)                                   = flatten{flat_eltype(a)}(a)
flatten{T}(a::NamedTuple   ) where{T}        = flatten{T}(values(a))
flatten{T}(a::Tuple        ) where{T}        = SVector{flat_length(a),T}(flatten{T}(first(a))..., flatten{T}(Base.tail(a))...)  
flatten{T}(a::Tuple{}      ) where{T}        = ()
flatten{T}(a::SArray       ) where{T}        = SVector{length(a),T}(a)  
flatten{T}(a::ℝ            ) where{T}        = T(a)

@inline longuestfirst(a,b) = length(a)>length(b) ? (a,b) : (b,a)

@inline function npartials(a::Tuple) 
    long,short = longuestfirst(npartials(first(a)),npartials(Base.tail(a)))
    Δ          = length(long)-length(short)
    for i      = 1:length(short)
        @assert long[Δ+i] == short[i] "Perturbation collision"
    end
    return long # Core.Const
end
@inline npartials(a::Tuple{}  )              = ()
@inline npartials(a::NamedTuple )            = npartials(values(a)) 
@inline npartials(a::SArray   )              = length(a)==0 ? () : npartials(first(a)) # edgecase empty array. Is that OK?
@inline npartials(a::∂ℝ{P,N,R}) where{P,N,R} = (N,npartials(a.x)...)
@inline npartials(a::ℝ        )              = ()

struct front{P} end
@inline front{P}(tup) where{P} = Base.front(front{P-1}(tup))
@inline front{0}(tup)          = tup

@inline zerovalue(v::∂ℝ{P,N,R}) where{P,N,R<:ℝ}  = ∂ℝ{P,N,R}(zerovalue(v.x),v.dx)
@inline zerovalue(v::       R ) where{    R<:ℝ}  = zero(R)

struct AllElements{S} # apply a scaling to a whole Tuple or SArray
    s::S
end

# type for static parametrization
""" 
    V = variate{O}(v[;context,scale])

variate `v` to the `O`-th order. 

    V = variate(v[;context,scale])

variates `v` to the first order.

The variable `v` is a nested structure of `NamedTuple`s, `Tuple`s and `SArrays` of 
`Real`s (possibly: `∂ℝ`s). `V` has the same structure as `v`. 

The output of `variate` must not be combined in computation
with variables obtained by `motion`, `revariate` or a separate call to `variate` (including the inputs `Λ`, `X`, `U`, `A` and `t`) to
`lagrangian` and `residual` (or variables computed from these).  

If output from `variate` (or variables computed from this ouput), are combined in computation
with variables obtained by `revariate`, `motion` or a separate call to `variate` (including the inputs `Λ`, `X`, `U`, `A` and `t`) to
`lagrangian` and `residual` (or variables computed from these), then the later variables
must be included in `context`. `context` can be a variable or a ``tuple` of variables.  

    V = variate(∂1{X};context=(Λ,X,U,t)) # A left out
    W = V*t+∂0(X)                        # safe
    Y = V*A[1]                           # not safe!!!

Not for use within elements: `scale` has the same structure as `v` except that `SArray`s and `NTuple`s can be replaced
by `AllElement(scaling)` to apply `scaling to all elements in the `SArray` or `NTuple` in `v`. 

See also: [`motion`](@ref), [`revariate`](@ref), [`chainrule`](@ref), [`variate_indices`](@ref)
"""
struct variate{O}   end
""" 
    V = variate0{O}(v[;context,scale])

The same as [`variate`](@ref), but while all derivatives are as from [`variate`](@ref), all
values are set to 0.
"""
struct variate0{O}  end # former δ: returns outputs with VALUE=0
""" 
    V_ = revariate{P}(V)

The variable `V` is a nested structure of `NamedTuple`s, `Tuple`s and `SArrays` of 
`Real`s (possibly: `∂ℝ`s).

`V` is stripped of its partials, an revariated to order `P`.

    V_ = revariate(V)

revariates to the order `precedence(V)`.  `TV` has the same structure as `V`. One typical usage is with a `Tuple`:

    a_,b_,c_ = revariate((a,b,c))

`revariate`, in conjunction with `chainrule` can be used to improve performance when the length of 
`V` is smaller than the length of its partials.

The output(s) of `revariate` must not be combined in computation
with variables obtained by `variate` or `revariate` (including the inputs `Λ`, `X`, `U`, `A` and `t`) to
`lagrangian` and `residual` (or variables computed from these).  To prevent bugs, mark variables "touched"
directly or indirectly by `revariate`'s output with e.g. an underscore to keep track of the separation.

See [`toolbox.EulerBeam3D`] for an example of usage.

See also: [`chainrule`](@ref), [`variate_indices`](@ref), [`variate`](@ref)
"""
struct revariate{P} end
struct variate_work{A,B,C} end
struct revariate_work{A} end
struct variate_nested{A,B,C,D,E} end

# API
@inline variate(     v;kwargs...)          = variate_work{1,identity ,:variate}(v;kwargs...)
@inline variate{O}(  v;kwargs...) where{O} = variate_work{O,identity ,:variate}(v;kwargs...)
@inline variate0(    v;kwargs...)          = variate_work{1,zerovalue,:variate}(v;kwargs...)
@inline variate0{O}( v;kwargs...) where{O} = variate_work{O,zerovalue,:variate}(v;kwargs...)
@inline revariate(   v;kwargs...)          = revariate_work{precedence(v)}(v;kwargs...)
@inline revariate{P}(v;kwargs...) where{P} = revariate_work{P            }(v;kwargs...)

# workhorses
@inline function variate_work{O,VZ,:variate}(v;context=0,scale=nothing) where{O,VZ}
    Nv = flat_length(v)    # 𝕫
    Nc = npartials((v,context))  # ntuple(𝕫,Pvc)
    P  = O+length(Nc)              # precedence of outputs
    return if scale==nothing  P,Nv,variate_nested{:variate,O,Nv,Nc,VZ}(v,1      )
    else                      P,Nv,variate_nested{:variate,O,Nv,Nc,VZ}(v,1,scale)
    end
end
@inline function revariate_work{P}(v;scale=nothing) where{P}     
    N = flat_length(v)
    return variate_nested{:revariate,P,N,Nothing,Nothing}(v,1) 
end

# nested types,recursion
@inline variate_nested{A,B,C,D,E}(v::NamedTuple,i               ) where{A,B,C,D,E  } = NamedTuple{keys(v)}(variate_nested{A,B,C,D,E}(values(v),i          )) 
@inline variate_nested{A,B,C,D,E}(v::NamedTuple,i,s::NamedTuple ) where{A,B,C,D,E  } = NamedTuple{keys(v)}(variate_nested{A,B,C,D,E}(values(v),i,values(s))) 
@inline variate_nested{A,B,C,D,E}(v::NamedTuple,i,s::AllElements) where{A,B,C,D,E  } = NamedTuple{keys(v)}(variate_nested{A,B,C,D,E}(values(v),i,s        )) 
@inline variate_nested{A,B,C,D,E}(v::Tuple     ,i               ) where{A,B,C,D,E  } = (variate_nested{A,B,C,D,E}(first(v),i         ),variate_nested{A,B,C,D,E}(Base.tail(v),i+flat_length(first(v))             )...)
@inline variate_nested{A,B,C,D,E}(v::Tuple     ,i,s::Tuple      ) where{A,B,C,D,E  } = (variate_nested{A,B,C,D,E}(first(v),i,first(s)),variate_nested{A,B,C,D,E}(Base.tail(v),i+flat_length(first(v)),Base.tail(s))...)
@inline variate_nested{A,B,C,D,E}(v::Tuple     ,i,s::AllElements) where{A,B,C,D,E  } = (variate_nested{A,B,C,D,E}(first(v),i,s.s     ),variate_nested{A,B,C,D,E}(Base.tail(v),i+flat_length(first(v)),s           )...)
@inline variate_nested{A,B,C,D,E}(v::Tuple{}   ,i               ) where{A,B,C,D,E  } = ()
@inline variate_nested{A,B,C,D,E}(v::Tuple{}   ,i,s::Tuple{}    ) where{A,B,C,D,E  } = ()
@inline variate_nested{A,B,C,D,E}(v::Tuple{}   ,i,s::AllElements) where{A,B,C,D,E  } = ()
@inline variate_nested{A,B,C,D,E}(v::SArray{S} ,i               ) where{A,B,C,D,E,S} = variate_nested{A,B,C,D,E}.(v,SArray{S,𝕫}(i-1+j for j∈eachindex(v))                                              )
@inline variate_nested{A,B,C,D,E}(v::SArray{S} ,i,s::SArray{S}  ) where{A,B,C,D,E,S} = variate_nested{A,B,C,D,E}.(v,SArray{S,𝕫}(i-1+j for j∈eachindex(v)),s                                            ) 
@inline variate_nested{A,B,C,D,E}(v::SArray{S} ,i,s::AllElements) where{A,B,C,D,E,S} = variate_nested{A,B,C,D,E}.(v,SArray{S,𝕫}(i-1+j for j∈eachindex(v)),SArray{S,typeof(s.s)}(s.s for j∈eachindex(v)))
@inline variate_nested{A,B,C,D,E}(v::ℝ         ,i               ) where{A,B,C,D,E  } = variate_nested{A,B,C,D,E}(v,i,1.)


# add O variations on the outside
@inline function variate_nested{:variate,O  ,Nv,Nc,VZ}(v::ℝ,i,s) where{O,Nv,Nc,VZ} 
    x  =         variate_nested{:variate,O-1,Nv,Nc,VZ}(v   ,i,s)
    R  = typeof(x)
    dx = SVector{Nv,R}(i==j ? R(s) : zero(R) for j=1:Nv)
    return ∂ℝ{O+length(Nc),Nv,R}(x,dx)
end 
# transition to getting v up to Nc
@inline function variate_nested{:variate, 0 ,Nv       ,Nc,VZ}(v::ℝ,i,s) where{Nv,Nc,VZ} 
    Nc_ = front{precedence(v)}(Nc) # cut the tail of Nc if v is not 𝕣
    return       variate_nested{:variate,-1  ,Nothing,Nc_,VZ}(v,i,s) 
end
# get v up to Nc
@inline          variate_nested{:variate,-1,Nothing,()           ,VZ}(v::ℝ,i,s) where{VZ} = VZ(v)
@inline function variate_nested{:variate,-1,Nothing,          Nc ,VZ}(v::ℝ,i,s) where{Nc,VZ}
    x     =      variate_nested{:variate,-1,Nothing,Base.tail(Nc),VZ}(v   ,i,s)
    P,N,R = precedence(v)+length(Nc),first(Nc),typeof(x)
    dx    = SVector{N,R}(zero(R) for j=1:N)
    return ∂ℝ{P,N,R}(x,dx)
end

# revariate
@inline variate_nested{         :revariate,0  ,N,D,E}(v::ℝ,i,s) where{  N,D,E} = VALUE(v)
@inline function variate_nested{:revariate,P  ,N,D,E}(v::ℝ,i,s) where{P,N,D,E} 
    x  = variate_nested{:revariate,P-1,N,D,E}(v   ,i,s)
    R  = typeof(x)
    dx = SVector{N,R}(i==j ? R(s) : zero(R) for j=1:N)
    return ∂ℝ{P,N,R}(x,dx)
end 


"""
    iV   = variate_indices(V)

For use in conjunction with
    
    TV   = revariate(V)

`iV` has the same structure as `V` and `TV` but contains integers: the indices into the partials of `TV` 

See also: [`revariate`](@ref)

"""
variate_indices( a              ) = revariate_indices_(a,0)
revariate_indices_(a::Tuple     ,i) = (revariate_indices_(first(a),i),revariate_indices_(Base.tail(a),i+flat_length(first(a)))...)
revariate_indices_(a::Tuple{}   ,i) = ()
revariate_indices_(a::NamedTuple,i) = NamedTuple{keys(a)}(revariate_indices_(values(a),i))
#revariate_indices_(a::SArray    ,i) = i+1:i+flat_length(a)
revariate_indices_(a::SArray    ,i) = SVector{flat_length(a),𝕫}(i+iₐ for iₐ∈1:flat_length(a))
revariate_indices_(a::𝕣         ,i) = i+1

"""
    to_order{P,N}(V)

Decrease (lossy) or increase (pad partials with zeros) the order of differentiation of `V`.
`V` is a nested structure of `NamedTuple`, `Tuple`, `SArray`, and the components
of `V` must be of type `∂ℝ` or `𝕣`.

IMPORTANT:
1) assumes all the nested differentiations are with respect to the same variable. 
2) to_order{P,N}(V::𝕣) = V
"""
struct to_order{P,N} end
to_order{P,N}(a::NamedTuple   ) where{P ,N} = NamedTuple{keys(a)}(to_order{P,N}(values(a)))
to_order{P,N}(a::Tuple        ) where{P ,N} = (to_order{P,N}(first(a)),to_order{P,N}(Base.tail(a))...)
to_order{P,N}(a::Tuple{}      ) where{P ,N} = ()
to_order{P,N}(a::AbstractArray) where{P ,N} = to_order{P,N}.(a)
to_order{0,N}(a::𝕣            ) where{   N} = a
to_order{P,N}(a::𝕣            ) where{P, N} = a  # deliberate!
to_order{0,N}(a::∂ℝ{Pa,N}     ) where{Pa,N} = VALUE(a) 
function to_order{P,N}(a::∂ℝ{Pa,N}) where{P,N,Pa}
    if     Pa==P
        a
    elseif Pa> P   
        to_order{P,N}(value{Pa}(a))
    elseif Pa< P
        ap  = to_order{P-1,N}(a) 
        ∂ℝ{P,N}(ap, SVector{N}(∂ℝ{P-1,N}(ap.dx[i]) for i=1:N) )
    end
end



#########################

struct McLaurin{d} end
"""
    v=McLaurin(Ty,Δx)

`Ty::∂ℝ` has partials to arbitrary order with respect to a variable `x`. These
partials define a McLaurin expansion, which `McLaurin` evaluates at value `Δx`.

Thus the elements in `v` have the same type as `Δx`

`McLaurin` handles nested structures of `Tuple`s and `SVector`s of `∂ℝ`, applying the
expansion to each element.

`McLaurin` is a utility function behind [`chainrule`](@ref) and [`Taylor`](@ref)

See also: [`chainrule`](@ref), [`Taylor`](@ref), [`revariate`](@ref), [`apply`](@ref)    
"""
McLaurin(y::Tuple,Δx)                          = tuple(McLaurin(first(y),Δx),McLaurin(Base.tail(y),Δx)...) 
McLaurin( ::Tuple{},Δx)                        = tuple() 
McLaurin(y::SArray{Sy,Ty,Dy,Ly},Δx::SVector{Sx,Tx}) where{Sy,Ty<:∂ℝ,Dy,Ly,Sx,Tx} = SArray{Sy,Tx,Dy,Ly}(McLaurin(yᵢ,Δx) for yᵢ∈y)
McLaurin(y::SArray{Sy,Ty,Dy,Ly},Δx::SVector{Sx,Tx}) where{Sy,Ty<: ℝ,Dy,Ly,Sx,Tx} = SArray{Sy,Ty,Dy,Ly}(McLaurin(yᵢ,Δx) for yᵢ∈y)
McLaurin(y::ℝ,Δx)                              = McLaurin{1}(y::ℝ,Δx)
McLaurin{D}(y::𝕣 ,Δx::SVector{N}) where{N,D}   = y
function McLaurin{D}(y::∂ℝ{P,N,R},Δx::SVector{N}) where{P,N,R,D} # Horner's algorithm
    @assert N>0 # if this edge case is relevant, create specialised method that converts VALUE(y) to eltype(Δx)
    o = McLaurin{D+1}(y.dx[1],Δx)*Δx[1]
    for i ∈ 2:N
        o += McLaurin{D+1}(y.dx[i],Δx)*Δx[i]
    end
    D>1 && (o *= inv(D))
    return VALUE(y)+o
end


"""
    Taylor(Ty,x₀,x)

`Ty::∂ℝ` has partials to arbitrary order evaluated at `x₀`. These
partials define a Taylor expansion, which `Taylor` evaluates at value `x`

`Taylor` handles nested structures of `Tuple`s and `SVector`s of `∂ℝ`, applying the
expansion to each element.

See also: [`chainrule`](@ref), [`McLaurin`](@ref), [`revariate`](@ref), [`apply`](@ref)    
"""
Taylor(y::Tuple,x₀,x) = McLaurin(y,x-x₀)

"""
    chainrule(Ty,x)

Apply a chain rule in automatic differentiation.  For example 
    Tx = revariate(x)
    Ty = f(Tx)
    y  = chainrule(Ty,x)    
is faster than
    y  = f(x)
if the length of `x` is smaller than the length of its partials.

See also: [`revariate`](@ref), [`apply`](@ref)    
"""
function chainrule(Ty,x) 
    fx = flatten(x)
    McLaurin(Ty,fx.-VALUE(fx))
end

"""
    y,... = apply{:chainrule}(f,x)
    y,... = apply{:direct   }(f,x)

In the context of forward automatic differentiation using `∂ℝ`, `apply{:chainrule}`
accelerates the evaluation of `y,...= f(x)` if the length of `x` is smaller than 
the length of its partials.

`apply{:direct}` simply executes `f(x)` (no chain rule is applied)

Also works where `x` is a nested structure of `Tuple`s and `NamedTuple`s where the leaves
are `ℝ` or `SArray{S,R} where {S,R<:ℝ}`.    

!!! warning
    If `f` is a closure, make sure that `f` does not capture variables of type `∂ℝ`.

!!! warning
    See [`chainrule`](@ref) about when it is advisable to use chainrule, and when not.

Wrapper function of [`revariate`](@ref) and [`chainrule`](@ref) 

"""
struct apply{Mode} end
function apply{:chainrule}(f,x) 
    x_ = revariate(x)
    y_ = f(x_)
    y  = chainrule(y_,x)
    return y
end
apply{:direct   }(f,x) = f(x)

"""
    chainrule_Jacobian{P}(Ty,X_)

Given `Ty` obtained using `revariate`, and `X_`, obtained using `motion{P}(X)` where `X` is a tuple
of `SVectors` and `P=precedence(X)`, compute `y`, a tuple of length `ND` of `AbstractArrays` of same `eltype` as vectors in `X,
and `y∂X₀`, the Jacobian of `∂0(y)` with respect to `∂0(X)`.

See also [`revariate`](@ref), [`motion`](@ref), [`motion⁻¹`](@ref) 
"""
struct chainrule_Jacobian{P} end
chainrule_Jacobian(Ty            ,X₀)  = chainrule(∂{precedence(Ty),npartial(Ty)}(Ty),X₀) # y∂X₀
chainrule_Jacobian(Ty::NamedTuple,X₀)  = NamedTuple{keys(Ty)}(chainrule_Jacobian(values(Ty),X₀))
chainrule_Jacobian(Ty::Tuple     ,X₀)  = (chainrule_Jacobian(first(Ty),X₀),chainrule_Jacobian(Base.tail(Ty),X₀)...)
chainrule_Jacobian(Ty::Tuple{}   ,X₀)  = ()


