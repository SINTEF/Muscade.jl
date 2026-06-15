# TODO
# implement sincos
# differentiate x^0 at x=0

struct T𝟘 <: Real end  # Bbbzero
const  𝟘 = T𝟘()
@inline Base.:(+)(a::𝕣 , ::T𝟘) = a
@inline Base.:(+)( ::T𝟘,b::𝕣 ) = b
@inline Base.:(+)( ::T𝟘, ::T𝟘) = 𝟘
@inline Base.:(-)(a::𝕣 , ::T𝟘) = a
@inline Base.:(-)( ::T𝟘,b::𝕣 ) = -b
@inline Base.:(-)( ::T𝟘, ::T𝟘) = 𝟘
@inline Base.:(*)(a::𝕣 , ::T𝟘) = 𝟘
@inline Base.:(*)( ::T𝟘,b::𝕣 ) = 𝟘
@inline Base.:(*)( ::T𝟘, ::T𝟘) = 𝟘
@inline Base.:(/)(a    , ::T𝟘) = error("division by 𝟘")
@inline Base.:(/)( ::T𝟘,b    ) = 𝟘
Base.show(io::IO, ::T𝟘) = print(io,"𝟘")

## Type and construction
const SV = SVector  
const SA = SArray 
const SM = SMatrix
# Types
# P precedence.  Newer, derivatives, outest in the adiff datastructure have higher numbers  
# N number of partials 
# R type of the variable  (and partials)
struct ∂ℝ{P,N,R} <:ℝ where{R<:ℝ}  # P for precedence, N number of partials, R type of the variable (∂ℝ can be nested)
    x  :: R
    dx :: SV{N,R}
end

# Constructors 
∂ℝ{P,N  }(x::R ,dx::SV{N,R}) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(  x ,SV{N,R}(dx                               ))
∂ℝ{P,N  }(x::R             ) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(  x ,SV{N,R}(zero(R)                 for j=1:N))
∂ℝ{P,N  }(x::R,i::ℤ        ) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(  x ,SV{N,R}(i==j ? one(R) : zero(R) for j=1:N))
∂ℝ{P,N  }(x::R,i::ℤ,s      ) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(  x ,SV{N,R}(i==j ? R(s)   : zero(R) for j=1:N))
∂ℝ{P,N,R}(x::𝕣             ) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(R(x),SV{N,R}(zero(R)                 for j=1:N))
function ∂ℝ{P,N}(x::Rx,dx::SV{N,Rdx}) where{P,N,Rx<:ℝ,Rdx<:ℝ}
    R = promote_type(Rx,Rdx)
    return ∂ℝ{P,N}(convert(R,x),convert.(R,dx))
end

# zeros, ones
Base.zero(::Type{∂ℝ{P,N,R}}) where{P,N,R<:ℝ}     = ∂ℝ{P,N,R}(zero(R), SV{N,R}(zero(R) for j=1:N))
Base.one( ::Type{∂ℝ{P,N,R}}) where{P,N,R<:ℝ}     = ∂ℝ{P,N,R}(one( R), SV{N,R}(zero(R) for j=1:N))
Base.isnan(   a::∂ℝ)                             = isnan(   VALUE(a))
Base.isone(   a::∂ℝ)                             = isone(   VALUE(a))
Base.iszero(  a::∂ℝ)                             = iszero(  VALUE(a))
Base.isinf(   a::∂ℝ)                             = isinf(   VALUE(a))
Base.isfinite(a::∂ℝ)                             = isfinite(VALUE(a))
Base.typemax( ::Type{∂ℝ{P,N,R}}) where{P,N,R<:ℝ} = typemax(R)
Base.typemin( ::Type{∂ℝ{P,N,R}}) where{P,N,R<:ℝ} = typemin(R)
Base.floatmax(::Type{∂ℝ{P,N,R}}) where{P,N,R<:ℝ} = floatmax(R)
Base.floatmin(::Type{∂ℝ{P,N,R}}) where{P,N,R<:ℝ} = floatmin(R)
Base.floatmax(::     ∂ℝ{P,N,R} ) where{P,N,R<:ℝ} = floatmax(R)  # because ℝ is Real, not AbstractFloat
Base.floatmin(::     ∂ℝ{P,N,R} ) where{P,N,R<:ℝ} = floatmin(R)  # because ℝ is Real, not AbstractFloat
Base.eps(     ::Type{∂ℝ{P,N,R}}) where{P,N,R<:ℝ} = eps(R)
Base.float(a::∂ℝ)                                = a

# promote rules
Base.promote_rule(::Type{∂ℝ{P ,N ,Ra}},::Type{∂ℝ{P,N,Rb}}) where{P ,N ,Ra<:ℝ,Rb<:ℝ} = ∂ℝ{P ,N ,promote_type(Ra,Rb)}
Base.promote_rule(::Type{∂ℝ{Pa,Na,Ra}},::Type{       Rb }) where{Pa,Na,Ra<:ℝ,Rb<:ℝ} = ∂ℝ{Pa,Na,promote_type(Ra,Rb)}
function Base.promote_rule(::Type{∂ℝ{Pa,Na,Ra}},::Type{∂ℝ{Pb,Nb,Rb}}) where{Pa,Pb,Na,Nb,Ra<:ℝ,Rb<:ℝ}
    if  Pa>Pb ∂ℝ{Pa,Nb,promote_type(      Ra    ,∂ℝ{Pb,Nb,Rb})}
    else      ∂ℝ{Pb,Nb,promote_type(∂ℝ{Pa,Na,Ra},      Rb    )}
    end
end

# conversions
Base.convert(::Type{∂ℝ{P,N,Ra}},b::∂ℝ{P,N,Rb}) where{P,N,Ra<:ℝ,Rb<:ℝ} = ∂ℝ{P ,N }(convert(Ra,b.x) ,convert.(Ra,b.dx))
Base.convert(::Type{∂ℝ{P,N,Ra}},b::ℝ         ) where{P,N,Ra<:ℝ      } = ∂ℝ{P ,N }(convert(Ra,b  ) ,SV{N,Ra}(zero(Ra) for j=1:N))
function Base.convert(::Type{∂ℝ{Pa,Na,Ra}},b::∂ℝ{Pb,Nb,Rb}) where{Pa,Pb,Na,Nb,Ra<:ℝ,Rb<:ℝ}
    if Pa> Pb return                                                    ∂ℝ{Pa,Na}(convert(Ra,b.x) ,convert.(Ra,b.dx))
    else      muscadeerror(print("Cannot convert precedence ",Pb," to ",Pa))
    end
end

# Pack and unpack
"""
    Muscade.precedence(a::∂ℝ{P,N,R}) → P
    Muscade.precedence(typeof(a)) → P

Also handles static arrays and tuples.    

See also: [`npartial`](@ref)

"""  
precedence  
"""
    Muscade.npartial(a::∂ℝ{P,N,R}) → N
    Muscade.npartial(typeof(a)) → N

Also handle static arrays and tuples.    

See also: [`precedence`](@ref)

""" 
npartial

precedence( ::Type{<:∂ℝ{P,N,R}}) where{P,N,R<:ℝ}          = P
npartial(   ::Type{<:∂ℝ{P,N,R}}) where{P,N,R<:ℝ}          = N
partialtype(::Type{<:∂ℝ{P,N,R}}) where{P,N,R<:ℝ}          = R
precedence( ::Type{Union{}})                              = 0
precedence( ::Type{<:ℝ})                                  = 0
npartial(   ::Type{<:ℝ})                                  = 0
partialtype(::Type{<:ℝ})                                  = Union{} 

precedence(a::ℝ)      = precedence(typeof(a))
npartial(  a::ℝ)      = npartial(typeof(a))
precedence( a,args...) = max(precedence(a),precedence(args...))
precedence(a::SA)     = precedence(eltype(a))
npartial(  a::SA)     = npartial(eltype(a))

function ∂ℝ_promote_type(::Ra,::Rb) where{Ra<:ℝ,Rb<:ℝ}
    Pa,Na = Muscade.precedence(Ra),Muscade.npartial(  Ra)
    Pb,Nb = Muscade.precedence(Rb),Muscade.npartial(  Rb)
    P     = Pa≥Pb ? Pa : Pb
    N     = Pa≥Pb ? Na : Nb
    R     = promote_type(partialtype(Ra),partialtype(Rb))
    return P,N,R
end

"""
    P = constants(a,b,c)

Generate a precedence `P` that is higher than the precedence of the arguments.   

See also: [`variate`](@ref), [`δ`](@ref), [`value`](@ref), [`∂`](@ref), [`VALUE`](@ref), [`value_∂`](@ref)
"""
constants(tup::Tuple)   = constants(tup...) 
constants( a,args...)   = max(constants(a),constants(args...))
constants( a)           = 1+precedence(a) 
constants(tup::Tuple{}) = 1 
constants( ::Nothing)   = 1

# variate
struct δ{P,N,R}                end # need dum, because syntax δ{P,N,R}() collides with default constructor
struct variate{P,N}            end
struct directional{P,N}        end 
struct ∂²ℝ{P,N}                end
"""
    X = δ{P,N,R}()

create a `SVector` of automatic differentiation objects of precedence `P` and value `zero`.    

    X = δ{P}()

Create automatic differentiation object of precedence `P` and value `zero`.  

See also: [`constants`](@ref), [`variate`](@ref), [`value`](@ref), [`∂`](@ref), [`VALUE`](@ref), [`value_∂`](@ref)
"""
@inline δ{P,N,R}(                          ) where{P,N,R<:ℝ} = SV{N,∂ℝ{P,N,R}}(∂ℝ{P,N  }(zero(R),i                                         ) for i=1:N)
@inline δ{P,N,R}(               δa::SV{N,𝕣}) where{P,N,R<:ℝ} = SV{N,∂ℝ{P,N,R}}(∂ℝ{P,N,R}(zero(R),SV{N,R}(i==j ? δa[i]  : zero(R) for i=1:N)) for j=1:N)
@inline δ{P    }(                          ) where{P       } =                 ∂ℝ{P,1,𝕣}(0.     ,SV{1,𝕣}(1.                               ))

"""
    X = variate{P,N}(x)

where `typeof(x)<:SVector{N}`, create a `SVector` of automatic differentiation objects of precedence `P`.    

    X = variate{P}(x)

where `typeof(x)<:Real`, create an object of precedence `P`.    

See also: [`constants`](@ref), [`δ`](@ref), [`value`](@ref), [`∂`](@ref), [`VALUE`](@ref), [`value_∂`](@ref)
"""
@inline variate{    P,N}(a::SV{N,R}            ) where{P,N,R<:ℝ} = SV{N,∂ℝ{P,N,R}}(∂ℝ{P,N  }(a[i],i) for i=1:N)
@inline variate{    P,N}(a::SV{N,R},δa::SV{N,𝕣}) where{P,N,R<:ℝ} = SV{N,∂ℝ{P,N,R}}(∂ℝ{P,N,R}(a[j]   ,SV{N,R}(i==j ? R(δa[i])  : zero(R) for i=1:N)) for j=1:N)
@inline variate{    P  }(a::R                  ) where{P,  R<:ℝ} =      ∂ℝ{P,1  }(a,SV{1,R}(one(R)))
@inline directional{P  }(a::SV{N,R},δa::SV{N,R}) where{P,N,R<:ℝ} = SV{N,∂ℝ{P,1,R}}(∂ℝ{P,1}(a[i],SV{1,R}(δa[i])) for i=1:N)
function ∂²ℝ{P,N}(x::R,i,s=one(R)) where{P,N,R<:ℝ}
    R1 = ∂ℝ{P,N,R}
    return ∂ℝ{P+1,N,R1}( ∂ℝ{P,N}(x,i,s) , SV{N,R1}(j==i ? R1(s) : zero(R1) for j=1:N ))
end
# Analyse
"""
    VALUE(Y)

Completely strip `Y` of partial derivatives.  Use only for debugging purpose.    

See also: [`constants`](@ref), [`variate`](@ref), [`δ`](@ref), [`value`](@ref), [`∂`](@ref), [`value_∂`](@ref)
"""
@inline VALUE(a::Nothing )                     =        nothing
@inline VALUE(a::ℝ )                           =        a
@inline VALUE(a::∂ℝ)                           = VALUE( a.x)
@inline VALUE(a::SA)                           = VALUE.(a)

struct ∂{P,N}                  end 
struct value{P,N}              end
struct value_∂{P,N}            end

"""
    y = value{P}(Y)

Extract the value of an automatic differentiation object, or `SArray` of such objects.    

See also: [`constants`](@ref), [`variate`](@ref), [`δ`](@ref), [`∂`](@ref), [`VALUE`](@ref), [`value_∂`](@ref)
"""
@inline value{P}(a::∂ℝ{P,N,R}) where{P,N,R   } = a.x
@inline value{P}(a::R        ) where{P  ,R<:ℝ} = a
@inline value{P}(a::SA       ) where{P       } = value{P}.(a)
@inline value{P}(a::Tuple    ) where{P       } = Tuple(value{P}(aᵢ) for aᵢ∈a) 

# ∂{P}(a) is handled as ∂{P,1}(a) and returns a scalar 
"""
    yₓ = ∂{P,N}(Y)

Extract the gradient of an automatic differentiation object.  If `Y` is a `SArray`, 
the index of the partial derivative is appended to the indices of `Y`.   

    y′ = ∂{P}(Y)

Extract the derivative of an automatic differentiation object (or `SArray` of such), where the variation
was created by the syntax `variate{P}`.

See also: [`constants`](@ref), [`variate`](@ref), [`δ`](@ref), [`value`](@ref), [`VALUE`](@ref), [`value_∂`](@ref)
"""
@inline ∂{P,N}(a::     ∂ℝ{P,N,R} ) where{  P,N,R   } = a.dx
@inline ∂{P,N}(a::            R  ) where{  P,N,R<:ℝ} = SV{  N,R}(zero(R)    for i=1:N      )
@inline ∂{P,N}(a::SV{M,∂ℝ{P,N,R}}) where{M,P,N,R   } = SM{M,N,R}(a[i].dx[j] for i=1:M,j∈1:N) # ∂(a,x)[i,j] = ∂a[i]/∂x[j]
@inline ∂{P,N}(a::SV{M,       R }) where{M,P,N,R   } = SM{M,N,R}(zero(R)    for i=1:M,j=1:N)
@inline ∂{P  }(a::            R  ) where{  P,  R<:ℝ} = zero(R)
@inline ∂{P  }(a::     ∂ℝ{P,1,R} ) where{  P,  R   } = a.dx[1]
@inline ∂{P  }(a::SV{N,∂ℝ{P,1,R}}) where{  P,N,R   } = SV{  N,R}(a[i].dx[1] for i=1:N     ) # ∂(a,x)[i]    = ∂a[i]/∂x

@inline ∂{P,N}(a::SM{      M1,M2       ,∂ℝ{P,N,R}}) where{M1,M2      ,P,N,R} = SA{Tuple{M1,M2,N       },R}(a[i].dx[j] for i∈eachindex(a),j∈1:N) # ∂(a,x)[i,...,j] = ∂a[i,...]/∂x[j]
@inline ∂{P,N}(a::SM{      M1,M2       ,       R }) where{M1,M2      ,P,N,R} = SA{Tuple{M1,M2,N       },R}(zero(R)    for i∈eachindex(a),j∈1:N)
@inline ∂{P,N}(a::SA{Tuple{M1,M2,M3   },∂ℝ{P,N,R}}) where{M1,M2,M3   ,P,N,R} = SA{Tuple{M1,M2,M3    ,N},R}(a[i].dx[j] for i∈eachindex(a),j∈1:N) # ∂(a,x)[i,...,j] = ∂a[i,...]/∂x[j]
@inline ∂{P,N}(a::SA{Tuple{M1,M2,M3   },       R }) where{M1,M2,M3   ,P,N,R} = SA{Tuple{M1,M2,M3    ,N},R}(zero(R)    for i∈eachindex(a),j∈1:N)
@inline ∂{P,N}(a::SA{Tuple{M1,M2,M3,M4},∂ℝ{P,N,R}}) where{M1,M2,M3,M4,P,N,R} = SA{Tuple{M1,M2,M3,M4 ,N},R}(a[i].dx[j] for i∈eachindex(a),j∈1:N) # ∂(a,x)[i,...,j] = ∂a[i,...]/∂x[j]
@inline ∂{P,N}(a::SA{Tuple{M1,M2,M3,M4},       R }) where{M1,M2,M3,M4,P,N,R} = SA{Tuple{M1,M2,M3,M4 ,N},R}(zero(R)    for i∈eachindex(a),j∈1:N)
@inline ∂{P,N}(a::Tuple                           ) where{            P,N  } = Tuple(∂{P,N}(aᵢ) for aᵢ∈a) 
"""
    y,yₓ = value_∂{P,N}(Y)
    y,y′ = value_∂{P  }(Y)
    
Get value and derivative in one operation.    

See also: [`constants`](@ref), [`variate`](@ref), [`δ`](@ref), [`value`](@ref), [`∂`](@ref), [`VALUE`](@ref)
"""
@inline value_∂{P,N}(a) where{  P,N}= value{P}(a),∂{P,N}(a)
@inline value_∂{P  }(a) where{  P  }= value{P}(a),∂{P  }(a)

# 
struct  ∂𝟘{P} end
@inline ∂𝟘{P}(a::∂ℝ{P}) where{P} = a.dx
@inline ∂𝟘{P}(a:: ℝ   ) where{P} = 𝟘
struct value_∂𝟘{P} end
@inline value_∂𝟘{P}(a::ℝ) where{P}= value{P}(a), ∂𝟘{P}(a)



## Binary operations
for OP∈(:(>),:(<),:(==),:(>=),:(<=),:(!=))
    @eval @inline Base.$OP(a::∂ℝ,b::∂ℝ)  = $OP(VALUE(a),VALUE(b))
    @eval @inline Base.$OP(a:: ℝ,b::∂ℝ)  = $OP(      a ,VALUE(b))
    @eval @inline Base.$OP(a::∂ℝ,b:: ℝ)  = $OP(VALUE(a),      b )
end
@inline Base.isapprox(a::∂ℝ,b::∂ℝ;kwargs...) = isapprox(VALUE(a),VALUE(b);kwargs...)  
@inline Base.isapprox(a:: ℝ,b::∂ℝ;kwargs...) = isapprox(      a ,VALUE(b);kwargs...)  
@inline Base.isapprox(a::∂ℝ,b:: ℝ;kwargs...) = isapprox(VALUE(a),      b ;kwargs...)  

macro DiffRule2(OP,AB,A,B)
    return esc(quote
        @inline $OP(a::∂ℝ{P,N,R},b::∂ℝ{P,N,R}) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}($OP(a.x,b.x),$AB)
        @inline $OP(a::∂ℝ{P,N,R},b::ℝ        ) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}($OP(a.x,b  ),$A )
        @inline $OP(a::ℝ        ,b::∂ℝ{P,N,R}) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}($OP(a  ,b.x),$B )
        @inline function $OP(a::∂ℝ{Pa,Na,Ra},b::∂ℝ{Pb,Nb,Rb}) where{Pa,Pb,Na,Nb,Ra<:ℝ,Rb<:ℝ}
            if Pa==Pb
                R = promote_type(Ra,Rb)
                return ∂ℝ{Pa,Na}(convert(R,$OP(a.x,b.x)),convert.(R,$AB))
            elseif Pa> Pb
                R = promote_type(Ra,typeof(b))
                return ∂ℝ{Pa,Na}(convert(R,$OP(a.x,b  )),convert.(R,$A ))
            else
                R = promote_type(typeof(a),Rb)
                return ∂ℝ{Pb,Nb}(convert(R,$OP(a  ,b.x)),convert.(R,$B ))
            end
        end
    end)
end
@DiffRule2(Base.atan,  (a.dx*b.x-b.dx*a.x)/(a.x^2+b.x^2),          (a.dx*b)/(a.x^2+b^2), -(b.dx*a)/(a^2+b.x^2) )   
@DiffRule2(Base.hypot, (a.dx*a.x+b.dx*b.x)/hypot(a.x,b.x),         a.dx*a.x/hypot(a.x,b), b.dx*b.x/hypot(a,b.x))   
@DiffRule2(Base.:(+),  a.dx+b.dx,                                  a.dx,                  b.dx                 )
@DiffRule2(Base.:(-),  a.dx-b.dx,                                  a.dx,                  -b.dx                )
@DiffRule2(Base.:(*),  a.dx*b.x+a.x*b.dx,                          a.dx*b,                a*b.dx               )
@DiffRule2(Base.:(/),  a.dx/b.x-a.x/b.x^2*b.dx,                    a.dx/b,                -a/b.x^2*b.dx        ) 
@DiffRule2(Base.:(^),  a.dx*b.x*a.x^(b.x-1)+log(a.x)*a.x^b.x*b.dx, a.dx*b*a.x^(b  -1),    log(a)*a ^b.x*b.dx   )  # for exponents ∈ ℝ
@inline Base.:(^)(a::∂ℝ{P,N,R},b::ℤ) where{P,N,R<:ℝ} = b==0 ? one(a) : ∂ℝ{P,N,R}(a.x^b ,a.dx*b*a.x^(b-1) )       # for exponents ∈ ℤ





## Functions
macro DiffRule1(OP,A)
    return esc(:(@inline $OP(a::∂ℝ{P,N}) where{P,N} = ∂ℝ{P,N}($OP(a.x),$A)))
end
@DiffRule1(Base.:(+),       a.dx                                                     )
@DiffRule1(Base.:(-),      -a.dx                                                     )
@DiffRule1(Base.abs  ,a.x==0.0 ? zero(a.dx) : (a.x>0.0 ? a.dx : -a.dx)               )
@DiffRule1(Base.conj ,      a.dx                                                     )
@DiffRule1(Base.sqrt,       a.dx / 2. / sqrt(a.x)                                    )
@DiffRule1(Base.cbrt,       a.dx / 3. / cbrt(a.x)^2                                  )
@DiffRule1(Base.abs2,       a.dx*2. * a.x                                            )
#@DiffRule1(Base.inv,       -a.dx * abs2(inv(a.x))                                    )
@inline function Base.inv(x::∂ℝ{P,N}) where{P,N}  
    x⁻¹ = inv(x.x)  # minimal gain from computing x⁻¹ once, not twice
    ∂ℝ{P,N}(x⁻¹,x.dx.*(-abs2(x⁻¹)))
end
@DiffRule1(Base.log,        a.dx / a.x                                               )
@DiffRule1(Base.log10,      a.dx / a.x / log(10.)                                    )
@DiffRule1(Base.log2,       a.dx / a.x / log(2.)                                     )
@DiffRule1(Base.log1p,      a.dx / (a.x + 1.)                                        )
@DiffRule1(Base.exp,         exp(a.x) * a.dx                                         )
@DiffRule1(Base.exp2,        log(2. ) * exp2( a.x) * a.dx                            )
@DiffRule1(Base.exp10,       log(10.) * exp10(a.x) * a.dx                            )
@DiffRule1(Base.expm1,       exp(a.x) * a.dx                                         )
@DiffRule1(Base.sin,         cos(a.x) * a.dx                                         )
@DiffRule1(Base.cos,        -sin(a.x) * a.dx                                         )
@DiffRule1(Base.tan,         (1. + tan(a.x)^2) * a.dx                                )
@DiffRule1(Base.sinpi,       π*cos(a.x) * a.dx                                       )
@DiffRule1(Base.cospi,      -π*sin(a.x) * a.dx                                       )
@DiffRule1(Base.sec,         sec(a.x) * tan(a.x) * a.dx                              )
@DiffRule1(Base.csc,        -csc(a.x) * cot(a.x) * a.dx                              )
@DiffRule1(Base.cot,        -(1. + cot(a.x)^2) * a.dx                                )
@DiffRule1(Base.sind,        π / 180. * cosd(a.x) * a.dx                             )
@DiffRule1(Base.cosd,       -π / 180. * sind(a.x) * a.dx                             )
@DiffRule1(Base.tand,        π / 180. * (1. + tand(a.x)^2) * a.dx                    )
@DiffRule1(Base.secd,        π / 180. * secd(a.x) * tand(a.x) * a.dx                 )
@DiffRule1(Base.cscd,       -π / 180. * cscd(a.x) * cotd(a.x) * a.dx                 )
@DiffRule1(Base.cotd,       -π / 180. * (1. + cotd(a.x)^2)  * a.dx                   )
@DiffRule1(Base.asin,        a.dx / sqrt(1. - a.x^2)                                 )
@DiffRule1(Base.acos,       -a.dx / sqrt(1. - a.x^2)                                 )
@DiffRule1(Base.atan,        a.dx / (1. + a.x^2)                                     )
@DiffRule1(Base.asec,        a.dx / abs(a.x) / sqrt(a.x^2 - 1.)                      )
@DiffRule1(Base.acsc,       -a.dx / abs(a.x) / sqrt(a.x^2 - 1.)                      )
@DiffRule1(Base.acot,       -a.dx / (1. + a.x^2)                                     )
@DiffRule1(Base.asind,       180. / π / sqrt(1. - a.x^2) * a.dx                      )
@DiffRule1(Base.acosd,      -180. / π / sqrt(1. - a.x^2) * a.dx                      )
@DiffRule1(Base.atand,       180. / π / (1. + a.x^2) * a.dx                          )
@DiffRule1(Base.asecd,       180. / π / abs(a.x) / sqrt(a.x^2- 1.) * a.dx            )
@DiffRule1(Base.acscd,      -180. / π / abs(a.x) / sqrt(a.x^2- 1.) * a.dx            )
@DiffRule1(Base.acotd,      -180. / π / (1. + a.x^2) * a.dx                          )
@DiffRule1(Base.sinh,        cosh(a.x) * a.dx                                        )
@DiffRule1(Base.cosh,        sinh(a.x) * a.dx                                        )
@DiffRule1(Base.tanh,        sech(a.x)^2 * a.dx                                      )
@DiffRule1(Base.sech,       -tanh(a.x) * sech(a.x) * a.dx                            )
@DiffRule1(Base.csch,       -coth(a.x) * csch(a.x) * a.dx                            )
@DiffRule1(Base.coth,       -csch(a.x)^2                                             )
@DiffRule1(Base.asinh,       a.dx / sqrt(a.x^2 + 1.)                                 )
@DiffRule1(Base.acosh,       a.dx / sqrt(a.x^2 - 1.)                                 )
@DiffRule1(Base.atanh,       a.dx / (1. - a.x^2)                                     )
@DiffRule1(Base.asech,      -a.dx / a.x / sqrt(1. - a.x^2)                           )
@DiffRule1(Base.acsch,      -a.dx / abs(a.x) / sqrt(1. + a.x^2)                      )
@DiffRule1(Base.acoth,       a.dx / (1. - a.x^2)                                     )
@DiffRule1(SpecialFunctions.erf,         2. * exp(-a.x^2) / sqrt(π) * a.dx           )
@DiffRule1(SpecialFunctions.erfc,       -2. * exp(-a.x^2) / sqrt(π) * a.dx           )
@DiffRule1(SpecialFunctions.erfi,        2. * exp( a.x^2) / sqrt(π) * a.dx           )
@DiffRule1(SpecialFunctions.gamma,       digamma(a.x) * gamma(a.x) * a.dx            )
@DiffRule1(SpecialFunctions.lgamma,      digamma(a.x) * a.dx                         )
@DiffRule1(SpecialFunctions.airy,        airyprime(a.x) * a.dx                       )  # note: only covers the 1-arg version
@DiffRule1(SpecialFunctions.airyprime,   airy(2., a.x) * a.dx                        )
@DiffRule1(SpecialFunctions.airyai,      airyaiprime(a.x) * a.dx                     )
@DiffRule1(SpecialFunctions.airybi,      airybiprime(a.x) * a.dx                     )
@DiffRule1(SpecialFunctions.airyaiprime, a.x * airyai(a.x) * a.dx                    )
@DiffRule1(SpecialFunctions.airybiprime, a.x * airybi(a.x) * a.dx                    )
@DiffRule1(SpecialFunctions.besselj0,   -besselj1(a.x) * a.dx                        )
@DiffRule1(SpecialFunctions.besselj1,   (besselj0(a.x) - besselj(2., a.x))/2. * a.dx )
@DiffRule1(SpecialFunctions.bessely0,   -bessely1(a.x) * a.dx                        )
@DiffRule1(SpecialFunctions.bessely1,   (bessely0(a.x) - bessely(2., a.x))/2. * a.dx )



## Find NaN in derivatives
@inline hasnan(a::ℝ   )              = isnan(a)
@inline hasnan(a::∂ℝ   )             = hasnan(a.x) || hasnan(a.dx)
@inline hasnan(a::Tuple)             = any(hasnan.(a))
@inline hasnan(a::NamedTuple)        = any(hasnan.(values(a)))
@inline hasnan(a...;)                = any(hasnan.(a))
@inline hasnan(a)                    = false
#hasnan(a::AbstractArray)     = any(hasnan.(a)) # slow
@inline function hasnan(a::AbstractArray) 
    for aᵢ ∈ a
        if hasnan(aᵢ)
            return true
        end
    end
    return false
end

# Print addifs
const subscripts = ('₁','₂','₃','₄','₅','₆','₇','₈','₉')
string_(a::Float64) = strip(@sprintf("%4.2g",a))
function string_(a::∂ℝ{P,N,R}) where{P,N,R}
    p = subscripts[P]
    x = string_(a.x)
    dx = N==0 ? "" : @sprintf("%s",string_(a.dx[1]))
    for i = 2:N
       dx =   @sprintf("%s,%s",dx,string_(a.dx[i])) 
    end
    return @sprintf("%s+∂%s⟨%s⟩ ",x,p,dx) # \partial \langle \rangle mathematicaly-explicit
end
Base.show(io::IO,x::∂ℝ) = print(io,string_(x))

"""
    Muscade.isapprox_dbg(a,b;kwargs...)

Compare two `∂ℝ` of the same type (or two conformant `AbstractArrays` of same `eltype`).
The comparison applies to the value and the partial derivatives.

This differs from 

    isapprox(a,b;kwargs...)

which only compares values, in keeping with the principle that `∂ℝ` must behave as `ℝ`.    

""" 
isapprox_dbg(a::R,b::R;kwargs...) where{R<:∂ℝ} = isapprox(a.x,b.x;kwargs...)  && all(isapprox.(a.dx,b.dx;kwargs...))
isapprox_dbg(a::AbstractArray{R},b::AbstractArray{R};kwargs...) where{R<:∂ℝ} = all(isapprox.(a,b;kwargs...))
