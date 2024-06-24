using   StaticArrays
using   SpecialFunctions
using   Printf


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
∂ℝ{P,N  }(x::R ,dx::SV{N,R}) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(x   ,SV{N,R}(dx))
∂ℝ{P,N  }(x::R             ) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(x   ,SV{N,R}(zero(R)                 for j=1:N))
∂ℝ{P,N  }(x::R,i::ℤ        ) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(x   ,SV{N,R}(i==j ? one(R) : zero(R) for j=1:N))
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
    else      muscadeerror(printf("Cannot convert precedence ",Pb," to ",Pa))
    end
end

# Pack and unpack
precedence( ::Type{<:∂ℝ{P,N,R}}) where{P,N,R<:ℝ}          = P
npartial(   ::Type{<:∂ℝ{P,N,R}}) where{P,N,R<:ℝ}          = N
precedence( ::Type{<:ℝ})                                  = 0
npartial(   ::Type{<:ℝ})                                  = 0
precedence(a::SA)     = precedence(eltype(a))
npartial(  a::SA)     = npartial(eltype(a))
precedence(a::ℝ)      = precedence(typeof(a))
npartial(  a::ℝ)      = npartial(typeof(a))
constants(tup::Tuple) = constants(tup...) 
constants( a,args...) = max(constants(a),constants(args...))
constants( a)         = 1+precedence(a) 
constants( ::Nothing) = 0

# variate
struct δ{P,N,R}                end # need dum, because syntax δ{P,N,R}() collides with default constructor
struct variate{P,N}            end
struct directional{P,N}        end 
δ{P,N,R}(                          ) where{P,N,R<:ℝ} = SV{N,∂ℝ{P,N,R}}(∂ℝ{P,N  }(zero(R),i                                         ) for i=1:N)
δ{P,N,R}(               δa::SV{N,𝕣}) where{P,N,R<:ℝ} = SV{N,∂ℝ{P,N,R}}(∂ℝ{P,N,R}(zero(R),SV{N,R}(i==j ? δa[i]  : zero(R) for i=1:N)) for j=1:N)


#variate{P,N}(a::SV{N,R}            ) where{P,N,R<:ℝ} = SV{N,∂ℝ{P,N,R}}(∂ℝ{P,N  }(a[i]   ,i                                         ) for i=1:N)
variate{P,N}(a::SV{N,R}            ) where{P,N,R<:ℝ} = SV{N,∂ℝ{P,N,R}}(∂ℝ{P,N  }(a[i],i) for i=1:N)


variate{P,N}(a::SV{N,R},δa::SV{N,𝕣}) where{P,N,R<:ℝ} = SV{N,∂ℝ{P,N,R}}(∂ℝ{P,N,R}(a[j]   ,SV{N,R}(i==j ? R(δa[i])  : zero(R) for i=1:N)) for j=1:N)

variate{P}(a::R) where{P,R<:ℝ} =  ∂ℝ{P,1}(a,SV{1,R}(one(R)))
directional{P}(a::SV{N,R},δa::SV{N,R}) where{P,N,R<:ℝ} = SV{N,∂ℝ{P,1,R}}(∂ℝ{P,1}(a[i],SV{1,R}(δa[i])) for i=1:N)

# Analyse
VALUE(a::Nothing )                     =        nothing
VALUE(a::ℝ )                           =        a
VALUE(a::∂ℝ)                           = VALUE( a.x)
VALUE(a::SA)                           = VALUE.(a)

struct ∂{P,N}                  end 
struct value{P,N}              end
struct value_∂{P,N}            end

value{P}(a::∂ℝ{P,N,R}) where{P,N,R   } = a.x
value{P}(a::R        ) where{P  ,R<:ℝ} = a
value{P}(a::SA       ) where{P       } = value{P}.(a)

# ∂{P}(a) is handled as ∂{P,1}(a) and returns a scalar 
∂{P,N}(a::     ∂ℝ{P,N,R} ) where{  P,N,R   } = a.dx
∂{P,N}(a::            R  ) where{  P,N,R<:ℝ} = SV{  N,R}(zero(R)    for i=1:N      )
∂{P,N}(a::SV{M,∂ℝ{P,N,R}}) where{M,P,N,R   } = SM{M,N,R}(a[i].dx[j] for i=1:M,j∈1:N) # ∂(a,x)[i,j] = ∂a[i]/∂x[j]
∂{P,N}(a::SV{M,       R }) where{M,P,N,R   } = SM{M,N,R}(zero(R)    for i=1:M,j=1:N)
∂{P  }(a::            R  ) where{  P,  R<:ℝ} = zero(R)
∂{P  }(a::     ∂ℝ{P,1,R} ) where{  P,  R   } = a.dx[1]
∂{P  }(a::SV{N,∂ℝ{P,1,R}}) where{  P,N,R   } = SV{  N,R}(a[i].dx[1] for i=1:N     ) # ∂(a,x)[i]    = ∂a[i]/∂x
#∂{P,N}(a::SA{M,∂ℝ{P,N,R}}) where{M,P,N,R}  = SA{(M...,N),R}(a[i].dx[j] for i∈eachindex(a),j∈1:N) # ∂(a,x)[i,...,j] = ∂a[i,...]/∂x[j]
#∂{P,N}(a::SA{M,       R }) where{M,P,N,R}  = SA{(M...,N),R}(zero(R)    for i∈eachindex(a),j∈1:N)

value_∂{P,N}(a) where{  P,N}= value{P}(a),∂{P,N}(a)
value_∂{P  }(a) where{  P  }= value{P}(a),∂{P  }(a)

## Binary operations
for OP∈(:(>),:(<),:(==),:(>=),:(<=),:(!=))
    @eval Base.$OP(a::∂ℝ,b::∂ℝ)  = $OP(VALUE(a),VALUE(b))
    @eval Base.$OP(a:: ℝ,b::∂ℝ)  = $OP(      a ,VALUE(b))
    @eval Base.$OP(a::∂ℝ,b:: ℝ)  = $OP(VALUE(a),      b )
end

macro Op2(OP,AB,A,B)
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

@Op2(Base.atan,  (a.dx*b.x+b.dx*a.x)/(a.x^2+b.x^2),          (a.dx*b)/(a.x^2+b^2),  (b.dx*a)/(a^2+b.x^2) )   
@Op2(Base.hypot, (a.dx*a.x+b.dx*b.x)/hypot(a.x,b.x),         a.dx*a.x/hypot(a.x,b), b.dx*b.x/hypot(a,b.x))   
@Op2(Base.:(+),  a.dx+b.dx,                                  a.dx,                  b.dx                 )
@Op2(Base.:(-),  a.dx-b.dx,                                  a.dx,                  -b.dx                )
@Op2(Base.:(*),  a.dx*b.x+a.x*b.dx,                          a.dx*b,                a*b.dx               )
@Op2(Base.:(/),  a.dx/b.x-a.x/b.x^2*b.dx,                    a.dx/b,                -a/b.x^2*b.dx        ) 
@Op2(Base.:(^),  a.dx*b.x*a.x^(b.x-1)+log(a.x)*a.x^b.x*b.dx, a.dx*b*a.x^(b  -1),    log(a)*a ^b.x*b.dx   )
@inline Base.:(^)(a::∂ℝ{P,N,R},b::Integer) where{P,N,R<:ℝ} = ∂ℝ{P,N,R}(a.x^b ,a.dx*b*a.x^(b-1) )

## Functions
macro Op1(OP,A)
    return esc(:(@inline $OP(a::∂ℝ{P,N}) where{P,N} = ∂ℝ{P,N}($OP(a.x),$A)))
end
@Op1(Base.:(+),       a.dx                                                     )
@Op1(Base.:(-),      -a.dx                                                     )
@Op1(Base.abs  ,a.x==0.0 ? zero(a.dx) : (a.x>0.0 ? a.dx : -a.dx)               )
@Op1(Base.conj ,      a.dx                                                     )
@Op1(Base.sqrt,       a.dx / 2. / sqrt(a.x)                                    )
@Op1(Base.cbrt,       a.dx / 3. / cbrt(a.x)^2                                  )
@Op1(Base.abs2,       a.dx*2. * a.x                                            )
@Op1(Base.inv,       -a.dx * abs2(inv(a.x))                                    )
@Op1(Base.log,        a.dx / a.x                                               )
@Op1(Base.log10,      a.dx / a.x / log(10.)                                    )
@Op1(Base.log2,       a.dx / a.x / log(2.)                                     )
@Op1(Base.log1p,      a.dx / (a.x + 1.)                                        )
@Op1(Base.exp,         exp(a.x) * a.dx                                         )
@Op1(Base.exp2,        log(2. ) * exp2( a.x) * a.dx                            )
@Op1(Base.exp10,       log(10.) * exp10(a.x) * a.dx                            )
@Op1(Base.expm1,       exp(a.x) * a.dx                                         )
@Op1(Base.sin,         cos(a.x) * a.dx                                         )
@Op1(Base.cos,        -sin(a.x) * a.dx                                         )
@Op1(Base.tan,         (1. + tan(a.x)^2) * a.dx                                )
@Op1(Base.sinpi,       π*cos(a.x) * a.dx                                       )
@Op1(Base.cospi,      -π*sin(a.x) * a.dx                                       )
@Op1(Base.sec,         sec(a.x) * tan(a.x) * a.dx                              )
@Op1(Base.csc,        -csc(a.x) * cot(a.x) * a.dx                              )
@Op1(Base.cot,        -(1. + cot(a.x)^2) * a.dx                                )
@Op1(Base.sind,        π / 180. * cosd(a.x) * a.dx                             )
@Op1(Base.cosd,       -π / 180. * sind(a.x) * a.dx                             )
@Op1(Base.tand,        π / 180. * (1. + tand(a.x)^2) * a.dx                    )
@Op1(Base.secd,        π / 180. * secd(a.x) * tand(a.x) * a.dx                 )
@Op1(Base.cscd,       -π / 180. * cscd(a.x) * cotd(a.x) * a.dx                 )
@Op1(Base.cotd,       -π / 180. * (1. + cotd(a.x)^2)  * a.dx                   )
@Op1(Base.asin,        a.dx / sqrt(1. - a.x^2)                                 )
@Op1(Base.acos,       -a.dx / sqrt(1. - a.x^2)                                 )
@Op1(Base.atan,        a.dx / (1. + a.x^2)                                     )
@Op1(Base.asec,        a.dx / abs(a.x) / sqrt(a.x^2 - 1.)                      )
@Op1(Base.acsc,       -a.dx / abs(a.x) / sqrt(a.x^2 - 1.)                      )
@Op1(Base.acot,       -a.dx / (1. + a.x^2)                                     )
@Op1(Base.asind,       180. / π / sqrt(1. - a.x^2) * a.dx                      )
@Op1(Base.acosd,      -180. / π / sqrt(1. - a.x^2) * a.dx                      )
@Op1(Base.atand,       180. / π / (1. + a.x^2) * a.dx                          )
@Op1(Base.asecd,       180. / π / abs(a.x) / sqrt(a.x^2- 1.) * a.dx            )
@Op1(Base.acscd,      -180. / π / abs(a.x) / sqrt(a.x^2- 1.) * a.dx            )
@Op1(Base.acotd,      -180. / π / (1. + a.x^2) * a.dx                          )
@Op1(Base.sinh,        cosh(a.x) * a.dx                                        )
@Op1(Base.cosh,        sinh(a.x) * a.dx                                        )
@Op1(Base.tanh,        sech(a.x)^2 * a.dx                                      )
@Op1(Base.sech,       -tanh(a.x) * sech(a.x) * a.dx                            )
@Op1(Base.csch,       -coth(a.x) * csch(a.x) * a.dx                            )
@Op1(Base.coth,       -csch(a.x)^2                                             )
@Op1(Base.asinh,       a.dx / sqrt(a.x^2 + 1.)                                 )
@Op1(Base.acosh,       a.dx / sqrt(a.x^2 - 1.)                                 )
@Op1(Base.atanh,       a.dx / (1. - a.x^2)                                     )
@Op1(Base.asech,      -a.dx / a.x / sqrt(1. - a.x^2)                           )
@Op1(Base.acsch,      -a.dx / abs(a.x) / sqrt(1. + a.x^2)                      )
@Op1(Base.acoth,       a.dx / (1. - a.x^2)                                     )
@Op1(SpecialFunctions.erf,         2. * exp(-a.x^2) / sqrt(π) * a.dx           )
@Op1(SpecialFunctions.erfc,       -2. * exp(-a.x^2) / sqrt(π) * a.dx           )
@Op1(SpecialFunctions.erfi,        2. * exp( a.x^2) / sqrt(π) * a.dx           )
@Op1(SpecialFunctions.gamma,       digamma(a.x) * gamma(a.x) * a.dx            )
@Op1(SpecialFunctions.lgamma,      digamma(a.x) * a.dx                         )
@Op1(SpecialFunctions.airy,        airyprime(a.x) * a.dx                       )  # note: only covers the 1-arg version
@Op1(SpecialFunctions.airyprime,   airy(2., a.x) * a.dx                        )
@Op1(SpecialFunctions.airyai,      airyaiprime(a.x) * a.dx                     )
@Op1(SpecialFunctions.airybi,      airybiprime(a.x) * a.dx                     )
@Op1(SpecialFunctions.airyaiprime, a.x * airyai(a.x) * a.dx                    )
@Op1(SpecialFunctions.airybiprime, a.x * airybi(a.x) * a.dx                    )
@Op1(SpecialFunctions.besselj0,   -besselj1(a.x) * a.dx                        )
@Op1(SpecialFunctions.besselj1,   (besselj0(a.x) - besselj(2., a.x))/2. * a.dx )
@Op1(SpecialFunctions.bessely0,   -bessely1(a.x) * a.dx                        )
@Op1(SpecialFunctions.bessely1,   (bessely0(a.x) - bessely(2., a.x))/2. * a.dx )

## Comparison for debug purposes
≗(a::ℝ,b::ℝ)                 = (typeof(a)==typeof(b)) && ((a-b) < 1e-10*max(1,a+b))
≗(a::SA,b::SA)               = (size(a)==size(b)) && all(a .≗ b)
≗(a::∂ℝ,b::∂ℝ)               = (typeof(a)==typeof(b)) && (a.x ≗ b.x) && (a.dx ≗ b.dx)

## Find NaN in derivatives
hasnan(a::ℝ   )              = isnan(a)
hasnan(a::∂ℝ   )             = hasnan(a.x) || hasnan(a.dx)
hasnan(a::AbstractArray)     = any(hasnan.(a))
hasnan(a::Tuple)             = any(hasnan.(a))
hasnan(a::NamedTuple)        = any(hasnan.(values(a)))
hasnan(a...;)                = any(hasnan.(a))
hasnan(a)                    = false

# cast: like `convert` but never throws an `inexact error` - and indeed willfully looses data if so asked
cast( ::Type{T}        ,a::T) where{T    } = a
cast(T::Type{∂ℝ{P,N,R}},a::𝕣) where{P,N,R} = ∂ℝ{P,N,R}(cast(R,a),SV{N,R}(zero(R) for j=1:N))
cast(T::Type{𝕣}        ,a::ℝ)              = VALUE(a)
function cast(T::Type{∂ℝ{PT,NT,RT}},a::∂ℝ{Pa,Na,Ra}) where{PT,NT,RT,Pa,Na,Ra}
    R = promote_type(RT,Ra)
    return if PT==Pa   ∂ℝ{Pa,Na}(cast(R ,a.x),cast.(R,a.dx)                 )
    elseif    PT> Pa   ∂ℝ{PT,NT}(cast(RT,a  ),SV{NT,RT}(zero(RT) for j=1:NT))
    else                         cast(T ,a.x)
    end
end

