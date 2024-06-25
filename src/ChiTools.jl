# Tools for transforming χ data structures


# lossy conversion of an adiff to exactly type T
cast( ::Type{T}        ,a::T) where{T    } = a
cast(T::Type{∂ℝ{P,N,R}},a::𝕣) where{P,N,R} = ∂ℝ{P,N,R}(cast(R,a),SV{N,R}(zero(R) for j=1:N))
cast(T::Type{𝕣}        ,a::ℝ)              = VALUE(a)
function cast(T::Type{∂ℝ{PT,NT,RT}},a::∂ℝ{Pa,Na,Ra}) where{PT,NT,RT,Pa,Na,Ra}
    R = promote_type(RT,Ra)
    return if PT==Pa ∂ℝ{Pa,Na}(cast(RT,a.x),cast.(RT,a.dx)                )
    elseif    PT> Pa ∂ℝ{PT,NT}(cast(RT,a  ),SV{NT,RT}(zero(RT) for j=1:NT))
    else                       cast(T ,a.x)
    end
end

# lossy conversion of an adiff to at most type T (never introduce new time derivatives)
castdown( ::Type{T}        ,a::T) where{T    } = a
castdown(T::Type{∂ℝ{P,N,R}},a::𝕣) where{P,N,R} = a
castdown(T::Type{𝕣}        ,a::ℝ)              = VALUE(a)
function castdown(T::Type{∂ℝ{PT,NT,RT}},a::∂ℝ{Pa,Na,Ra}) where{PT,NT,RT,Pa,Na,Ra}
    R = promote_type(RT,Ra)
    return if PT==Pa ∂ℝ{Pa,Na}(castdown(RT,a.x),castdown.(RT,a.dx))
    elseif    PT> Pa           castdown(RT,a  )
    else                       castdown(T ,a.x)
    end
end

# lossless promotion to an adiff of at least type T 
castup( ::Type{T}        ,a::T) where{T    } = a
castup(T::Type{∂ℝ{P,N,R}},a::𝕣) where{P,N,R} = ∂ℝ{P,N,R}(castup(R,a),SV{N,R}(zero(R) for j=1:N))
castup(T::Type{𝕣}        ,a::ℝ)              = a
function castup(T::Type{∂ℝ{PT,NT,RT}},a::∂ℝ{Pa,Na,Ra}) where{PT,NT,RT,Pa,Na,Ra}
    R = promote_type(RT,Ra)
    return if PT==Pa ∂ℝ{Pa,Na}(castup(RT,a.x),castup.(RT,a.dx)              )
    elseif    PT> Pa ∂ℝ{PT,NT}(castup(RT,a  ),SV{NT,RT}(zero(RT) for j=1:NT))
    else             ∂ℝ{Pa,Na}(castup(T ,a.x),castup.(T,a.dx)               )
    end
end

# recursively apply a χ-cleaning function f to a data structure χ
χrecurse(f,χ::ℝ)             = f(χ)
χrecurse(f,χ::AbstractArray) = χrecurse.(f,χ)
χrecurse(f,χ::Tuple)         = χrecurse.(f,χ)
χrecurse(f,χ::NamedTuple)    = NamedTuple{keys(χ)}(χrecurse(f,values(χ)))
χrecurse(f,χ)                = χ  
χrecurse(f,χ::ℤ)             = χ  
χrecurse(f,χ::𝔹)             = χ  
