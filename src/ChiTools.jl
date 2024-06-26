# Tools for transforming χ data structures

χinit(model::Model)                           = χinit.(model.eleobj) # For each element type
χinit(e::Vector{E}) where{E<:AbstractElement} = χinit.(e)            # For each element in a type

χalloc(T∂,χ) = [χalloc_(T∂,χᵢ) for χᵢ∈χ] # for each element type
function χalloc_(T∂,χᵢ) # for one element type
    nel = length(χᵢ)
    Tχ  = casttype(T∂,first(χᵢ))
    return Vector{Tχ}(undef,nel)  
end

# lossy conversion of ℝ to T in a χ-structure
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
cast(::Type{T},a::AbstractArray) where{T} = cast.(T,a)
cast(::Type{T},a::Tuple)         where{T} = cast.(T,a)
cast(::Type{T},a::NamedTuple)    where{T} = NamedTuple{keys(a)}(cast.(T,values(a)))
cast(::Type{T},a::ℤ)             where{T} = a  
cast(::Type{T},a::𝔹)             where{T} = a  
cast(::Type{T},a::Symbol)        where{T} = a  

# type of the lossy conversion of an adiff to exactly type T
casttype( ::Type{T}        ,::Type{T}) where{T    } = T
casttype(T::Type{∂ℝ{P,N,R}},::Type{𝕣}) where{P,N,R} = ∂ℝ{P,N,R}
casttype(T::Type{𝕣}        ,::Type{R}) where{R<:ℝ}  = 𝕣
function casttype(T::Type{∂ℝ{PT,NT,RT}},::Type{∂ℝ{Pa,Na,Ra}}) where{PT,NT,RT,Pa,Na,Ra}
    R = promote_type(RT,Ra)
    return if PT==Pa ∂ℝ{Pa,Na,casttype(RT,Ra)}
    elseif    PT> Pa ∂ℝ{PT,NT,casttype(RT,∂ℝ{Pa,Na,Ra}  )}
    else                      casttype(T ,Ra)
    end
end
casttype(::Type{T},::Type{V }) where{T,V<:Tuple               } =      Tuple{        casttype.(T,fieldtypes(V))... }
casttype(::Type{T},::Type{NT}) where{T,K,V,NT<:NamedTuple{K,V}} = NamedTuple{K,Tuple{casttype.(T,fieldtypes(V))...}}
casttype(::Type{T},::Type{Z})          where{T,Z<:ℤ} = Z  
casttype(::Type{T},::Type{B})          where{T,B<:𝔹} = B  
casttype(::Type{T},::Type{Symbol})     where{T}      = Symbol  



