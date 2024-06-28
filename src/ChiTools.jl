# Tools for transforming χ data structures

# get as-meshed χ by asking all elements
χinit(model::Model)                           = χinit.(model.eleobj) # For each element type
χinit(e::Vector{E}) where{E<:AbstractElement} = χinit.(e)            # For each element in a type

# given χ (typicaly: as meshed), and a type (e.g. "store A-sensitivities), allocate memory for State.χ
χalloc(T∂,χ)                         =  [χalloc_(T∂,χᵢ) for χᵢ∈χ] 
χalloc_(T∂,χ::Vector{Tχ }) where{Tχ} = Vector{χcasttype(T∂,Tχ)}(undef,length(χ))

# lossy conversion of ℝ to T∂ in a χ-structure (incl. in Arrays)
χcast(  ::Type{T∂}       ,a::T∂) where{T∂   } = a
χcast(T∂::Type{∂ℝ{P,N,R}},a::𝕣 ) where{P,N,R} = ∂ℝ{P,N,R}(χcast(R,a),SV{N,R}(zero(R) for j=1:N))
χcast(T∂::Type{𝕣}        ,a::ℝ )              = VALUE(a)
χcast(T∂::Type{𝕣}        ,a::ℤ )              = a
function χcast(T∂::Type{∂ℝ{PT,NT,RT}},a::∂ℝ{Pa,Na,Ra}) where{PT,NT,RT,Pa,Na,Ra}
    R = promote_type(RT,Ra)
    return if PT==Pa ∂ℝ{Pa,Na}(χcast(RT,a.x),χcast.(RT,a.dx)                )
    elseif    PT> Pa ∂ℝ{PT,NT}(χcast(RT,a  ),SV{NT,RT}(zero(RT) for j=1:NT))
    else                       χcast(T∂ ,a.x)
    end
end
χcast(::Type{T∂},a::AbstractArray) where{T∂} = χcast.(T∂,a) # covers Array-over-elements (Allocates!), but also SArray within element
χcast(::Type{T∂},a::Tuple)         where{T∂} = χcast.(T∂,a)
χcast(::Type{T∂},a::NamedTuple)    where{T∂} = NamedTuple{keys(a)}(χcast.(T∂,values(a)))
χcast(::Type{T∂},a::ℤ)             where{T∂} = a  
χcast(::Type{T∂},a::𝔹)             where{T∂} = a  
χcast(::Type{T∂},a::Symbol)        where{T∂} = a  
χcast(::Type{T∂},a::Nothing)       where{T∂} = a  

# type of the lossy conversion of ℝ to T∂ in a χ-structure (excludes Arrays, but includes SArrays)
χcasttype(  ::Type{T∂}       ,::Type{T∂}) where{T∂   } = T∂
χcasttype(T∂::Type{∂ℝ{P,N,R}},::Type{𝕣 }) where{P,N,R} = ∂ℝ{P,N,R}
χcasttype(T∂::Type{𝕣}        ,::Type{R }) where{R<:ℝ}  = 𝕣       ####
χcasttype(T∂::Type{𝕣}        ,::Type{Z }) where{Z<:ℤ}  = Z       ####
function χcasttype(T∂::Type{∂ℝ{PT,NT,RT}},::Type{∂ℝ{Pa,Na,Ra}}) where{PT,NT,RT,Pa,Na,Ra}
    R = promote_type(RT,Ra)
    return if PT==Pa ∂ℝ{Pa,Na,χcasttype(RT,Ra)}
    elseif    PT> Pa ∂ℝ{PT,NT,χcasttype(RT,∂ℝ{Pa,Na,Ra})}
    else                      χcasttype(T∂,Ra)
    end
end
χcasttype(::Type{T∂},::Type{SVector{L,  T }}) where{T∂,L,  T} = SVector{L,  χcasttype(T∂,T)}
χcasttype(::Type{T∂},::Type{SMatrix{M,N,T }}) where{T∂,M,N,T} = SMatrix{M,N,χcasttype(T∂,T)}
χcasttype(::Type{T∂},::Type{SArray{ S,  T }}) where{T∂,S,  T} = SArray{ S,  χcasttype(T∂,T)}
χcasttype(::Type{T∂},::Type{V              }) where{T∂,V<:Tuple} = Tuple{χcasttype.(T∂,fieldtypes(V))...}
χcasttype(::Type{T∂},::Type{NamedTuple{K,V}}) where{T∂,K,V     } = NamedTuple{K,χcasttype.(T∂,V)}
χcasttype(::Type{T∂},::Type{Z})          where{T∂,Z<:ℤ} = Z    ####
χcasttype(::Type{T∂},::Type{B})          where{T∂,B<:𝔹} = B  
χcasttype(::Type{T∂},::Type{Symbol})     where{T∂}      = Symbol  
χcasttype(::Type{T∂},::Type{Nothing})    where{T∂}      = Nothing  



