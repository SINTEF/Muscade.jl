
 #\circ \otimes

∘₀(a,b) = dots(a,b,Val(0))
"""
    c = a∘₁b

Compute the single-dot product of two arrays, so that `cᵢⱼ=Σₖ aᵢₖ bₖⱼ` where `i` and `j` can be multiple indices.

See also: [`⊗`](@ref),[`∘₂`](@ref)
"""     
∘₁(a,b) = dots(a,b,Val(1))
"""
    c = a∘₂b

Compute the double-dot product of two arrays, so that `cᵢⱼ=Σₖₗ aᵢₖₗ bₖₗⱼ` where `i` and `j` can be multiple indices.

See also: [`∘₁`](@ref),[`⊗`](@ref)
"""     
∘₂(a,b) = dots(a,b,Val(2))

"""
    c = a⊗b

Compute the exterior product of two arrays, so that `cᵢⱼ=aᵢ bⱼ` where `i` and `j` can be multiple indices.

See also: [`∘₁`](@ref),[`∘₂`](@ref)
"""     
⊗(a,b) = dots(a,b,Val(0))

@generated function dots(a::AbstractArray{Ta,Na},b::AbstractArray{Tb,Nb},::Val{ndot}) where {Ta,Tb,Na,Nb,ndot}
    # "Elrod", on Discourse/Julia, is acknowledged for the forerunner to this code
    Nar    = Na - ndot
    Nbr    = Nb - ndot
    Nc     = Nar+ Nbr
    Nloop  = Nc + ndot
    Tc     = promote_type(Ta,Tb)

    if Nc > 0  # TODO: output is an Array, not StaticArray
        return quote
            sc = Base.@ntuple $Nc i -> begin
                i <= $Nar ? size(a)[i] : size(b)[i - $Nar + $ndot]
            end
            c  = zeros($Tc, sc)
            @nloops $Nloop i j -> begin
                j <= $Na ? axes(a,j) : axes(b,j - $Na + $ndot)
            end begin
                (@nref $Nc c j -> j <= $Nar ? i_{j} : i_{j+$ndot} ) += (@nref $Na a i) * (@nref $Nb b j -> i_{j+$Nar})
            end
            return c
        end
    else # scalar products
        Tc = promote_type(Ta,Tb)
        return quote
            c  = zero($Tc)
            @nloops $Na i a begin
                c += (@nref $Na a i) * (@nref $Na b i)
            end
            return c
        end
    end
end
dots(a::              Ta    ,b::AbstractArray{Tb,Nb},::Val{0}) where {Ta<:Number,Tb<:Number   ,Nb} = a*b
dots(a::AbstractArray{Ta,Na},b::              Tb    ,::Val{0}) where {Ta<:Number,Tb<:Number,Na   } = a*b
dots(a::              Ta    ,b::              Tb    ,::Val{0}) where {Ta<:Number,Tb<:Number      } = a*b

dots(a::SArray{Tuple{M1,M2,M3   }},b::SVector{M3},::Val{1}) where{M1,M2,M3   } = 
      SMatrix{M1,M2}(sum(a[i,j,k]*b[k] for k∈1:M3) for i∈1:M1,j∈1:M2)
dots(a::SArray{Tuple{M1,M2,M3,M4}},b::SVector{M4},::Val{1}) where{M1,M2,M3,M4} = 
      SArray{Tuple{M1,M2,M3}}(sum(a[i,j,k,ℓ]*b[ℓ] for ℓ∈1:M4) for i∈1:M1,j∈1:M2,k∈1:M3)

# Accelerators
dots(a::AbstractMatrix{Ta},b::AbstractMatrix{Tb},::Val{1}) where {Ta<:Number,Tb<:Number      } = a*b
dots(a::AbstractVector{Ta},b::AbstractMatrix{Tb},::Val{1}) where {Ta<:Number,Tb<:Number      } = b'*a
dots(a::AbstractMatrix{Ta},b::AbstractVector{Tb},::Val{1}) where {Ta<:Number,Tb<:Number      } = a*b
dots(a::AbstractVector{Ta},b::AbstractVector{Tb},::Val{0}) where {Ta<:Number,Tb<:Number      } = a*b'
