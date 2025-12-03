
 #\circ \otimes

"""
    c = a⊗b

Compute the exterior product of two arrays, so that `cᵢⱼ=aᵢ bⱼ` where `i` and `j` can be multiple indices.

See also: [`∘₁`](@ref),[`∘₂`](@ref)
"""     
⊗(a,b) = dots(a,b,Val(0))
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


# Arrays, Views, etc. but not StaticArrays
@generated function dots(a::AbstractArray{Ta,Na},b::AbstractArray{Tb,Nb},::Val{ndot}) where {Ta,Tb,Na,Nb,ndot}
    # "Elrod", on Discourse/Julia, is acknowledged for the forerunner to this code
    Nar    = Na - ndot
    Nbr    = Nb - ndot
    Nc     = Nar+ Nbr
    Nloop  = Nc + ndot
    Tc     = promote_type(Ta,Tb)

    if Nc > 0  
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

# Accelerators (won't work on views, transpose etc.)
dots(a::Matrix{Ta},b::Matrix{Tb},::Val{1}) where {Ta<:Number,Tb<:Number      } = a*b
dots(a::Vector{Ta},b::Matrix{Tb},::Val{1}) where {Ta<:Number,Tb<:Number      } = b'*a
dots(a::Matrix{Ta},b::Vector{Tb},::Val{1}) where {Ta<:Number,Tb<:Number      } = a*b
dots(a::Vector{Ta},b::Vector{Tb},::Val{0}) where {Ta<:Number,Tb<:Number      } = a*b'

# StaticArrays
@generated function dots(a::StaticArray,b::StaticArray,::Val{ndots}) where{ndots}
    sa,sb       = size(a),size(b)
    ha,ta,hb,tb = sa[1:end-ndots], sa[end-ndots+1:end], sb[1:ndots], sb[ndots+1:end]
    lha,lta,ltb = prod(ha), prod(ta), prod(tb)
    so          = Tuple{ha...,tb...}
    zlo         = ntuple(i->0,lha*ltb)
    @assert length(sa)≥ndots
    @assert length(sb)≥ndots
    @assert ta==hb
    return if ndots==0                             quote SArray{$so}(SVector{$lha}(a)*SVector{$ltb}(b)')            end # external product
    elseif length(sa)==ndots && length(sb)==ndots  quote transpose(SVector{$lta}(a))*SVector{$lta}(b)               end # scalar product
    elseif lta==0                                  quote SArray{$so}($(zlo...))                                     end # sum over zero terms
    else                                           quote SArray{$so}(SMatrix{$lha,$lta}(a)*SMatrix{$lta,$ltb}(b))   end # general product with summation and array output
    end
end


# Scalar times array (it's a taste)
dots(a::              Ta    ,b::AbstractArray{Tb,Nb},::Val{0}) where {Ta<:Number,Tb<:Number   ,Nb} = a*b
dots(a::AbstractArray{Ta,Na},b::              Tb    ,::Val{0}) where {Ta<:Number,Tb<:Number,Na   } = a*b
dots(a::              Ta    ,b::              Tb    ,::Val{0}) where {Ta<:Number,Tb<:Number      } = a*b
