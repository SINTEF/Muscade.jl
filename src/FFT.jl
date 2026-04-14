#https://en.wikipedia.org/wiki/Bit-reversal_permutation
function bitreversalpermutation(p) 
    n                           = 2^p
    brp                         = Vector{𝕫}(undef,n)
    brp[1]                      = 0
    ek                          = 1
    for k                       = 0:p-1  
        @simd for j             = 1:ek
            @inbounds brp[j]   *= 2
            @inbounds brp[j+ek] =  brp[j]+1
        end
        ek                     *=2
    end 
    @simd for j                 = 1:n # base-1 indexing of arrays in Julia
        @inbounds brp[j]       += 1
    end
    return brp
end
function applybrp!(a,brp)
    @assert length(a) == length(brp)
    @simd for i = 1:length(a)
        @inbounds j = brp[i]
        if j≤i
            @inbounds a[i],a[j] = a[j],a[i]
        end
    end
end
function getiW(nc) # = 𝑖 * expπ𝑖(-(0:nc-1)/nc) - the twiddles
    ωₘ         = expπ𝑖(-1/nc) 
    iW         = Vector{Complex{Float64}}(undef,nc)
    iW[1]      = 𝑖
    @simd for i      = 1:nc-1
        @inbounds iW[i+1]= iW[i]*ωₘ
        if mod(i,1024) == 15  # once in a while, renormalize iW 
            @inbounds iW[i+1] *= (3-ℜ(iW[i+1])^2-ℑ(iW[i+1])^2)/2  # but do so fast, do it using a 1st order Taylor expansion at 1 of of 1/|ω|
        end
    end
    return iW
end

# iterative radix-2 FFT algorithm implemented using bit-reversal permutation.
# https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm 
# Cormen, Thomas H.; Leiserson, Charles; Rivest, Ronald; Stein, Clifford (2009). 
# Introduction to algorithms (3rd ed.). Cambridge, Mass.: MIT Press. pp. 915–918. ISBN 978-0-262-03384-8
#
# A (complex, length 2^p): mutable memory for input and output
# p (integer): log2 of the length of A and a
# z (integer): -1 for forward and 1 for inverse transform
function basic_fft!(A::AbstractVector{Complex{R}},brp,p,z) where{R<:Real}   # Wikipedia convention
    # Assumes n == 2^p
    n                              = length(A)
    @assert 2^p == n
    applybrp!(A,brp)
    m                              = 1
    for s ∈ 0:p-1 # s-1 in Wikipedia          Scale of blocks m = 1,2,4,8,..,n/2
        ωₘ                         = expπ𝑖(z/m) 
        for k ∈ 1:2m:n # k+1 in Wikipedia     Block index
            ω                      = complex(1.)
            @simd for j ∈ 0:m-1              #      Within block 
                @inbounds t        = ω*A[k+j+m]  # if a∈ℝ, i>1 then A[i]=A[n-i+2]'
                @inbounds u        =   A[k+j  ]  # hence for m=n/2, A[n/2+k+j] = A[n/2+2-k-j]'
                @inbounds A[k+j]   = u+t
                @inbounds A[k+j+m] = u-t
                ω                 *= ωₘ   # TODO use a precomputed twiddle instead
                if mod(j,1024) == 15  # once in a while, renormalize ω 
                    ω             *= (3-ℜ(ω)^2-ℑ(ω)^2)/2  # but do so fast, do it using a 1st order Taylor expansion at 1 of of 1/|ω|
                end
            end
        end
        m                         *= 2
    end
end

## ℝ → ℂ transforms
# See "fft real.pdf" for theory
# forward 2^pc complex transform A of 2*2^pc real vector a
#   A = reinterpret(Complex{Float64},a) # just a view
#   basic_rfft!(A,brp,iW,pc)   # TODO reinterpret is not an AbstractVector
# A is the half Fourrier transform
function basic_rfft!(A::AbstractVector{Complex{R}},brp,iW,pc) where{R<:Real}
    nc     = 2^pc # N in the theory
    @assert length(A  )== nc
    @assert length(brp)== nc
    
    basic_fft!(A,brp,pc,-1)
    @inbounds A[1]      *= Complex(1,-1)
    @simd for i      = 1:(div(nc,2))-1
        j      = mod(nc-i,nc)   
        @inbounds α      = 1-iW[j+1]
        @inbounds β      = 1+iW[j+1]
        @inbounds Aᵢ     = (     A[i+1] *conj(α) + conj(A[j+1])*conj(β))/2
        @inbounds Aⱼ     = (conj(A[i+1])*     β  +      A[j+1] *     α )/2
        @inbounds A[i+1] = Aᵢ
        @inbounds A[j+1] = Aⱼ
    end
end
# A is the half Fourrier transform
#   basic_irfft!(A,brp,iW,pc)
#   a = reinterpret(Float64,A) # just a view
# A is the half Fourrier transform
function basic_irfft!(A::AbstractVector{Complex{R}},brp,iW,pc) where{R<:Real}
    nc     = 2^pc # N in the theory
    @assert length(A  )== nc
    @assert length(brp)== nc
    @inbounds A[1]      /= Complex(1,-1)
    @simd for i      = 1:(div(nc,2))-1
        j      = mod(nc-i,nc)   
        @inbounds Aᵢ     = A[i+1]
        @inbounds Aⱼ     = A[j+1]
        @inbounds α      = 1-iW[j+1]
        @inbounds β      = 1+iW[j+1]
        det    = α^2-β^2
        @inbounds A[i+1] = conj(( α*conj(Aᵢ) - β*Aⱼ)/det) *2
        @inbounds A[j+1] =     ((-β*conj(Aᵢ) + α*Aⱼ)/det) *2
    end
    basic_fft!(A,brp,pc,1)  
end
"""
    X = Muscade.𝔉(x,δt)  # typeset with \\mfrakF\\Bbbr

Fourrier transform of a real time series x stored at time steps `δt` and length `2N = 2*2^p`
into a complex spectre X stored at frequency intervals `δω=getδω(2N,δt)=2π/(2N*δt)`.  
The length of the spectre is `N`: only positive frequencies are stored (the Fourrier 
transform of real functions are Hermitian).

This provides a discretization of the unitary Fourrier transform, 
    
    G(ω) = 𝔉(g)(ω) = 1/√(2π) ∫exp(-𝑖ωt) g(t) dt

𝔉 is unitary, in the sense that

    sum(abs2.(x))*δt ≈ 2*(sum(abs2.(X)) - abs2.(X[1])/2)*δω 
    
(since the discrete spectre is provided for ω≥0, it contains only half the energy)

# Arguments
- `x` a vector of real numbers representing a time series.  Its length must be a power of two.
- `δt` the time step of the time series

# Example

```
X   = 𝔉(x,δt) 
δω  = getδω(length(x),δt)
x′  = 𝔉⁻¹(X,δω) # ≈ x
```

See also: [`𝔉⁻¹`](@ref), [`getδω`](@ref), [`getδt`](@ref),
"""
function 𝔉(a::AbstractVector{R},δt::ℝ) where{R<:Real} #\mfrakF
    nr      = length(a)
    nc      = div(nr,2)
    pc      = 𝕫log2(nc)
    A       = copy(reinterpret(Complex{R},a))
    iW      = getiW(nc)
    brp     = bitreversalpermutation(pc)
    basic_rfft!(A,brp,iW,pc)
    A     .*= δt/√(2π)   
    return A
end
# A = reinterpret(Complex{R},a)  # a[it,idof] length 2N
# 𝔉!(A,δt)                       # A[iω,idof] length N
function 𝔉!(A::AbstractMatrix{Complex{R}},δt::ℝ) where{R<:Real} #\mfrakF
    nc,ndof  = size(A)
    pc      = 𝕫log2(nc)
    iW      = getiW(nc)
    brp     = bitreversalpermutation(pc)
    Threads.@threads for idof = 1:ndof  
        basic_rfft!(view(A,:,idof),brp,iW,pc)
    end
    A     .*= δt/√(2π)   
    return A
end
"""
    x = Muscade.𝔉⁻¹(X,δω)  # typeset with \\mfrakF\\^-\\^1

See [`𝔉`](@ref)    

# Arguments
- `X` a vector of complex numbers representing one side of a spectra. Its length must be a power of two.
- `δω`, the angular frequency step of spectra

# Example

```
X   = 𝔉(x,δt) 
δω  = getδω(length(x),δt)
x′  = 𝔉⁻¹(X,δω) # ≈ x
```

See also: [`𝔉⁻¹`](@ref), [`getδω`](@ref), [`getδt`](@ref),

"""
function 𝔉⁻¹(A::AbstractVector{Complex{R}},δω::ℝ) where{R<:Real} #\mfrakF
    nc      = length(A)
    nr      = 2nc
    pc      = 𝕫log2(nc)
    a       = Vector{R}(undef,nr)
    iW      = getiW(nc)
    brp     = bitreversalpermutation(pc)
    a       = copy(A)
    basic_irfft!(a,brp,iW,pc)
    a     .*= δω*√(2/π)
    return reinterpret(R,a) 
end
# 𝔉⁻¹!(A,δω)                     # A[iω,idof] length N
# a = reinterpret(Complex{R},A)  # a[it,idof] length 2N
function 𝔉⁻¹!(A::AbstractMatrix{Complex{R}},δω::ℝ) where{R<:Real} #\mfrakF
    nc,ndof  = size(A)
    pc      = 𝕫log2(nc)
    iW      = getiW(nc)
    brp     = bitreversalpermutation(pc)
    Threads.@threads for idof = 1:ndof  
        basic_irfft!(view(A,:,idof),brp,iW,pc)
    end
    A     .*= δω*√(2/π) 
    return A
end

"""
    δω=Muscade.getδω(n,δt) = 2π/(n*δt)

    `n` is the length of the time series
"""
getδω(n,δt)    =  2π/(n*δt)
"""
    Muscade.getδt(n,δω) = 2π/(n*δω)
 
    `n` is the length of the time series

"""
getδt(n,δω) = 2π/(n*δω)
