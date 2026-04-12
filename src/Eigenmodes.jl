
# A pseudo random number generator that will always generate the same sequence of Float64 ∈ [0,1[, hopefully on any machine
mutable struct PRNG
    state::UInt64
end
PRNG() = PRNG(0x60d6a817f531f835)
Base.rand(prng::PRNG) = prng.state = 0x4820824402284023 * prng.state + 0x0000000000000001
Base.rand(::Type{𝕣},prng::PRNG) =   rand(prng)/0xffffffffffffffff
Base.rand(::Type{𝕔},prng::PRNG) = 𝕔(rand(prng)/0xffffffffffffffff,
                               rand(prng)/0xffffffffffffffff)
function Base.rand(::Type{T}, siz, prng) where{T} 
    out = Array{T,length(siz)}(undef,siz...)
    for i∈eachindex(out)
        out[i] = rand(T,prng)
    end
    return out
end

function normalize∞!(vec) # ensures that the largest term is 1 (with zero ℑ part)
    imax  = argmax(abs2(vecᵢ) for vecᵢ∈vec)
    vec ./= vec[imax]
end

"""
    λ,v,ncv = Muscade.geneig{ALGO}(A,B,neig=5)

Solves `(A-λ*B)*v=0`, finding the `neig` lowest eigenvalues `λ` (in absolute value) and the corresponding eigenvectors `v` (a `Vector{Vector}`)

# Input
`ALGO` can be
- :SDP       if `A` is symmetric definite positive and `B` is symmetric.  Will return real `λ` and `v`.
- :Hermitian if `A` is symmetric indefinite and `B` is symmetric. Will return real `λ` and `v`.
- :Complex   otherwise, will return complex `λ` and `v`. 
Optional keyword arguments:
- `maxiter     = 300`
- `verbosity   = 0 ∈ {0,1,2,3}`
- `krylovdim   = 2neig+6`

Uses KrylovKit.jl. Freely based on VibrationGEPHelpers.jl and input from PetrKryslUCSD and stevengj.  
See GIThub-blame for bug-credits.
"""
struct geneig{ALGO} end
function geneig{:complex}(luA::SparseArrays.UMFPACK.UmfpackLU,B,neig=5;maxiter=300,verbosity=0,seed=rand(𝕔,size(luA,1),PRNG()),normalize=true,kwargs...) 
    val, vec, info = eigsolve(x->B*(luA\x), seed,neig,:LR; maxiter,verbosity,ishermitian=false,kwargs...)
    for vecᵢ ∈ vec  
        vecᵢ .= luA\vecᵢ
    end
    val .= 1 ./val
    normalize && normalize∞!.(vec)

    return val, vec, info.converged
end
function geneig{:Hermitian}(luA::SparseArrays.UMFPACK.UmfpackLU,B,neig=5;kwargs...) 
    val, vec, info = geneig{:complex}(luA,B,neig;kwargs...) 
    return ℜ.(val), ℜ.(vec), info
end
function geneig{:symmetric}(luA::SparseArrays.UMFPACK.UmfpackLU,B,neig=5;maxiter=300,verbosity=0,seed=rand(𝕣,size(luA,1),PRNG()),normalize=true,kwargs...) 
    val, vec, info = eigsolve(x->B*(luA\x), seed,neig,:LR; maxiter,verbosity,issymmetric=false,kwargs...)
    for vecᵢ ∈ vec  
        vecᵢ .= luA\vecᵢ
    end
    val .= 1 ./val
    normalize && normalize∞!.(vec)
    return ℜ.(val), ℜ.(vec), info.converged
end
function geneig{:SDP}(L::SparseArrays.CHOLMOD.FactorComponent,B=I,neig=5;maxiter=300,verbosity=0,seed=rand(𝕣,size(L,1),PRNG()),normalize=true,kwargs...)
    val, vec, info = eigsolve(x->L\(B*(L'\x)),seed,neig,:LR; maxiter,verbosity,ishermitian=true,kwargs...)
    for vecᵢ ∈ vec
        vecᵢ .= ℜ.(L'\vecᵢ)
    end
    val .= 1 ./val
    normalize && normalize∞!.(vec)
    return val, vec, info.converged
end

geneig{:complex  }(A::SparseMatrixCSC,B,neig=5;kwargs...) = geneig{:complex  }(lu(A)                     ,B,neig;                     kwargs...)
geneig{:Hermitian}(A::SparseMatrixCSC,B,neig=5;kwargs...) = geneig{:Hermitian}(lu(A)                     ,B,neig;                     kwargs...)
geneig{:symmetric}(A::SparseMatrixCSC,B,neig=5;kwargs...) = geneig{:symmetric}(lu(A)                     ,B,neig;                     kwargs...)
geneig{:SDP      }(A::SparseMatrixCSC,B,neig=5;kwargs...) = geneig{:SDP      }(cholesky(Symmetric(A)).PtL,B,neig;                     kwargs...)
