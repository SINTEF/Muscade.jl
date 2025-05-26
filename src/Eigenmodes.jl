
"""
    λ,v,ncv = geneig{ALGO}(A,B,neig=5)

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
function geneig{:SDP}(L::SparseArrays.CHOLMOD.FactorComponent,B=I,neig=5;maxiter=300,verbosity=0,krylovdim=2neig+6,seed=rand(size(A,1)),kwargs...)
    val, vec, info = eigsolve(x->L\(B*(L'\x)),seed,neig,:LR; maxiter,verbosity,ishermitian=true,krylovdim,kwargs...)
    for vecᵢ ∈ vec
        vecᵢ .= ℜ.(L'\vecᵢ)
    end
    val .= 1 ./val
    normalize!.(vec)
    return val, vec, info.converged
end
function geneig{:Complex}(luA::SparseArrays.UMFPACK.UmfpackLU,B,neig=5;maxiter=300,verbosity=0,krylovdim=2neig+6,seed=rand(𝕔,size(luA,1)),kwargs...) 
    val, vec, info = eigsolve(x->B*(luA\x), seed,neig,:LR; maxiter,verbosity,ishermitian=false,krylovdim,kwargs...)
    for vecᵢ ∈ vec  
        vecᵢ .= luA\vecᵢ
    end
    val .= 1 ./val
    normalize!.(vec)
    return val, vec, info.converged
end
function geneig{:Hermitian}(luA::SparseArrays.UMFPACK.UmfpackLU,B,neig=5;kwargs...) 
    val, vec, info = geneig{:Complex}(luA,B,neig;kwargs...) 
    return ℜ.(val), ℜ.(vec), info
end
geneig{:Hermitian}(A::SparseMatrixCSC,B,neig=5;kwargs...) = geneig{:Hermitian}(lu(A)                     ,B,neig;seed=rand(size(A,1)),kwargs...)
geneig{:Complex  }(A::SparseMatrixCSC,B,neig=5;kwargs...) = geneig{:Complex  }(lu(A)                     ,B,neig;                     kwargs...)
geneig{:SDP      }(A::SparseMatrixCSC,B,neig=5;kwargs...) = geneig{:SDP      }(cholesky(Symmetric(A)).PtL,B,neig;                     kwargs...)
