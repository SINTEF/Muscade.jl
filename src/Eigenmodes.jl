using KrylovKit: KrylovKit,eigsolve


"""
    λ,v,ncv = geneig{ALGO}(A,B,neig=5)

Solves `(A-λ*B)*v=0`, finding the `neig` lowest eigenvalues `λ` (in absolute value) and the corresponding eigenvectors `v` (a `Vector{Vector}`)

# Input
`ALGO` can be
- :SDP       if `A` is symmetric definite positive and `B` is symmetric.  Will return real `λ` and `v`.
- :Hermitian if `A` is symmetric indefinite and `B` is symmetric. Will return real `λ` and `v`.
- :Complex   otherwise
Optional keyword arguments:
- `maxiter     = 300`
- `verbosity   = 0 ∈ {0,1,2,3}`
- `krylovdim   = 2neig+6`

Uses KrylovKit.jl. Freely based on VibrationGEPHelpers.jl and input from PetrKryslUCSD and stevengj.  
See GIThub-blame for bug-credits.
"""
struct geneig{ALGO} end
function geneig{:SDP}(A::SparseMatrixCSC,B=I,neig=5;maxiter=300,verbosity=0,krylovdim=2neig+6,seed=rand(size(A,1)),kwargs...)
    L = cholesky(Symmetric(A)).PtL
    val, vec, info = eigsolve(x->L\(B*(L'\x)),seed,neig,:LR; maxiter,verbosity,ishermitian=true,krylovdim,kwargs...)
    info.converged ≥ neig || @sprintf("eigensolver only converged for %i out of %i modes",info.converged,neig)
    val = val[1:neig]
    vec = vec[1:neig]
    for vecᵢ ∈ vec
        vecᵢ .= ℜ.(L'\vecᵢ)
    end
    val .= 1 ./val
    normalize!.(vec)
    ix   = sortperm(abs.(val))
    return val[ix], vec[ix], info.converged
end
function geneig{:Complex}(A::SparseMatrixCSC,B,neig=5;maxiter=300,verbosity=0,krylovdim=2neig+6,seed=rand(𝕔,size(A,1)),kwargs...) 
    luA = lu(A)
    val, vec, info = eigsolve(x->B*(luA\x), seed,neig,:LR; maxiter,verbosity,ishermitian=false,krylovdim,kwargs...)
    info.converged ≥ neig || @sprintf("eigensolver only converged for %i out of %i modes",info.converged,neig)
    val = val[1:neig]
    vec = vec[1:neig]
    for vecᵢ ∈ vec  
        vecᵢ .= luA\vecᵢ
    end
    val .= 1 ./val
    normalize!.(vec)
    ix   = sortperm(abs.(val))
    return val[ix], vec[ix], info.converged
end
function geneig{:Hermitian}(A::SparseMatrixCSC,B,neig=5;seed=rand(size(A,1)),kwargs...) 
    val, vec, info = geneig{:Complex}(A,B,neig;seed,kwargs...) 
    return ℜ.(val), ℜ.(vec), info
end
