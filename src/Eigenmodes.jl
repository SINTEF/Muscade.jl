using KrylovKit: KrylovKit


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
function geneig{:SDP}(A,B=I,neig=5;maxiter=300,verbosity=0,krylovdim=2neig+6,seed=rand(size(A,1)))
    L = cholesky(Symmetric(A)).PtL
    vals, vecs, info = KrylovKit.eigsolve(x->L\(B*(L'\x)),               seed,neig,:LR; maxiter,verbosity,ishermitian=true,krylovdim)
    for vecsᵢ ∈ vecs
        vecsᵢ .= ℜ.(L'\vecsᵢ)
    end
    vals .= 1 ./vals
    vecs = normalize.(vecs)
    ix   = sortperm(abs.(vals))
    return vals[ix], vecs[ix], info.converged
end
function geneig{ALGO}(A,B,neig=5;maxiter=300,verbosity=0,krylovdim=2neig+6,seed=rand(size(A,1))) where{ALGO}
    luK = lu(A)
    L,U,s,p,q = luK.L,luK.U,luK.Rs,luK.p,luK.q
    vals, vecs, info = KrylovKit.eigsolve(x->L\((s.*(B*((U\x)[q])))[p]), seed,neig,:LR; maxiter,verbosity,ishermitian=ALGO==:Hermitian,krylovdim)
    for vecsᵢ ∈ vecs  
        vecsᵢ .= (U\vecsᵢ)[q]
    end
    vals .= 1 ./vals
    normalize!.(vecs)
    ix   = sortperm(abs.(vals))
    return vals[ix], vecs[ix], info.converged
end




