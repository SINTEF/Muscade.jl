######## Takes a matrix of matrices (a pattern of blocks), and assembles it into a bigmat

# i index
# n length
# p ptr (using the colptr compression convention - the index of start of a swath within a vector)
#
# b block  - indexing blocks within pattern
# l local  - indexing entries inside a block
# g global - indexing the whole bigmat
#
# r row
# c col
# v non-zero value
#
# sparse.rowval → ilr
# sparse.colptr → pilr  
# ilv    = pilr[ilc]+ i-1
# irow   = ilr[ilv]  


"""
    bigmat,bigmatasm,bigvecasm,bigvecdis = Muscade.prepare(pattern)

Prepare for the assembly of sparse blocks into a large sparse matrix. 
`bigmat` is allocated, with the correct sparsity structure, but its `nzval` undef'ed.    
Where some blocks share the same sparsity structure, `blocks` in `pattern` can have `===` elements.

`pattern` is a `SparseMatrixCSC{<:SparseMatrixCSC}`, where empty blocks are structuraly zero

See also: [`addin!`](@ref)
""" 
function prepare(pattern::SparseMatrixCSC{SparseMatrixCSC{Tv,𝕫},𝕫}) where{Tv} 
    nbr,nbc                       = size(pattern)
    nbr>0 && nbc>0 || muscadeerror("must have size(pattern,i)>0 ∀i∈[1,2]")
    nlr                           = [-1 for ibr=1:nbr+1] # number of local rows in each block-row of pattern
    nlc                           = [-1 for ibc=1:nbc+1] # number of local cols in each block-col of pattern
    ngv                           = 0
    asm_igv                       = Vector{𝕫1}(undef,nnz(pattern))
    for ibc                       = 1:nbc
        for ibv                   = pattern.colptr[ibc]:pattern.colptr[ibc+1]-1
            ibr                   = pattern.rowval[ibv]
            block                 = pattern.nzval[ ibv]
            nlr[ibr+1]==-1 || nlr[ibr+1]==block.m || muscadeerror(@printf "block row %i of the pattern has blocks with different numbers of rows")
            nlc[ibc+1]==-1 || nlc[ibc+1]==block.n || muscadeerror(@printf "block column %i of the pattern has blocks with different numbers of columns")
            nlr[ibr+1]            = block.m
            nlc[ibc+1]            = block.n
            nlv                   = length(block.nzval)
            ngv                  += nlv
            asm_igv[ibv]          = 𝕫1(undef,nlv) # allocate: for each block in pattern, where in bigmat value to put block value 
        end
    end
    nlr[1]                        = 1
    nlc[1]                        = 1
    if any(nlr.==-1) || any(nlc.==-1) 
        for i = findall(nlr.==-1)
            @printf("    row    %i of the pattern has only empty blocks\n",i)
        end
        for i = findall(nlc.==-1)
            @printf("    column %i of the pattern has only empty blocks\n",i)
        end
        muscadeerror("invalid sparse-of-sparse pattern")
    end
    pgr                           = cumsum(nlr)  # pgr[ibr]→igr global row corresponding to the first local row of each block
    pgc                           = cumsum(nlc)  # pgc[ibc]→igc global column corresponding to the first local row of each block
    ngr                           = pgr[end]-1
    ngc                           = pgc[end]-1

    # create asm and global matrix (gv, aka nzval is undef in global matrix)
    pigr                          = 𝕫1(undef,ngc+1)        # aka global.colptr
    igr                           = 𝕫1(undef,ngv  )        # aka global.rowval
    gv                            = Vector{Tv}(undef,ngv)  # aka global.nzval
    asm                           = Vector{𝕫2}(undef,nbr)  # asm[ibc][ibr,ilc] → igv for a given block, and local column, where does the storage start?
    pigr[1]                       = 1
    igv                           = 1
    for ibc                       = 1:nbc                  # for each block column
        for ilc                   = 1:nlc[ibc+1]           # for each local column
            igc                   = pgc[ibc]-1 + ilc 
            for ibv               = pattern.colptr[ibc]:pattern.colptr[ibc+1]-1 # for each block row
                ibr               = pattern.rowval[ibv]
                block             = pattern.nzval[ ibv]
                pilr,ilr          = block.colptr, block.rowval
                for ilv           = pilr[ilc]:pilr[ilc+1]-1 # for each value in block row
                    igr[igv]      = pgr[ibr]-1 + ilr[ilv]   # row of corresponding global value
                    asm_igv[ibv][ilv] = igv                 # index of global value for block value and local value 
                    igv          += 1    
                end
            end
            pigr[igc+1]           = igv   
        end
    end
    bigmat    = SparseMatrixCSC(ngr,ngc,pigr,          igr,           gv     )
    bigmatasm = SparseMatrixCSC(nbr,nbc,pattern.colptr,pattern.rowval,asm_igv) # The assembler is a sparse, with same sparsity as pattern, which values are vectors telling where the values of a block go in the values of bigmat
    return    bigmat,bigmatasm,pgr,pgc    
end
"""
    Muscade.addin!(asm,global,block,ibr,ibc,factor=1.)

Add a sparse `block` into a large `out` sparse matrix, at block-row and -column `ibr` and `ibc`.  
   Use [`prepare`](@ref) to allocate memory for `global` and build the assembler `asm`.
""" 
function addin!(asm::SparseMatrixCSC{𝕫1,𝕫},out::SparseMatrixCSC{Tv,Ti},block::SparseMatrixCSC{Tv,Ti},ibr::𝕫,ibc::𝕫,factor=idmult) where{Tv,Ti<:Integer}
    # dichotomy to find ibv (index into asm.nzval)
    lo   = asm.colptr[ibc]         
    hi   = asm.colptr[ibc+1]-1
    local ibv # declare ibv so that ibv in the while scope survives the end of the loop 
    while true
        ibv  = div(lo+hi,2)
        aibr = asm.rowval[ibv]
        if ibr == aibr       break # correct ibv
        else
            hi ≤ lo && muscadeerror(@sprintf("BlockSparseAssembler pattern has no block [%i,%i]",ibr,ibc))
            if ibr > aibr    lo = ibv+1
            else             hi = ibv-1  
            end
        end
    end
    # addin the block
    aigv = asm.nzval[ibv]     
    for ilv ∈ eachindex(aigv)  
        out.nzval[aigv[ilv]] += block.nzval[ilv] * factor
    end
end


"""
    addin!(asm,outvec,blockvec,ibr)

Add a full `block` vector into a large `outvec` full vector.  at block-row `ibr`.
Use [`prepare`](@ref) to create `asm`.

See also: [`prepare`](@ref)
""" 
function addin!(pgr::𝕫1,out::AbstractVector{Tv},block::Vector{Tv},ibr::𝕫,factor=idmult) where{Tv}
    for (ilv,igv)∈enumerate(pgr[ibr]:pgr[ibr+1]-1) 
        out[igv] += block[ilv] * factor
    end
end

# disassemble a block from a big-vector
disblock(pgc::𝕫1,v::Vector,ibc::𝕫) = view(v,pgc[ibc]:(pgc[ibc+1]-1))


"""
    sparser!(S::SparseMatrixCSC,keep::Function)

Eliminate terms that do not satisfy a criteria from the storage of a sparse matrix.
`S` will be mutated, and the size of its internal storage modified.
`keep` is a `Function` which to an index into `S.nzval` associate `true` if storage is to be kept for this term
and `false` otherwise. Alternatively, `keep` can be a `Vector{Bool}`

# Examples

    sparser!([T],S,i->abs(S.nzval[i])>tol)
    sparser!([T],S,keep::Vector{Boolean}) 

If `T` is provided (initialised as `T = copy(S)`), then result is set in `T`, `S` is unchanged, other wise `S` is mutated in place.
The input `keep::Vector{Boolean}` must be of length `nnz(S)`.

    sparser!([S1,S2,...],rtol=1e-9)
    
Operates in place, reducing the sparses `S1`, `S2` etc... to a common sparsity pattern.    

!!! warning
    In the first example, the `keep` *function* accesses `S.nzval[i]`, and the term is then mutated by `sparser!`. 
    Any criteria requiring multiple access to *nzval* must build a `Vector` before calling `sparser!`.

!!! warning
    Note that `assemble!` computes the `nzval` of a sparse, assumning that its sparsity structure `colptr` and `rowval`
    is unchanged since sparse storage was allocated by `asmmat` in `prepare`.  In other words, if applying `sparser!`
    directly to a `sparse` returned by `assemble!`, `assemble!` can no longer be called for this matrix. 
"""
function sparser!(T::SparseMatrixCSC,S::SparseMatrixCSC,keep::Function) 
    # it is assumed that T has same matrix-shape as S.
    # works also with S===T
    resize!(T.nzval ,length(S.nzval )) 
    resize!(T.rowval,length(S.rowval))
    ndrop               = 0
    for icol  = 1:S.n
        colptr          = S.colptr[icol]
        T.colptr[icol]  = colptr-ndrop
        for inz         = colptr:S.colptr[icol+1]-1
            if keep(inz)
                T.nzval[ inz-ndrop] = S.nzval[ inz]
                T.rowval[inz-ndrop] = S.rowval[inz]
            else    
                ndrop  += 1
            end
        end
    end
    T.colptr[T.n+1] = S.colptr[S.n+1] - ndrop
    nnz                 = length(S.nzval)
    resize!(T.nzval ,nnz-ndrop) 
    resize!(T.rowval,nnz-ndrop)
end

sparser!(T::SparseMatrixCSC,S::SparseMatrixCSC,keep::Vector{Bool})  = sparser!(T,S,i->keep[i])
function sparser!(T::SparseMatrixCSC,S::SparseMatrixCSC,rtol=1e-9) 
    atol  = rtol*maximum(abs,S)
    sparser!(T,S,i->abs(S.nzval[i])≥atol)
end
sparser!(S::SparseMatrixCSC,args...) = sparser!(S,S,args...) 
function sparser!(S::AbstractVector{SMat},rtol=1e-9::𝕣) where{SMat<:SparseMatrixCSC}
    atol  = [rtol*maximum(abs,Sᵢ) for Sᵢ∈S]
    keep  = [any(abs(S[i].nzval[j]) >atol[i] for i∈eachindex(S)) for j∈1:nnz(S[1])]
    for Sᵢ ∈ S
        sparser!(Sᵢ,keep)
    end
    return keep
end
# Undocumented, untested, unused
function sparser(S::SparseMatrixCSC,args...)
    T = copy(S)
    sparser!(T,S,args...)
end



 
