######## Takes a matrix of matrices (a pattern of blocks), and assembles it into a bigsparse

# i index
# n length
# p ptr (using the colptr compression convention - the index of start of a swath within a vector)
#
# b block  - indexing blocks within pattern
# l local  - indexing entries inside a block
# g global - indexing the whole bigsparse
#
# r row
# c col
# v non-zero value
#
# sparse.rowval → ilr
# sparse.colptr → pilr  
# ilv    = pilr[ilc]+ i-1
# irow   = ilr[ilv]  

struct BlockSparseAssembler
    pigr::Vector{𝕫2}   # pigr[ibc][ibr,ilc] → igv   Variable built up as "asm" in "prepare".  Given, index into block column & row and index into local colum, get index into gobal nz
    pgc ::𝕫1           # pgc[ibc]           → igc   Given index of block column, get index of first global column of the block 
end

Base.zero(SparseMatrixCSC) = nothing # so that, for a sparse pattern, indexing the pattern at a a structuraly zero block returns nothing
# provide same syntax for indexing into the pattern (full or sparse), returning `nothing` for empty blocks 
getblock(pattern::SparseMatrixCSC{SparseMatrixCSC{Tv,𝕫}},row,col) where{Tv} = pattern[row,col]
getblock(pattern::Matrix{SparseMatrixCSC{Tv,𝕫}},row,col) where{Tv} = isassigned(pattern,row,col) ? pattern[row,col] : nothing

"""
    bigsparse,asm = prepare(pattern)

Prepare for the assembly of sparse blocks into a large sparse matrix. Unassigned elements in `pattern` (undef'ed) will be treated as all-zero.
`bigsparse` is allocated, with the correct sparsity structure, but its `nzval` undef'ed.    
Where some blocks share the same sparsity structure, `blocks` can have `===` elements.

`pattern` can be a `Matrix{<:SparseMatrixCSC}`, where empty blocks are unassigned.
`pattern` can be a `SparseMatrixCSC{<:SparseMatrixCSC}`, where empty blocks are structuraly zero

See also: [`addin!`](@ref)
""" 
function prepare(pattern::AbstractMatrix{SparseMatrixCSC{Tv,𝕫}}) where{Tv} 
    nbr,nbc                 = size(pattern)  
    nbr>0 && nbc>0 || muscadeerror("must have length(pattern)>0")

    # determine the number rows in each row of blocks, store in pgr
    pgr                     = 𝕫1(undef,nbr+1)         # pgr[ibr]→igr pointers to the start of each block in global solution vector, where global*solution=rhs
    pgr[1]                  = 1
    for ibr                 = 1:nbr
        nlr                 = 0
        for ibc             = 1:nbc
            b = getblock(pattern,ibr,ibc)
            if ~isnothing(b)
                nlr         = size(b,1)
                break
            end
            ibc<nbc || muscadeerror("pattern has an empty row")
        end
        nlr > 0 || muscadeerror("every block-row must contain at least one assigned Sparse")
        pgr[ibr+1]          = pgr[ibr]+nlr
    end 
    ngr                     = pgr[end]-1

    # determine the number columns in each column of blocks, store in pgc
    pgc                     = 𝕫1(undef,nbc+1)         # pgc[ibc]→igc pointers to the start of each block in global rhs vector
    pgc[1]                  = 1
    for ibc                 = 1:nbc
        nlc                 = 0
        for ibr             = 1:nbr
            b = getblock(pattern,ibr,ibc)
            if ~isnothing(b)
                nlc         = size(b,2)
                break
            end
            ibr<nbr || muscadeerror("pattern has an empty column")
        end
        nlc > 0 || muscadeerror("every block-column must contain at least one assigned Sparse")
        pgc[ibc+1]          = pgc[ibc]+nlc
    end 
    ngc                     = pgc[end]-1

    # allocate arrays for global matrix
    ngv                     = 0                            
    for ibc                 = 1:nbc
        for ibr             = 1:nbr
            b = getblock(pattern,ibr,ibc)
            if ~isnothing(b)
                ngv        += nnz(b)
            end
        end
    end
    pigr                    = 𝕫1(undef,ngc+1)       # aka colptr
    igr                     = 𝕫1(undef,ngv  )       # aka rowval
    gv                      = Vector{Tv}(undef,ngv) # aka nzval

    # create asm and global matrix (gv, aka nzval is undef in global matrix)
    asm                               = Vector{𝕫2}(undef,nbr)  # asm[ibc][ibr,ilc] → igv for a given block, and local column, where does the storage start?
    pigr[1]                           = 1
    igv                               = 1
    for ibc                           = 1:nbc
       nlc                            = pgc[ibc+1]-pgc[ibc]
       asm[ibc]                       = zeros(𝕫,nbr,nlc)
       for ilc                        = 1:nlc
            igc                       = ilc + pgc[ibc]-1
            for ibr                   = 1:nbr
                b                     = getblock(pattern,ibr,ibc)
                if ~isnothing(b)
                    asm[ibc][ibr,ilc] = igv
                    pilr,ilr          = b.colptr, b.rowval
                    for ilv           = pilr[ilc]:pilr[ilc+1]-1 
                        igr[igv]      = pgr[ibr]-1 + ilr[ilv]  
                        igv          += 1    
                    end
                end
            end
            pigr[igc+1] = igv   
        end
    end
    return SparseMatrixCSC(ngr,ngc,pigr,igr,gv),BlockSparseAssembler(asm,pgc)    
end
"""
    addin!(asm,bigsparse,block::SparseMatrixCSC,ibr,ibc,factor=1.)

Add a sparse into one of the blocks of a large sparse matrix.  Will fail silently or throw an error unless
`bigsparse` has the correct sparsity structure for the given `blocks`. Use [`prepare`](@ref) to
    create `bigsparse` and `asm`.
""" 
function addin!(asm::BlockSparseAssembler,out::SparseMatrixCSC{Tv,Ti},block::SparseMatrixCSC{Tv,Ti},ibr::𝕫,ibc::𝕫,factor::ℝ=1.) where{Tv,Ti<:Integer}
    gv              = out.nzval
    asm.pigr[ibc][ibr,1] > 0 || muscadeerror("Trying to addin! into an empty block")
    for ilc         = 1:size(block,2)
        igv         = asm.pigr[ibc][ibr,ilc]
        pilr,lv     = block.colptr,block.nzval 
        for ilv     = pilr[ilc]:pilr[ilc+1]-1 
            gv[igv]+= lv[ilv]*factor 
            igv    += 1
        end
    end
end
"""
    addin!(bigvec,block::Vector,asm,ibc)

Add a vector into one of the blocks of a large full vector.  Use [`prepare`](@ref) to create `asm`.

See also: [`prepare`](@ref)
""" 
function addin!(asm::BlockSparseAssembler,out::Vector{Tv},block::Vector{Tv},ibc::𝕫) where{Tv}
    for ilc         = 1:length(block)
        out[asm.pgc[ibc]-1+ilc] += block[ilc]
    end
end
# disassemble a block from a big-vector
disblock(asm::BlockSparseAssembler,out::Vector,ibc::𝕫) = view(out,asm.pgc[ibc]:(asm.pgc[ibc+1]-1))



