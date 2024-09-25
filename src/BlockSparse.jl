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
function prepare(pattern::SparseMatrixCSC{SparseMatrixCSC{Tv,𝕫},𝕫}) where{Tv} 
    nbr,nbc                       = size(pattern)
    nbr>0 && nbc>0 || muscadeerror("must have length(pattern)>0")
    nlr                           = [-1 for ibr=1:nbr+1]
    nlc                           = [-1 for ibc=1:nbc+1]
    ngv                           = 0
    for ibc                       = 1:nbc
        for pbr                   = pattern.colptr[ibc]:pattern.colptr[ibc+1]-1
            ibr                   = pattern.rowval[pbr]
            block                 = pattern.nzval[pbr]
            nlr[ibr+1]            = block.m
            nlc[ibc+1]            = block.n
            ngv                  += length(block.nzval)
        end
    end
    nlr[1]                        = 1
    nlc[1]                        = 1
    all(nlr.>0) || muscadeerror("every row of the pattern must contain at least one non-zero block")
    all(nlc.>0) || muscadeerror("every column of the pattern contain at least one non-zero block")
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
        asm[ibc]                  = zeros(𝕫,nbr,nlc[ibc+1])    # TODO CPU 
        for ilc                   = 1:nlc[ibc+1]           # for each local column
            igc                   = pgc[ibc]-1 + ilc 
            for pbr               = pattern.colptr[ibc]:pattern.colptr[ibc+1]-1 # for each block row
                ibr               = pattern.rowval[pbr]
                block             = pattern.nzval[pbr]
                pilr,ilr          = block.colptr, block.rowval
                asm[ibc][ibr,ilc] = igv
                for ilv           = pilr[ilc]:pilr[ilc+1]-1 
                    igr[igv]      = pgr[ibr]-1 + ilr[ilv]  
                    igv          += 1    
                end
            end
            pigr[igc+1]           = igv   
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
    asm.pigr[ibc][ibr,1] > 0 || muscadeerror(@printf("Trying to addin! into an empty block: [%i,%i]\n",ibr,ibc))
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
function addin!(asm::BlockSparseAssembler,out::Vector{Tv},block::Vector{Tv},ibc::𝕫,factor::ℝ=1.) where{Tv}
    for ilc         = 1:length(block)
        out[asm.pgc[ibc]-1+ilc] += block[ilc]*factor
    end
end
# disassemble a block from a big-vector
disblock(asm::BlockSparseAssembler,v::Vector,ibc::𝕫) = view(v,asm.pgc[ibc]:(asm.pgc[ibc+1]-1))

