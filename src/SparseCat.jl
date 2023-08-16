# i index
# n length
# p ptr (using the colptr compression convention)
#
# b block  - indexing blocks
# l local  - indexing inside a block
# g global - indexing the whole matrix
#
# r row
# c col
# v non zero value
#
# sparse.rowval → ilr
# sparse.colptr → pilr  
# ilv    = pilr[ilc]+ i-1
# irow   = ilr[ilv]  

struct BlockSparseAssembler
    pigr::Vector{𝕫2}
    pgc ::𝕫1
end

"""
    bigsparse,asm = blocksparse(blocks::Matrix{<:SparseMatrixCSC})

Prepare for the assembly of sparse blocks into a large sparse matrix. Unassigned blocks (undef'ed) will be treated as all-zero.
`bigsparse` is allocated, with the correct sparsity structure, but its `nzval` undef'ed.    
Where some blocks share the same sparsity structure, `blocks` can have `===` elements.

See also: [`zero!`](@ref),[`addin!`](@ref),[`cat!`](@ref)
""" 
function blocksparse(B::Matrix{SparseMatrixCSC{Tv,Ti}}) where{Tv,Ti<:Integer}
    nbr,nbc                 = size(B)  
    nbr>0 && nbc>0 || muscadeerror("must have length(B)>0")

    pgr                     = 𝕫1(undef,nbr+1)         # pointers into solution vector, with global*solution=rhs
    pgr[1]                  = 1
    for ibr                 = 1:nbr
        nlr                 = 0
        for ibc             = 1:nbc
            if isassigned(B,ibr,ibc)
                nlr         = size(B[ibr,ibc],1)
                break
            end
            ibc<nbc || muscadeerror("B has an empty row")
        end
        nlr > 0 || muscadeerror("every block-row must contain at least one assigned Sparse")
        pgr[ibr+1]          = pgr[ibr]+nlr
    end 
    ngr                     = pgr[end]-1

    pgc                     = 𝕫1(undef,nbc+1)         # pointers into rhs vector
    pgc[1]                  = 1
    for ibc                 = 1:nbc
        nlc                 = 0
        for ibr             = 1:nbr
            if isassigned(B,ibr,ibc)
                nlc         = size(B[ibr,ibc],2)
                break
            end
            ibc<nbc || muscadeerror("B has an empty column")
        end
        nlc > 0 || muscadeerror("every block-column must contain at least one assigned Sparse")
        pgc[ibc+1]          = pgc[ibc]+nlc
    end 
    ngc                     = pgc[end]-1

    ngv                     = 0 
    for ibc                 = 1:nbc
        for ibr             = 1:nbr
            if isassigned(B,ibr,ibc)
                ngv        += nnz(B[ibr,ibc])
            end
        end
    end

    igr                     = 𝕫1(undef,ngv  ) 
    gv                      = 𝕣1(undef,ngv  )
    pigr                    = 𝕫1(undef,ngc+1) 
    pigr[1]                 = 1

    asm                     = Vector{𝕫2}(undef,nbr)  
    igv                     = 1
    for ibc                 = 1:nbc
       nlc                  = pgc[ibc+1]-pgc[ibc]
       asm[ibc]             = zeros(𝕫,nbr,nlc)
       for ilc              = 1:nlc
            igc             = ilc + pgc[ibc]-1
            for ibr         = 1:nbr
                if isassigned(B,ibr,ibc)
                    asm[ibc][ibr,ilc] = igv
                    pilr,ilr    = B[ibr,ibc].colptr, B[ibr,ibc].rowval
                    for ilv     = pilr[ilc]:pilr[ilc+1]-1 
                        igr[igv]= pgr[ibr]-1 + ilr[ilv]  
                        igv    += 1    
                    end
                # else
                #     pigr[igc] = xx  # REPRISE
                end
                pigr[igc+1] = igv   
            end
        end
    end
    return SparseMatrixCSC(ngr,ngc,pigr,igr,gv),BlockSparseAssembler(asm,pgc)    
end
"""
    cat!(bigsparse,blocks::Matrix{<:SparseMatrixCSC})

Assemble sparse blocks into a preallocated sparse matrix.  Will fail silently or throw an error unless
`bigsparse` has the correct sparsity structure for the given `blocks`. 

See also: [`blocksparse`](@ref),[`zero!`](@ref),[`addin!`](@ref)
""" 
function cat!(out::SparseMatrixCSC{Tv,Ti},B::Matrix{SparseMatrixCSC{Tv,Ti}},asm::BlockSparseAssembler) where{Tv,Ti<:Integer}
    nbr,nbc                 = size(B)  
    igv,gv                  = 1,out.nzval
    for ibc                 = 1:nbc
        nlc                 = asm.pgc[ibc+1]-asm.pgc[ibc]
        for ilc             = 1:nlc
            for ibr         = 1:nbr
                if isassigned(B,ibr,ibc)
                    pilr,lv     = B[ibr,ibc].colptr,B[ibr,ibc].nzval 
                    for ilv     = pilr[ilc]:pilr[ilc+1]-1 
                        gv[igv] = lv[ilv]
                        igv    += 1    
                    end
                end
            end
        end
    end
end
"""
    addin!(bigsparse,block::SparseMatrixCSC,asm,ibr,ibc)

Add a sparse into one of the blocks of a large sparse matrix.  Will fail silently or throw an error unless
`bigsparse` has the correct sparsity structure for the given `blocks`. Use [`blocksparse`](@ref) to
    create `bigsparse` and `asm`.

See also: [`blocksparse`](@ref),[`zero!`](@ref),[`cat!`](@ref)
""" 
function addin!(out::SparseMatrixCSC{Tv,Ti},B::SparseMatrixCSC{Tv,Ti},asm::BlockSparseAssembler,ibr::𝕫,ibc::𝕫) where{Tv,Ti<:Integer}
    gv              = out.nzval
    asm.pigr[ibc][ibr,1] > 0 || muscadeerror("Trying to cat! into an empty block")
    for ilc         = 1:size(B,2)
        igv         = asm.pigr[ibc][ibr,ilc]
        pilr,lv     = B.colptr,B.nzval 
        for ilv     = pilr[ilc]:pilr[ilc+1]-1 
            gv[igv]+= lv[ilv] 
            igv    += 1
        end
    end
end
"""
    addin!(bigvec,block::Vector,asm,ibr,ibc)

Add a vector into one of the blocks of a large full vector.  Use [`blocksparse`](@ref) to create `asm`.

See also: [`blocksparse`](@ref),[`zero!`](@ref),[`cat!`](@ref)
""" 
function addin!(out::Vector{Tv},B::Vector{Tv},asm::BlockSparseAssembler,ibc::𝕫) where{Tv}
    for ilc         = 1:length(B)
        out[asm.pgc[ibc]-1+ilc] += B[ilc]
    end
end
getblock(out::Vector,asm::BlockSparseAssembler,ibc::𝕫) = view(out,asm.pgc[ibc]:(asm.pgc[ibc+1]-1))



