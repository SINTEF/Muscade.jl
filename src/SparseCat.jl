# i index
# n length
# p ptr (using the colptr compression convention)
#
# b block  - indexing blocs
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

"""
    bigsparse = cat(blocks::Matrix{<:SparseMatrixCSC})

Assemble sparse blocs into a large sparse matrix.

See also: [`cat!`](@ref)
""" 
function Base.cat(B::Matrix{SparseMatrixCSC{Tv,Ti}}) where{Tv,Ti<:Integer}
    nbr,nbc                 = size(B)  
    nbr>0 && nbc>0 || muscadeerror("must have length(B)>0")

    pgr                     = 𝕫1(undef,nbr+1)         # pointers into solution vector, with global*solution=rhs
    pgr[1]                  = 1
    for ibr                 = 1:nbr
        pgr[ibr+1]          = pgr[ibr]+size(B[ibr,1],1)
    end 
    ngr                     = pgr[end]-1

    pgc                     = 𝕫1(undef,nbc+1)         # pointers into rhs vector
    pgc[1]                  = 1
    for ibc                 = 1:nbc
        pgc[ibc+1]          = pgc[ibc]+size(B[1,ibc],2)
    end 
    ngc                     = pgc[end]-1

    ngv                     = 0 
    for ibc                 = 1:nbc
        for ibr             = 1:nbr
            ngv            += nnz(B[ibr,ibc])
        end
    end

    igr                     = 𝕫1(undef,ngv) 
    gv                      = 𝕣1(undef,ngv)
    pigr                    = 𝕫1(undef,ngc+1) 
    pigr[1]                 = 1

    igv                     = 1
    for ibc                 = 1:nbc
        for ilc             = 1:size(B[1,ibc],2)
            igc             = ilc + pgc[ibc]-1
            for ibr         = 1:nbr
                pilr,ilr,lv = B[ibr,ibc].colptr, B[ibr,ibc].rowval,B[ibr,ibc].nzval 
                for ilv     = pilr[ilc]:pilr[ilc+1]-1 
                    igr[igv]= pgr[ibr]-1 + ilr[ilv]  
                    gv[igv] = lv[ilv]
                    igv    += 1    
                end
                pigr[igc+1] = igv   
            end
        end
    end
    return SparseMatrixCSC(ngr,ngc,pigr,igr,gv)    
end
"""
    bigsparse = cat!(bigsparse,blocs::Matrix{<:SparseMatrixCSC})

Assemble sparse blocs into a large sparse matrix.  Will fail silently or throw an error unless
`bigsparse` has the correct sparsity structure for the given `blocs`. 

See also: [`cat`](@ref)
""" 
function cat!(out::SparseMatrixCSC{Tv,Ti},B::Matrix{SparseMatrixCSC{Tv,Ti}}) where{Tv,Ti<:Integer}
    nbr,nbc                 = size(B)  
    igv,gv                  = 1,out.nzval
    for ibc                 = 1:nbc
        for ilc             = 1:size(B[1,ibc],2)
            for ibr         = 1:nbr
                pilr,lv     = B[ibr,ibc].colptr,B[ibr,ibc].nzval 
                for ilv     = pilr[ilc]:pilr[ilc+1]-1 
                    gv[igv] = lv[ilv]
                    igv    += 1    
                end
            end
        end
    end
end
