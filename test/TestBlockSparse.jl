module TestBlockSparse

using Test,SparseArrays
using Muscade

using LinearAlgebra,SparseArrays

nrow = 3
ncol = 2
B = Matrix{SparseMatrixCSC{𝕣,𝕫}}(undef,nrow,ncol)
for irow = 1:nrow
    for icol = 1:ncol
        if irow>1 || icol>1
            B[irow,icol] = sparse([1,1,2,3,3],[1,2,2,2,3],randn(5))
        end
    end
end
m,blkasm = prepare(B)

zero!(m)
for irow = 1:nrow
    for icol = 1:ncol
        if irow>1 || icol>1
            addin!(m,B[irow,icol],blkasm,irow,icol)
        end
    end
end
m2 = Matrix(m)

@testset "BlockSparse" begin
    @test m2[4:6,1:3] == B[2,1]
end
end