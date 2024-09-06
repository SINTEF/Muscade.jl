module TestBlockSparse

using Test,SparseArrays
using Muscade

# Sparse pattern of blocks

nrow = 3
ncol = 2
block        = sparse([1,1,2,3,3],[1,2,2,2,3],randn(5))  # block
pattern      = sparse([1,2,2,3,3],[2,1,2,1,2],[block,block,block,block,block]) # pattern of the blocks in bigsparse
big,bigasm   = prepare(pattern)

zero!(big)
for irow = 1:nrow, icol = 1:ncol
    if irow>1 || icol>1
        addin!(big,block,bigasm,irow,icol)
    end
end

big2 = Matrix(big)
@testset "BlockSparseFromSparse" begin
    @test big2[4:6,1:3] == block
    @test big2[4:6,4:6] == block
end

# Full pattern of blocks

pattern = Matrix{SparseMatrixCSC{𝕣,𝕫}}(undef,nrow,ncol)
for irow = 1:nrow,icol = 1:ncol
    if irow>1 || icol>1
        pattern[irow,icol] = sparse([1,1,2,3,3],[1,2,2,2,3],randn(5))
    end
end
big,bigasm = prepare(pattern)  # uses only the structure of pattern, not the values

zero!(big)
for irow = 1:nrow, icol = 1:ncol
    if irow>1 || icol>1
        block = pattern[irow,icol] # for convenience, in the test, we use the values of pattern as blocks, but that's not typical usage
        addin!(big,block,bigasm,irow,icol) 
    end
end

big2 = Matrix(big)
@testset "BlockSparseFromMatrix" begin
    @test big2[4:6,1:3] == pattern[2,1]
    @test big2[4:6,4:6] == pattern[2,2]
end

end