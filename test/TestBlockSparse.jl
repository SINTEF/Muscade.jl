#module TestBlockSparse

# cd("C:\\Users\\philippem\\.julia\\dev\\Muscade")
# using Pkg 
# Pkg.activate(".")

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
        addin!(bigasm,big,block,irow,icol)
    end
end

big2 = Matrix(big)
@testset "BlockSparseFromSparse" begin
    @test big2[4:6,1:3] == block
    @test big2[4:6,4:6] == block
end

# Sparse pattern of blocks #2

nrow = 3
ncol = 3
block        = sparse([1,1,1,2,3,4],[1,2,4,2,3,4],ones(6))  # block
pattern      = sparse([1,2,3],[1,2,3],[block,block,block]) # pattern of the blocks in bigsparse
big,bigasm   = prepare(pattern)

@testset "bigasm 2" begin
    @test bigasm.pibr == [1,2,3,4]
    @test bigasm.ibr == [1,2,3]
    @test bigasm.igv == [[1, 2, 3, 4, 5, 6],[7, 8, 9, 10, 11, 12],[13, 14, 15, 16, 17, 18]]
    @test bigasm.pgr == [1,5,9,13]
    @test bigasm.pgc == [1,5,9,13]
end


zero!(big)
(i,j,v) = findnz(pattern)
for k = 1:nnz(pattern) 
    addin!(bigasm,big,block,i[k],j[k])
end

big2 = Matrix(big)
@testset "BlockSparseFromSparse 2" begin
    @test big2[1:4,1:4] == block
    @test big2[5:8,5:8] == block
    @test big2[9:12,9:12] == block
end

#end