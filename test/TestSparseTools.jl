module TestSparseTools

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
bigsparse,bigsparseasm,bigveca   = prepare(pattern)


zero!(bigsparse)
for irow = 1:nrow, icol = 1:ncol
    if irow>1 || icol>1
        addin!(bigsparseasm,bigsparse,block,irow,icol)
    end
end

big2 = Matrix(bigsparse)
@testset "BlockSparseFromSparse" begin
    @test big2[4:6,1:3] == block
    @test big2[4:6,4:6] == block
end

# Sparse pattern of blocks #2

nrow = 3
ncol = 3
block        = sparse([1,1,1,2,3,4],[1,2,4,2,3,4],ones(6))  # block
pattern      = sparse([1,2,3],[1,2,3],[block,block,block]) # pattern of the blocks in bigsparse
bigsparse,bigsparseasm,bigvecasm,bigvecdis   = prepare(pattern)

@testset "bigsparseasm 2" begin
    @test bigsparseasm.colptr == [1,2,3,4]
    @test bigsparseasm.rowval == [1,2,3]
    @test bigsparseasm.nzval == [[1, 2, 3, 4, 5, 6],[7, 8, 9, 10, 11, 12],[13, 14, 15, 16, 17, 18]]
    @test bigvecasm          == [1,5,9,13]
    @test bigvecdis          == [1,5,9,13]
end


zero!(bigsparse)
(i,j,v) = findnz(pattern)
for k = 1:nnz(pattern) 
    addin!(bigsparseasm,bigsparse,block,i[k],j[k])
end

big2 = Matrix(bigsparse)
@testset "BlockSparseFromSparse 2" begin
    @test big2[1:4,1:4] == block
    @test big2[5:8,5:8] == block
    @test big2[9:12,9:12] == block
end

i = [3, 7, 2, 3, 6, 2, 7, 9, 2, 6, 4, 5, 9, 3, 9, 1, 7, 8, 10, 4, 9, 7]
j = [2, 2, 3, 3, 3, 4, 4, 4, 6, 6, 7, 7, 7, 8, 8, 9, 9, 9, 9, 10, 10, 11]
v = [0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]

s = sparse(i,j,v)
Muscade.sparser!(s,i->s.nzval[i]>0.5)

@testset "sparser!" begin
    @test s.colptr == [1, 1, 2, 3, 5, 5, 5, 6, 6, 7, 8, 8]
    @test s.rowval == [7, 6, 2, 7, 4, 1, 9]
    @test s.nzval ≈ [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    @test size(s) == (10,11)
end


s = sprandn(10,11,0.2)
nnzs = 0
for (j,nz) ∈ enumerate(s.nzval)
    if nz>0
        s.nzval[j] = 1. 
        global nnzs +=1
    else
        s.nzval[j] = 0.
    end
end
Muscade.sparser!(s,i->s.nzval[i]>0.5)

@testset "sparser!2" begin
    @test nnzs == nnz(s)
end


end