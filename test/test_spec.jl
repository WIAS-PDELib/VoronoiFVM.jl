module test_spec
using VoronoiFVM
using SparseArrays
using Test

function test_addrows(m, n, mx)
    sa = sprand(m, n, 0.5)
    da = Matrix(sa)
    SparseMatrixCSC(VoronoiFVM.addzrows(da, mx)) == VoronoiFVM.addzrows(sa, mx)
    return sa == da
end

function runtests()
    @test test_addrows(0, 100, 1)
    @test test_addrows(1, 100, 1)
    @test test_addrows(1, 100, 2)
    @test test_addrows(1, 100, 4)
    @test test_addrows(10, 100, 11)
    @test test_addrows(10, 100, 1)
    return true
end

end
