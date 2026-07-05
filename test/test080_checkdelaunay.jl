module test_checkdelaunay
using Test
using ExtendableGrids: Coordinates, simplexgrid
using VoronoiFVM: nondelaunay

function runtests()
    X = 0:0.1:10
    g = simplexgrid(X, X)
    @test length(nondelaunay(g)) == 0


    coord = g[Coordinates]
    for i in 1:size(coord, 2)
        coord[1, i] += 0.01 * (rand() - 0.5)
        coord[2, i] += 0.01 * (rand() - 0.5)
    end

    @test length(nondelaunay(g)) > 0
    return nothing
end

end
