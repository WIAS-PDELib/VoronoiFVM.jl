module test_formfactors
using Test
using ExtendableGrids
using VoronoiFVM: cellfactors!, bfacefactors!

randpoint = rand(-10:0.01:10, 2)
function ttri(; ntest = 100)
    cellnodes = [1 2 3;]'
    icell = 1
    epar2d = zeros(3)
    npar2d = zeros(3)
    epar3d = zeros(3)
    npar3d = zeros(3)
    for i in 1:100
        coord2d = rand(-10:0.01:10, 2, 3)
        coord3d = vcat(coord2d, [0.0, 0, 0]')

        cellfactors!(Triangle2D, Cartesian2D, coord2d, cellnodes, 1, npar2d, epar2d)
        bfacefactors!(Triangle2D, Cartesian3D, coord3d, cellnodes, 1, npar3d, epar3d)

        @test npar3d ≈ npar2d
        @test epar3d ≈ epar2d
    end
    return
end


end
