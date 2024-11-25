#=

# 201: 2D Laplace equation 
([source code](@__SOURCE_URL__))

=#

module Example201_Laplace2D

using VoronoiFVM, ExtendableGrids
using GridVisualize
using LinearAlgebra
import Metis

## Flux function which describes the flux
## between neighboring control volumes
function g!(f, u, edge, data)
    f[1] = u[1, 1] - u[1, 2]
end

function main(; Plotter = nothing, n = 5, is_linear = true, assembly = :edgewise)
    nspecies = 1
    ispec = 1
    X = collect(0:(1.0 / n):1)
    grid = simplexgrid(X, X)
    grid=partition(grid, PlainMetisPartitioning(npart=20))
    @show grid
    physics = VoronoiFVM.Physics(; flux = g!)
    sys = VoronoiFVM.System(grid, physics; is_linear = is_linear, assembly = assembly)
    enable_species!(sys, ispec, [1])
    boundary_dirichlet!(sys, ispec, 1, 0.0)
    boundary_dirichlet!(sys, ispec, 3, 1.0)
    solution = solve(sys; inival = 0)
    nf = nodeflux(sys, solution)
    vis = GridVisualizer(; Plotter = Plotter)
    scalarplot!(vis, grid, solution[1, :]; clear = true, colormap = :summer)
    vectorplot!(vis, grid, nf[:, 1, :]; clear = false, vscale = 0.5, rasterpoints=10)
    reveal(vis)
    return norm(solution) + norm(nf)
end

## Called by unit test

using Test
function runtests()
    testval = 9.63318042491699

    @test main(; assembly = :edgewise) ≈ testval &&
          main(; assembly = :cellwise) ≈ testval
end

end
