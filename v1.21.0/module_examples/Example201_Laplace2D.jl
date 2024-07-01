#=

# 201: 2D Laplace equation 
([source code](@__SOURCE_URL__))

=#

module Example201_Laplace2D

using VoronoiFVM, ExtendableGrids
using GridVisualize
using LinearAlgebra
## Flux function which describes the flux
## between neighboring control volumes
function g!(f, u, edge)
    f[1] = u[1, 1] - u[1, 2]
end

function main(; Plotter = nothing, n = 5, is_linear = true, assembly = :edgewise)
    nspecies = 1
    ispec = 1
    X = collect(0:(1.0 / n):1)
    grid = simplexgrid(X, X)

    physics = VoronoiFVM.Physics(; flux = g!)
    sys = VoronoiFVM.System(grid, physics; is_linear = is_linear, assembly = assembly)
    enable_species!(sys, ispec, [1])
    boundary_dirichlet!(sys, ispec, 1, 0.0)
    boundary_dirichlet!(sys, ispec, 3, 1.0)
    solution = solve(sys; inival = 0)
    nf = nodeflux(sys, solution)
    vis = GridVisualizer(; Plotter = Plotter)
    scalarplot!(vis, grid, solution[1, :]; clear = true, colormap = :summer)
    vectorplot!(vis, grid, nf[:, 1, :]; clear = false, spacing = 0.1, vscale = 0.5)
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
