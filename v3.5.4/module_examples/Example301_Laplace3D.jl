#=

# 301: 3D Laplace equation 
([source code](@__SOURCE_URL__))

=#

module Example301_Laplace3D

using VoronoiFVM, ExtendableGrids
using GridVisualize

## Flux function which describes the flux
## between neighboring control volumes
function g!(f, u, edge, data)
    f[1] = u[1, 1] - u[1, 2]
    return nothing
end

function s(f, node, data)
    n = view(node.coord, :, node.index)
    f[1] = n[1] * sin(5.0 * n[2]) * exp(n[3])
    return nothing
end

function main(; Plotter = nothing, n = 5, assembly = :edgewise)
    nspecies = 1
    ispec = 1
    X = collect(0:(1 / n):1)
    grid = simplexgrid(X, X, X)
    physics = VoronoiFVM.Physics(; flux = g!, source = s)
    sys = VoronoiFVM.System(grid, physics; assembly = assembly)
    enable_species!(sys, ispec, [1])
    boundary_dirichlet!(sys, ispec, 5, 0.0)
    boundary_dirichlet!(sys, ispec, 6, 0.0)
    solution = solve(sys)
    scalarplot(grid, solution[1, :]; Plotter = Plotter)
    return solution[43]
end

## Called by unit test

using Test
function runtests()
    testval = 0.012234524449380824
    @test main(; assembly = :edgewise) ≈ testval &&
        main(; assembly = :cellwise) ≈ testval
    return nothing
end

end
