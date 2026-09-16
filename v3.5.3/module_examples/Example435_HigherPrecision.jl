#=

# 434: Higher precision
 ([source code](@__SOURCE_URL__))
 
=#
module Example435_HigherPrecision

using VoronoiFVM, ExtendableGrids
using GridVisualize
using ExtendableSparse
using LinearSolve


## Flux function which describes the flux
## between neighboring control volumes
function g!(f, u, edge, data)
    f[1] = u[1, 1]^2 - u[1, 2]^2
    return nothing
end

function main(; Plotter = nothing, n = 5, assembly = :edgewisem, valuetype = Float64)
    nspecies = 1
    ispec = 1
    X = collect(0:(1.0 / n):1)
    grid = simplexgrid(X, X)
    physics = VoronoiFVM.Physics(; flux = g!)
    sys = VoronoiFVM.System(grid, physics; valuetype, assembly = assembly)
    enable_species!(sys, ispec, [1])
    boundary_dirichlet!(sys, ispec, 1, 0.1)
    boundary_dirichlet!(sys, ispec, 3, 1.0)
    return solution = solve(sys; inival = 0.1, verbose = "n")
    #    vis = GridVisualizer(; Plotter = Plotter)
    #    scalarplot!(vis, grid, solution[1, :]; clear = true, colormap = :summer)
    #    reveal(vis)
end

end
