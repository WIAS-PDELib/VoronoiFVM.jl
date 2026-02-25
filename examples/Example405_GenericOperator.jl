#=

# 405: Generic operator
 ([source code](@__SOURCE_URL__))

Handle an operator which does not fit into the storage/flux/reaction API.
This uses automatic sparsity detection.
=#

module Example405_GenericOperator
using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize

function main(; n = 10, Plotter = nothing, verbose = false, unknown_storage = :sparse)
    ## Same as Example102 with upwind

    ## Create a one-dimensional discretization
    h = 1.0 / convert(Float64, n)
    X = collect(0:h:1)
    grid = simplexgrid(X)

    ## A parameter which is "passed" to the flux function via scope
    D = 1.0e-2
    v = 1.0

    ## This generic operator works on the full solution seen as linear vector, and indexing
    ## shall be done  by reshaping it into a solution vector of the system.
    ## Here, instead of the flux function we provide a "generic operator"
    ## which provides the stiffness part of the problem. Its sparsity is detected automatically
    ## using Symbolics.jl
    function generic_operator!(f0, u0, sys, data)
        f = reshape(f0, sys)
        u = reshape(u0, sys)
        for i in 1:(length(X) - 1)
            du = D * (u[1, i] - u[1, i + 1]) / (X[i + 1] - X[i]) +
                v * (v > 0 ? u[1, i] : u[1, i + 1])
            f[1, i] += du
            f[1, i + 1] -= du
        end
        return nothing
    end

    ## Create a physics structure
    physics = VoronoiFVM.Physics(; generic = generic_operator!)

    ## Create a finite volume system - either
    ## in the dense or  the sparse version.
    ## The difference is in the way the solution object
    ## is stored - as dense or as sparse matrix

    sys = VoronoiFVM.System(grid, physics; unknown_storage = unknown_storage)

    ## Add species 1 to region 1
    enable_species!(sys, 1, [1])

    ## Set boundary conditions
    boundary_dirichlet!(sys, 1, 1, 0.0)
    boundary_dirichlet!(sys, 1, 2, 1.0)

    ## Create a solution array
    inival = unknowns(sys; inival = 0.5)
    solution = unknowns(sys)


    ## Stationary solution of the problem
    solution = solve(sys; inival = 0.5, verbose)

    scalarplot(
        grid, solution[1, :]; title = "Nonlinear Poisson", Plotter = Plotter,
        resolution = (300, 300)
    )
    return sum(solution)
end

using Test
function runtests()
    testval = 1.099999999614456
    @test main(; unknown_storage = :sparse) ≈ testval &&
        main(; unknown_storage = :dense) ≈ testval
    return nothing
end

end
