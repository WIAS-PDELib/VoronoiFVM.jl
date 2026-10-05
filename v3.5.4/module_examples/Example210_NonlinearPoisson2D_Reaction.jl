# # 210: 2D Nonlinear Poisson with reaction
# ([source code](@__SOURCE_URL__))

module Example210_NonlinearPoisson2D_Reaction

using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize
using LinearSolve, ExtendableSparse
import Metis

function main(;
        n = 10, Plotter = nothing, verbose = false,
        unknown_storage = :sparse, assembly = :edgewise, tstep = 0.01
    )

    X = range(0, 1, length = n + 1)
    Y = range(0, 1, length = n + 1)

    grid = simplexgrid(X, Y)

    #    grid = partition(grid, PlainMetisPartitioning(npart = 10))
    @show grid
    data = (eps = 1.0e-2, k = 1.0)

    function reaction!(f, u, node, data)
        f[1] = data.k * (u[1] - u[2])
        f[2] = data.k * (u[2] - u[1])
        return nothing

    end

    function flux!(f, u, edge, data)
        f[1] = data.eps * (u[1, 1] - u[1, 2])
        f[2] = data.eps * (u[2, 1] - u[2, 2])
        return nothing
    end

    function source!(f, node, data)
        x1 = node[1] - 0.5
        x2 = node[2] - 0.5
        f[1] = exp(-20 * (x1^2 + x2^2))
        return nothing
    end

    function storage!(f, u, node, data)
        f[1] = u[1]
        f[2] = u[2]
        return nothing
    end

    sys = VoronoiFVM.System(
        grid;
        data,
        flux = flux!,
        storage = storage!,
        reaction = reaction!,
        source = source!,
        unknown_storage = unknown_storage, assembly = assembly
    )

    enable_species!(sys, 1, [1])
    enable_species!(sys, 2, [1])

    inival = unknowns(sys)
    inival .= 0.0

    control = VoronoiFVM.SolverControl(;
        verbose,
        Δt = tstep,
        Δt_min = tstep,
        Δt_max = tstep,
        Δu_opt = 1.0e5,
        method_linear = KrylovJL_BICGSTAB(precs = LinearSolvePreconBuilder(UMFPACKFactorization())),
        factorize_every_timestep = 2
    )

    p = GridVisualizer(; Plotter = Plotter, layout = (2, 1))
    function post(U, Uold, t, Δt)
        scalarplot!(p[1, 1], grid, U[1, :]; clear = true, limits = (0, 0.5))
        scalarplot!(p[2, 1], grid, U[2, :]; show = true, limits = (0, 0.5))
        return
    end

    tsol = solve(sys; inival, times = (0, 1), control, post)
    testval = sum(tsol.u[end])
    return testval
end

using Test
function runtests()
    testval = 16.01812472041518
    @test main(; unknown_storage = :sparse, assembly = :edgewise) ≈ testval &&
        main(; unknown_storage = :dense, assembly = :edgewise) ≈ testval &&
        main(; unknown_storage = :sparse, assembly = :cellwise) ≈ testval &&
        main(; unknown_storage = :dense, assembly = :cellwise) ≈ testval
    return nothing
end

end
