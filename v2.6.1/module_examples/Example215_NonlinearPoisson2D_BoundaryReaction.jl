# # 215: 2D Nonlinear Poisson with boundary reaction
# ([source code](@__SOURCE_URL__))

module Example215_NonlinearPoisson2D_BoundaryReaction

using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize
using ExtendableSparse

function main(;
        n = 10, Plotter = nothing, verbose = false, unknown_storage = :sparse, assembly = :edgewise,
        tend = 100
    )
    h = 1.0 / convert(Float64, n)
    X = collect(0.0:h:1.0)
    Y = collect(0.0:h:1.0)

    grid = simplexgrid(X, Y)

    eps = 1.0e-2
    physics = VoronoiFVM.Physics(;
        breaction = function (f, u, node, data)
            if node.region == 2
                f[1] = 1 * (u[1] - u[2])
                f[2] = 1 * (u[2] - u[1])
            else
                f[1] = 0
                f[2] = 0
            end
            return nothing
        end, flux = function (f, u, edge, data)
            f[1] = eps * (u[1, 1] - u[1, 2])
            f[2] = eps * (u[2, 1] - u[2, 2])
            return nothing
        end, storage = function (f, u, node, data)
            f[1] = u[1]
            f[2] = u[2]
            return nothing
        end
    )

    sys = VoronoiFVM.System(grid, physics; unknown_storage = unknown_storage, assembly = assembly)
    enable_species!(sys, 1, [1])
    enable_species!(sys, 2, [1])

    inival = unknowns(sys)
    inival[1, :] .= map((x, y) -> exp(-5.0 * ((x - 0.5)^2 + (y - 0.5)^2)), grid)
    inival[2, :] .= 0

    control = VoronoiFVM.NewtonControl()
    control.verbose = verbose
    control.reltol_linear = 1.0e-5

    tstep = 0.01
    time = 0.0
    istep = 0
    u25 = 0

    p = GridVisualizer(; Plotter = Plotter, layout = (2, 1))
    while time < tend
        time = time + tstep
        U = solve(sys; inival, control, tstep)
        inival .= U
        if verbose
            @printf("time=%g\n", time)
        end
        I = integrate(sys, physics.storage, U)
        Uall = sum(I)
        tstep *= 1.2
        istep = istep + 1
        u25 = U[25]
        scalarplot!(
            p[1, 1], grid, U[1, :];
            title = @sprintf("U1: %.3g U1+U2:%8.3g", I[1, 1], Uall),
            flimits = (0, 1)
        )
        scalarplot!(
            p[2, 1], grid, U[2, :]; title = @sprintf("U2: %.3g", I[2, 1]),
            flimits = (0, 1)
        )
        reveal(p)
    end
    return u25
end

using Test
function runtests()
    testval = 0.2760603343272377
    @test main(; unknown_storage = :dense, assembly = :edgewise) ≈ testval &&
        main(; unknown_storage = :sparse, assembly = :cellwise) ≈ testval &&
        main(; unknown_storage = :dense, assembly = :cellwise) ≈ testval
    return nothing
end
end
