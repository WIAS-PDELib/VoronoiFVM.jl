# # 220: 2D Nonlinear Poisson with boundary reaction and boundary species
# ([source code](@__SOURCE_URL__))

module Example220_NonlinearPoisson2D_BoundarySpecies

using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize

function main(; n = 10, Plotter = nothing, verbose = false, unknown_storage = :sparse)
    h = 1.0 / convert(Float64, n)
    X = collect(0.0:h:1.0)
    Y = collect(0.0:h:1.0)

    grid = simplexgrid(X, Y)

    k = 1.0
    eps::Float64 = 1.0
    physics = VoronoiFVM.Physics(;
        breaction = function (f, u, node, data)
            if node.region == 2
                f[1] = k * (u[1] - u[3])
                f[3] = k * (u[3] - u[1]) + k * (u[3] - u[2])
                f[2] = k * (u[2] - u[3])
            end
            return nothing
        end, bstorage = function (f, u, node, data)
            if node.region == 2
                f[3] = u[3]
            end
            return nothing
        end, flux = function (f, u, edge, data)
            f[1] = eps * (u[1, 1] - u[1, 2])
            f[2] = eps * (u[2, 1] - u[2, 2])
            return nothing
        end, source = function (f, node, data)
            x1 = node[1] - 0.5
            x2 = node[2] - 0.5
            f[1] = exp(-20.0 * (x1^2 + x2^2))
            return nothing
        end, storage = function (f, u, node, data)
            f[1] = u[1]
            f[2] = u[2]
            return nothing
        end
    )

    sys = VoronoiFVM.System(grid, physics; unknown_storage = unknown_storage)

    enable_species!(sys, 1, [1])
    enable_species!(sys, 2, [1])
    enable_boundary_species!(sys, 3, [2])

    function tran32!(a, b)
        return a[1] = b[2]
    end

    bgrid2 = subgrid(grid, [2]; boundary = true, transform = tran32!)

    inival = unknowns(sys)
    inival .= 0.0

    eps = 1.0e-2

    control = VoronoiFVM.NewtonControl()
    control.verbose = verbose
    control.reltol_linear = 1.0e-5
    control.reltol = 1.0e-5
    tstep = 0.01
    time = 0.0
    istep = 0
    u5 = 0
    p = GridVisualizer(; Plotter = Plotter, layout = (3, 1))
    while time < 1
        time = time + tstep
        U = solve(sys; inival, control, tstep)
        inival .= U
        if verbose
            @printf("time=%g\n", time)
        end
        tstep *= 1.0
        istep = istep + 1
        U_bound = view(U[3, :], bgrid2)
        u5 = U_bound[5]
        scalarplot!(p[1, 1], grid, U[1, :]; clear = true)
        scalarplot!(p[2, 1], grid, U[2, :])
        scalarplot!(p[3, 1], bgrid2, U_bound; show = true, flimits = (0, 0.0025))
    end
    return u5
end

using Test
function runtests()
    @test main(; unknown_storage = :sparse) ≈ 0.0020781361856598
    main(; unknown_storage = :dense) ≈ 0.0020781361856598
    return nothing
end
end
