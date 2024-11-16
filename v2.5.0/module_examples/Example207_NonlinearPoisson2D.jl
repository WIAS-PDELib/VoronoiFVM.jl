# # 207: 2D Nonlinear Poisson equation
# ([source code](@__SOURCE_URL__))

module Example207_NonlinearPoisson2D

using Printf
using VoronoiFVM
using ExtendableGrids
using ExtendableSparse: ILUZeroPreconBuilder
using GridVisualize
using LinearSolve
using ILUZero

function main(; n = 10, Plotter = nothing, verbose = false, unknown_storage = :sparse,
              method_linear = nothing, assembly = :edgewise)
    h = 1.0 / convert(Float64, n)
    X = collect(0.0:h:1.0)
    Y = collect(0.0:h:1.0)

    grid = simplexgrid(X, Y)

    eps = 1.0e-2

    physics = VoronoiFVM.Physics(; reaction = function (f, u, node, data)
                                     f[1] = u[1]^2
                                 end, flux = function (f, u, edge, data)
                                     f[1] = eps * (u[1, 1]^2 - u[1, 2]^2)
                                 end, source = function (f, node, data)
                                     x1 = node[1] - 0.5
                                     x2 = node[2] - 0.5
                                     f[1] = exp(-20.0 * (x1^2 + x2^2))
                                 end, storage = function (f, u, node, data)
                                     f[1] = u[1]
                                 end)
    sys = VoronoiFVM.System(grid, physics; unknown_storage, assembly = assembly)
    enable_species!(sys, 1, [1])

    boundary_dirichlet!(sys, 1, 2, 0.1)
    boundary_dirichlet!(sys, 1, 4, 0.1)

    inival = unknowns(sys)
    inival .= 0.5

    control = VoronoiFVM.NewtonControl()
    control.verbose = verbose
    control.reltol_linear = 1.0e-5
    control.method_linear = method_linear
    tstep = 0.01
    time = 0.0
    u15 = 0
    p = GridVisualizer(; Plotter = Plotter)
    while time < 1.0
        time = time + tstep
        U = solve(sys; inival, control, tstep)
        u15 = U[15]
        inival .= U

        scalarplot!(p[1, 1], grid, U[1, :]; Plotter = Plotter, clear = true, show = true)
        tstep *= 1.0
    end
    return u15
end

using Test
function runtests()
    # test at once for iterative solution here
    testval = 0.3554284760906605
    @test main(; unknown_storage = :sparse, assembly = :edgewise) ≈ testval &&
          main(; unknown_storage = :dense, assembly = :edgewise) ≈ testval &&
          main(; unknown_storage = :sparse, method_linear = KrylovJL_CG(precs=ILUZeroPreconBuilder()),
               assembly = :edgewise) ≈ testval &&
          main(; unknown_storage = :dense, method_linear =  KrylovJL_CG(precs=ILUZeroPreconBuilder()),
               assembly = :edgewise) ≈ testval &&
          main(; unknown_storage = :sparse, assembly = :cellwise) ≈ testval &&
          main(; unknown_storage = :dense, assembly = :cellwise) ≈ testval &&
          main(; unknown_storage = :sparse, method_linear =  KrylovJL_CG(precs=ILUZeroPreconBuilder()),
               assembly = :cellwise) ≈ testval &&
          main(; unknown_storage = :dense, method_linear =  KrylovJL_CG(precs=ILUZeroPreconBuilder()),
               assembly = :cellwise) ≈ testval
end
end
