# # 120: Differing species sets in regions, 1D
# ([source code](@__SOURCE_URL__))

module Example120_ThreeRegions1D

using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize
using LinearSolve
using OrdinaryDiffEqRosenbrock
using SciMLBase: NoInit

function reaction(f, u, node, data)
    k = data.k
    if node.region == 1
        f[1] = k[1] * u[1]
        f[2] = -k[1] * u[1]
    elseif node.region == 3
        f[2] = k[3] * u[2]
        f[3] = -k[3] * u[2]
    else
        f[1] = 0
    end
    return nothing
end

function source(f, node, data)
    if node.region == 1
        f[1] = 1.0e-4 * (3.0 - node[1])
    end
    return nothing
end

## Since 0.17.0 one can
## write into the result also where
## the corresponding species has not been enabled
## Species information is used to prevent the assembly.
function correctionflux(f, u, edge, data)
    eps = data.eps
    for i in 1:3
        f[i] = eps[i] * (u[i, 1] - u[i, 2])
    end
    return nothing
end

function correctionstorage(f, u, node, data)
    f .= u
    return nothing
end

## This is the "old" way:
## Write into result only where
## the corresponding species has been enabled
function pickyflux(f, u, edge, data)
    eps = data.eps
    if edge.region == 1
        f[1] = eps[1] * (u[1, 1] - u[1, 2])
        f[2] = eps[2] * (u[2, 1] - u[2, 2])
    elseif edge.region == 2
        f[2] = eps[2] * (u[2, 1] - u[2, 2])
    elseif edge.region == 3
        f[2] = eps[2] * (u[2, 1] - u[2, 2])
        f[3] = eps[3] * (u[3, 1] - u[3, 2])
    end
    return nothing
end

function pickystorage(f, u, node, data)
    if node.region == 1
        f[1] = u[1]
        f[2] = u[2]
    elseif node.region == 2
        f[2] = u[2]
    elseif node.region == 3
        f[2] = u[2]
        f[3] = u[3]
    end
    return nothing
end


function main(;
        n = 30, Plotter = nothing, plot_grid = false, verbose = false,
        unknown_storage = :sparse, tend = 10,
        diffeq = false,
        rely_on_corrections = false, assembly = :edgewise
    )

    X = range(0, 3, length = n)
    grid = simplexgrid(X)
    cellmask!(grid, [0.0], [1.0], 1)
    cellmask!(grid, [1.0], [2.1], 2)
    cellmask!(grid, [1.9], [3.0], 3)

    subgrid1 = subgrid(grid, [1])
    subgrid2 = subgrid(grid, [1, 2, 3])
    subgrid3 = subgrid(grid, [3])

    if plot_grid
        plotgrid(grid; Plotter = Plotter)
        return
    end

    data = (eps = [1, 1, 1], k = [1, 1, 1])

    flux = rely_on_corrections ? correctionflux : pickyflux
    storage = rely_on_corrections ? correctionstorage : pickystorage

    sys = VoronoiFVM.System(
        grid; data,
        flux, reaction, storage, source,
        unknown_storage, assembly
    )

    enable_species!(sys, 1, [1])
    enable_species!(sys, 2, [1, 2, 3])
    enable_species!(sys, 3, [3])

    boundary_dirichlet!(sys, 3, 2, 0.0)

    testval = 0
    p = GridVisualizer(; Plotter = Plotter, layout = (1, 1))

    function plot_timestep(U, time)
        U1 = view(U[1, :], subgrid1)
        U2 = view(U[2, :], subgrid2)
        U3 = view(U[3, :], subgrid3)

        scalarplot!(
            p[1, 1], subgrid1, U1; label = "spec1", color = :darkred,
            xlimits = (0, 3), flimits = (0, 1.0e-3),
            title = @sprintf("three regions t=%.3g", time)
        )
        scalarplot!(
            p[1, 1], subgrid2, U2; label = "spec2", color = :green,
            clear = false
        )
        scalarplot!(
            p[1, 1], subgrid3, U3; label = "spec3", color = :navyblue,
            clear = false, show = true
        )
        return if ismakie(Plotter)
            sleep(0.02)
        end
    end

    if diffeq
        inival = unknowns(sys, inival = 0)
        problem = ODEProblem(sys, inival, (0, tend))
        ## use fixed timesteps just for the purpose of CI
        odesol = solve(problem, Rosenbrock23(), initializealg = NoInit(), dt = 1.0e-2, adaptive = false)
        tsol = reshape(odesol, sys)
    else
        tsol = solve(
            sys; inival = 0, times = (0, tend),
            verbose, Δu_opt = 1.0e-5,
            method_linear = KLUFactorization()
        )
    end

    testval = 0.0
    for i in 2:length(tsol.t)
        ui = view(tsol, 2, :, i)
        Δt = tsol.t[i] - tsol.t[i - 1]
        testval += sum(view(ui, subgrid2)) * Δt
    end

    if !isnothing(Plotter)
        for i in 2:length(tsol.t)
            plot_timestep(tsol.u[i], tsol.t[i])
        end
    end
    return testval
end

using Test

function runtests()
    testval = 0.06922262169719146
    testvaldiffeq = 0.06889809741891571
    @test main(; unknown_storage = :sparse, rely_on_corrections = false, assembly = :edgewise) ≈ testval
    @test main(; unknown_storage = :dense, rely_on_corrections = false, assembly = :edgewise) ≈ testval
    @test main(; unknown_storage = :sparse, rely_on_corrections = true, assembly = :edgewise) ≈ testval
    @test main(; unknown_storage = :dense, rely_on_corrections = true, assembly = :edgewise) ≈ testval
    @test main(; unknown_storage = :sparse, rely_on_corrections = false, assembly = :cellwise) ≈ testval
    @test main(; unknown_storage = :dense, rely_on_corrections = false, assembly = :cellwise) ≈ testval
    @test main(; unknown_storage = :sparse, rely_on_corrections = true, assembly = :cellwise) ≈ testval
    @test main(; unknown_storage = :dense, rely_on_corrections = true, assembly = :cellwise) ≈ testval


    @test main(; diffeq = true, unknown_storage = :sparse, rely_on_corrections = false, assembly = :edgewise) ≈ testvaldiffeq
    @test main(; diffeq = true, unknown_storage = :dense, rely_on_corrections = false, assembly = :edgewise) ≈ testvaldiffeq
    @test main(; diffeq = true, unknown_storage = :sparse, rely_on_corrections = true, assembly = :edgewise) ≈ testvaldiffeq
    @test main(; diffeq = true, unknown_storage = :dense, rely_on_corrections = true, assembly = :edgewise) ≈ testvaldiffeq
    @test main(; diffeq = true, unknown_storage = :sparse, rely_on_corrections = false, assembly = :cellwise) ≈ testvaldiffeq
    @test main(; diffeq = true, unknown_storage = :dense, rely_on_corrections = false, assembly = :cellwise) ≈ testvaldiffeq
    @test main(; diffeq = true, unknown_storage = :sparse, rely_on_corrections = true, assembly = :cellwise) ≈ testvaldiffeq
    @test main(; diffeq = true, unknown_storage = :dense, rely_on_corrections = true, assembly = :cellwise) ≈ testvaldiffeq
    return nothing
end

end
