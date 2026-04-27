#=

# 006: Newton Solver History
([source code](@__SOURCE_URL__))

This example demonstrates how to inspect the linear solver timing
information stored in the Newton solver history.

=#

module DevEx006_NewtonSolverHistory

using VoronoiFVM
using ExtendableGrids
using LinearSolve
using ExtendableSparse: ILUZeroPreconBuilder

using Test


function main(; n = 10, assembly = :edgewise, kwargs...)
    h = 1.0 / convert(Float64, n)
    X = collect(0.0:h:1.0)
    Y = collect(0.0:h:1.0)

    grid = simplexgrid(X, Y)

    eps = 1.0e-2

    function reaction(f, u, node, data)
        f[1] = u[1]^2
        return nothing
    end

    function flux(f, u, edge, data)
        f[1] = eps * (u[1, 1]^2 - u[1, 2]^2)
        return nothing
    end

    function source(f, node, data)
        x1 = node[1] - 0.5
        x2 = node[2] - 0.5
        f[1] = exp(-20.0 * (x1^2 + x2^2))
        return nothing
    end

    function bcondition(f, u, node, data)
        boundary_dirichlet!(f, u, node; species = 1, region = 2, value = 1.0)
        boundary_dirichlet!(f, u, node; species = 1, region = 4, value = 1.0)
        return nothing
    end

    sys = VoronoiFVM.System(
        grid; reaction, flux, source, bcondition, assembly,
        species = [1]
    )

    @info "UMFPACK:"
    sol = solve(
        sys;
        inival = 0.5,
        method_linear = LinearSolve.UMFPACKFactorization(),
        log = true,
        kwargs...
    )
    sol_history = history(sol)
    tlinsolve = sol_history.tlinsolve
    tlinsolve_setup = sol_history.tlinsolve_setup
    tlinsolve_solve = sol_history.tlinsolve_solve
    @test tlinsolve ≥ tlinsolve_setup + tlinsolve_solve

    @info "Krylov-ilu0:"
    sol = solve(
        sys;
        inival = 0.5,
        method_linear = KrylovJL_BICGSTAB(precs = ILUZeroPreconBuilder()),
        log = true,
        kwargs...
    )
    sol_history = history(sol)
    tlinsolve = sol_history.tlinsolve
    tlinsolve_setup = sol_history.tlinsolve_setup
    tlinsolve_solve = sol_history.tlinsolve_solve
    @test tlinsolve ≥ tlinsolve_setup + tlinsolve_solve

    return @test tlinsolve ≥ tlinsolve_setup + tlinsolve_solve
end

function runtests()
    @testset "edgewise" begin
        main(; assembly = :edgewise)
    end
    @testset "cellwise" begin
        main(; assembly = :cellwise)
    end
    return nothing
end

end
