module test_cachesol
using ExtendableGrids
using VoronoiFVM
using LinearSolve
using ExtendableSparse
using LinearAlgebra
using SparseArrays
using Test

function flux!(f, u, edge, data)
    f[1] = u[1, 1] - u[1, 2]
    return nothing
end

function bcondition!(f, u, bnode, data)
    boundary_dirichlet!(f, u, bnode; species = 1, region = 1, value = 0.0)
    boundary_dirichlet!(f, u, bnode; species = 1, region = 2, value = 1.0)
    return nothing
end

function setup_linear_solve(n)
    # define system, initialize system matrix
    X = linspace(0.0, 1.0, n)
    system = VoronoiFVM.System(X; flux = flux!, bcondition = bcondition!, species = [1])
    F = unknowns(system; inival = 0.0)
    # assemble system matrix A and the residual F of A*U=RHS in U=0, i.e.
    # F = A*0 - RHS = -RHS
    F, A = VoronoiFVM.evaluate_residual_and_jacobian(system, F)

    p = LinearProblem(A, -VoronoiFVM.dofs(F))
    linear_cache = init(p, LinearSolve.UMFPACKFactorization())

    return X, linear_cache
end

function solve_problem(cache)
    sol = LinearSolve.solve!(cache)
    return sol
end

function main(; n = 10)
    X, cache = setup_linear_solve(n)
    approx_sol = solve_problem(cache)
    true_sol = map((x) -> (x), X)
    @test norm(approx_sol .- true_sol, Inf) ≤ 3.0e-16
    return nothing
end

function runtests()
    main()
    return nothing
end
end
