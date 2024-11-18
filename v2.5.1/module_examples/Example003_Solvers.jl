#=
# 003: New linear solver API
([source code](@__SOURCE_URL__))
=#

module Example003_Solvers

## under development

using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize
using LinearSolve
using ExtendableSparse
using ExtendableSparse: ILUZeroPreconBuilder, JacobiPreconBuilder, SmoothedAggregationPreconBuilder
using SparseArrays
using AMGCLWrap
using AlgebraicMultigrid
using LinearAlgebra


using Test


function main(; n = 10, Plotter = nothing, assembly = :edgwwise, kwargs...)
    h = 1.0 / convert(Float64, n)
    X = collect(0.0:h:1.0)
    Y = collect(0.0:h:1.0)

    grid = simplexgrid(X, Y)
    nn = num_nodes(grid)

    eps = 1.0e-2

    function reaction(f, u, node, data)
        f[1] = u[1]^2
    end

    function flux(f, u, edge, data)
        f[1] = eps * (u[1, 1]^2 - u[1, 2]^2)
    end

    function source(f, node, data)
        x1 = node[1] - 0.5
        x2 = node[2] - 0.5
        f[1] = exp(-20.0 * (x1^2 + x2^2))
    end

    function storage(f, u, node, data)
        f[1] = u[1]
    end

    function bcondition(f, u, node, data)
        boundary_dirichlet!(f,
                            u,
                            node;
                            species = 1,
                            region = 2,
                            value = ramp(node.time; dt = (0, 0.1), du = (0, 1)))
        boundary_dirichlet!(f,
                            u,
                            node;
                            species = 1,
                            region = 4,
                            value = ramp(node.time; dt = (0, 0.1), du = (0, 1)))
    end

    sys = VoronoiFVM.System(grid; reaction, flux, source, storage, bcondition, assembly,
                            species = [1])
    @info "UMFPACK:"
    umf_sol = solve(sys; inival = 0.5, method_linear = LinearSolve.UMFPACKFactorization(), kwargs...)

    @info "KLU:"
    sol = solve(sys; inival = 0.5, method_linear = LinearSolve.KLUFactorization(), kwargs...)
    @test norm(sol - umf_sol, Inf)<1.0e-7
    
    @info "Sparspak:"
    sol = solve(sys; inival = 0.5, method_linear = LinearSolve.SparspakFactorization(), kwargs...)
    @test norm(sol - umf_sol, Inf)<1.0e-7

    @info "Krylov-ilu0:"
    sol = solve(sys;
                inival = 0.5,
                method_linear = KrylovJL_BICGSTAB(precs=ILUZeroPreconBuilder()),
                kwargs...)
    @test norm(sol - umf_sol, Inf)<1.0e-7

    @info "Krylov-block1"
    precs=BlockPreconBuilder(;precs=ILUZeroPreconBuilder(), partitioning= A->  [1:(size(A,1) ÷ 2), (size(A,1) ÷ 2 + 1):size(A,1)])
    sol = solve(sys;
                inival = 0.5,
                method_linear = KrylovJL_BICGSTAB(;precs),
                kwargs...)
    @test norm(sol - umf_sol, Inf)<1.0e-7

    @info "Krylov-block2"
    precs=BlockPreconBuilder(;precs=ILUZeroPreconBuilder(), partitioning= A->   [1:2:size(A,1), 2:2:size(A,1)])
    sol = solve(sys;
                inival = 0.5,
                method_linear = KrylovJL_BICGSTAB(;precs),
                log=true,
                kwargs...)
    @test norm(sol - umf_sol, Inf)<1.0e-7
    
    @info "Krylov - delayed factorization:"
    sol = solve(sys;
                inival = 0.5,
                method_linear = KrylovJL_BICGSTAB(;precs=LinearSolvePreconBuilder(SparspakFactorization())),
                keepcurrent_linear =false,
                log=true,
                kwargs...)
    @test norm(sol - umf_sol, Inf)<1.0e-7
    @test summary(history(sol)).factorizations == 1
    
    @info "Krylov - jacobi:"
    sol = solve(sys;
                inival = 0.5,
                method_linear = KrylovJL_BICGSTAB(;precs=JacobiPreconBuilder()),
                keepcurrent_linear = true, log=true,
                kwargs...)
    @test norm(sol - umf_sol, Inf)<1.0e-7
    @test summary(history(sol)).factorizations > 1
    
    @info "Krylov - SA_AMG:"
    sol = solve(sys;
                inival = 0.5,
                method_linear = KrylovJL_BICGSTAB(;precs=SmoothedAggregationPreconBuilder()),
                keepcurrent_linear = true,
                kwargs...)
    @test norm(sol - umf_sol, Inf)<1.0e-7

    @info "Krylov - AMGCL_AMG:"
    sol = solve(sys;
                inival = 0.5,
                method_linear = KrylovJL_BICGSTAB(;precs=AMGPreconBuilder()),
                keepcurrent_linear = true,
                kwargs...)
    @test norm(sol - umf_sol, Inf)<1.0e-7

    @info "Krylov - AMGnCL_RLX:"
    sol = solve(sys;
                inival = 0.5,
                method_linear = KrylovJL_BICGSTAB(;precs=RLXPreconBuilder()),
                keepcurrent_linear = true,
                kwargs...)
    @test norm(sol - umf_sol, Inf)<1.0e-7
end

function runtests()
    @testset "edgewise" begin
        main(; assembly = :edgewise)
    end
    @testset "cellwise" begin
        main(; assembly = :cellwise)
    end
end
end
