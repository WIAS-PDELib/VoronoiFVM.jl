module test_norms
using VoronoiFVM
using ExtendableGrids
using Test
using LinearAlgebra

function grid(X, dim)
    if dim == 1
        return simplexgrid(X)
    elseif dim == 2
        return simplexgrid(X, X)
    else
        return simplexgrid(X, X, X)
    end
    return nothing
end

function test_solint(; dim = 2, c = 1.0, assembly = :edgewise, h = 0.1)
    X = 0:h:1
    g = grid(X, dim)
    sys = VoronoiFVM.System(g; species = [1], assembly)
    VoronoiFVM._complete!(sys)
    u = map(c, sys)
    return VoronoiFVM.integrate(sys, u)[1, 1]
end


function test_edgeint(; dim = 2, c = 1.0, assembly = :edgewise, h = 0.1)
    X = 0:h:1
    g = grid(X, dim)
    sys = VoronoiFVM.System(g; species = [1], assembly)
    VoronoiFVM._complete!(sys)
    u = map(c, sys)
    function f(y, u, edge, data = nothing)
        y[1] = 0.5 * (u[1, 1] + u[2, 1])
        return nothing
    end
    return VoronoiFVM.edgeintegrate(sys, f, u)[1, 1]
end

function test_const(; dim = 2, c = 1.0, nrm = l2norm, assembly = :edgewise, h = 0.1)
    X = 0:h:1
    g = grid(X, dim)
    sys = VoronoiFVM.System(g; species = [1], assembly)
    VoronoiFVM._complete!(sys)
    F = map(c, sys)
    return nrm(sys, F)
end

function test_lin(; dim = 2, c = 1.0, nrm = h1seminorm, assembly = :edgewise, h = 0.1)
    X = 0:h:1
    g = grid(X, dim)
    sys = VoronoiFVM.System(g; species = [1], assembly)
    VoronoiFVM._complete!(sys)
    F = unknowns(sys)
    F[1, :] .= map((x...) -> (c * sum(x) / sqrt(dim)), g)
    return nrm(sys, F)
end

function test_transient(;
        dim = 2,
        c = 1.0,
        nrm = l2h1seminorm,
        assembly = :edgewise,
        h = 0.1,
    )
    X = 0:h:1
    g = grid(X, dim)
    sys = VoronoiFVM.System(g; species = [1], assembly)
    VoronoiFVM._complete!(sys)
    F = unknowns(sys)
    F[1, :] .= map((x...) -> (c * sum(x) / sqrt(dim)), g)
    U = TransientSolution(X[1], F)
    for i in 2:length(X)
        append!(U, X[i], F)
    end

    return nrm(sys, U)
end

function test_transient_exp(;
        dim = 2,
        c = 1.0,
        nrm = l2h1seminorm,
        time_dep = t -> 1,
        assembly = :edgewise,
        h = 0.1,
    )
    X = 0:h:1
    g = grid(X, dim)
    sys = VoronoiFVM.System(g; species = [1], assembly)
    VoronoiFVM._complete!(sys)
    F = unknowns(sys)
    F[1, :] .= map((x...) -> c * sum(x) / sqrt(dim), g)
    U = TransientSolution(X[1], F)
    for i in 2:length(X)
        append!(U, X[i], time_dep(X[i]) * F)
    end

    return nrm(sys, U)
end

function test_nodevol(; dim = 1, h = 0.1, assembly = :edgewise)
    X = 0:h:1
    g = grid(X, dim)
    sys = VoronoiFVM.System(g; species = [1], assembly)
    return sum(nodevolumes(sys))
end

function runtests()
    for assembly in (:edgewise, :cellwise)
        for dim in 1:3
            for c in [0.5, 1, 2.0]
                @test test_solint(; dim, c, assembly) ≈ c
                @test test_edgeint(; dim, c, assembly) ≈ c
                @test test_const(; nrm = l2norm, dim, c, assembly) ≈ c
                @test test_const(; nrm = h1seminorm, dim, c, assembly) ≈ 0
                @test test_lin(; nrm = h1seminorm, dim, c, assembly) ≈ c
                @test test_transient(; nrm = l2h1seminorm, dim, c, assembly) ≈ c

                # thx  Theo Belin (@Theleb97) foe these:
                nrm = VoronoiFVM.l2h1seminorm
                @test test_transient_exp(; time_dep = t -> 1, nrm, dim, c, assembly) ≈ c
                @test test_transient_exp(; time_dep = t -> sqrt(t), nrm, dim, c, assembly) ≈ c * 0.7416198487095663 # c/sqrt(2)
                @test test_transient_exp(; time_dep = t -> t^(-1 / 4), nrm, dim, c, assembly) ≈ c * 1.2600710094548464 # c*sqrt(2)
                @test test_transient_exp(; time_dep = t -> t, nrm, dim, c, assembly) ≈ c * 0.620483682299542   # c/sqrt(3)
                @test test_transient_exp(; time_dep = t -> exp(-t / 2), nrm, dim, c, assembly) ≈ c * 0.7953912484980525 # c* sqrt((1-exp(-1)))
            end
            @test test_nodevol(; dim, assembly) ≈ 1.0
        end
    end
    return true
end

end
