#=

# 430: Parameter Derivatives (stationary)
 ([source code](@__SOURCE_URL__))

Explore different ways to calculate sensitivities.
This is still experimental. 
=#
module Example430_ParameterDerivativesStationary

using VoronoiFVM, ExtendableGrids
using GridVisualize
using ExtendableSparse
using ExtendableSparse: ILUZeroPreconBuilder
using ForwardDiff, DiffResults
using SparseDiffTools, SparseArrays
using ILUZero, LinearSolve

"""
    f(P)

Parameter dependent function which creates system and solves it
"""
function f(P; n = 10)
    p = P[1]

    valuetype = typeof(p)
    nspecies = 1
    ispec = 1

    function flux!(f, u, edge, data)
        f[1] = (1 + p) * (u[1, 1]^2 - u[1, 2]^2)
        return nothing
    end

    function r!(f, u, edge, data)
        f[1] = p * u[1]^5
        return nothing
    end

    function bc!(f, u, node, data)
        boundary_dirichlet!(f, u, node, ispec, 1, 0.0)
        boundary_dirichlet!(f, u, node, ispec, 3, p)
        return nothing
    end

    X = collect(0:(1.0 / n):1)
    grid = simplexgrid(X, X)
    sys = VoronoiFVM.System(
        grid; valuetype, species = [1], flux = flux!, reaction = r!,
        bcondition = bc!
    )
    tff = VoronoiFVM.TestFunctionFactory(sys)
    tfc = testfunction(tff, [1], [3])
    sol = solve(sys; inival = 0.5)
    return [integrate(sys, tfc, sol)[1]]
end

"""
    runf(;Plotter, n=10)

Run parameter series, plot f(p), df(p).
For each p,create a new system. Use VoronoiFVM with dual numbers. Pass parameters via closure.
"""
function runf(; Plotter = nothing, n = 10)
    P = 0.1:0.05:2
    dresult = DiffResults.JacobianResult(ones(1))
    F = zeros(0)
    DF = zeros(0)
    ff(p) = f(p; n)
    @time for p in P
        ForwardDiff.jacobian!(dresult, ff, [p])
        push!(F, DiffResults.value(dresult)[1])
        push!(DF, DiffResults.jacobian(dresult)[1])
    end
    vis = GridVisualizer(; Plotter, legend = :lt)
    scalarplot!(vis, P, F; color = :red, label = "f")
    scalarplot!(vis, P, DF; color = :blue, label = "df", clear = false, show = true)
    return sum(DF)
end

function fluxg!(f, u, edge, data)
    f[1] = (1 + data.p) * (u[1, 1]^2 - u[1, 2]^2)
    return nothing
end

function rg!(f, u, edge, data)
    f[1] = data.p * u[1]^5
    return nothing
end

function bcg!(f, u, node, data)
    boundary_dirichlet!(f, u, node, 1, 1, 0.0)
    boundary_dirichlet!(f, u, node, 1, 3, data.p)
    return nothing
end

Base.@kwdef mutable struct MyData{Tv}
    p::Tv = 1.0
end

"""
    rung(;Plotter, n=10)

Same as runf, but keep one system pass parameters via data.
"""
function rung(; Plotter = nothing, method_linear = SparspakFactorization(), n = 10)
    X = collect(0:(1.0 / n):1)
    grid = simplexgrid(X, X)

    # ugly but simple. By KISS we should first provide this way.

    sys = nothing
    data = nothing
    tfc = nothing

    function g(P)
        Tv = eltype(P)
        if isnothing(sys)
            data = MyData(one(Tv))
            sys = VoronoiFVM.System(
                grid; valuetype = Tv, species = [1], flux = fluxg!,
                reaction = rg!, bcondition = bcg!, data,
                unknown_storage = :dense
            )
            tff = VoronoiFVM.TestFunctionFactory(sys)
            tfc = testfunction(tff, [1], [3])
        end
        data.p = P[1]
        sol = solve(sys; inival = 0.5, method_linear)
        return [integrate(sys, tfc, sol)[1]]
    end

    dresult = DiffResults.JacobianResult(ones(1))

    P = 0.1:0.05:2
    G = zeros(0)
    DG = zeros(0)
    @time for p in P
        ForwardDiff.jacobian!(dresult, g, [p])
        push!(G, DiffResults.value(dresult)[1])
        push!(DG, DiffResults.jacobian(dresult)[1])
    end

    vis = GridVisualizer(; Plotter, legend = :lt)
    scalarplot!(vis, P, G; color = :red, label = "g")
    scalarplot!(vis, P, DG; color = :blue, label = "dg", clear = false, show = true)
    return sum(DG)
end

#########################################################################

function fluxh!(f, u, edge, data)
    p = parameters(u)[1]
    f[1] = (1 + p) * (u[1, 1]^2 - u[1, 2]^2)
    return nothing
end

function rh!(f, u, edge, data)
    p = parameters(u)[1]
    f[1] = p * u[1]^5
    return nothing
end

function bch!(f, u, node, data)
    p = parameters(u)[1]
    boundary_dirichlet!(f, u, node, 1, 1, 0.0)
    boundary_dirichlet!(f, u, node, 1, 3, p)
    return nothing
end

"""
    runh(;Plotter, n=10)

Same as runf, but use "normal" calculation (don't solve in dual numbers), and calculate dudp during
main assembly loop. 

This needs quite a bit of additional implementation + corresponding API and still lacks local assembly of the 
measurement derivative (when using testfunction based calculation) when calculating current.
"""
function runh(; Plotter = nothing, n = 10)
    X = collect(0:(1.0 / n):1)
    grid = simplexgrid(X, X)

    sys = VoronoiFVM.System(
        grid; species = [1], flux = fluxh!, reaction = rh!,
        bcondition = bch!, unknown_storage = :dense, nparams = 1
    )
    tff = VoronoiFVM.TestFunctionFactory(sys)
    tfc = testfunction(tff, [1], [3])

    function measp(params, u)
        Tp = eltype(params)
        up = Tp.(u)
        return integrate(sys, tfc, up; params = params)[1]
    end

    params = [0.0]

    function mymeas!(meas, U)
        u = reshape(U, sys)
        meas[1] = integrate(sys, tfc, u; params)[1]
        return nothing
    end

    dp = 0.05
    P = 0.1:dp:2
    state = VoronoiFVM.SystemState(sys)
    U0 = solve!(state; inival = 0.5, params = [P[1]])

    ndof = num_dof(sys)
    colptr = [i for i in 1:(ndof + 1)]
    rowval = [1 for i in 1:ndof]
    nzval = [1.0 for in in 1:ndof]
    ∂m∂u = zeros(1, ndof)
    colors = matrix_colors(∂m∂u)

    H = zeros(0)
    DH = zeros(0)
    DHx = zeros(0)
    m = zeros(1)

    @time for p in P
        params[1] = p
        sol = solve!(state; inival = 0.5, params)

        mymeas!(m, sol)
        push!(H, m[1])

        # this one is expensive - we would need to assemble this jacobian via local calls
        forwarddiff_color_jacobian!(∂m∂u, mymeas!, vec(sol); colorvec = colors)

        # need to have the full derivative of m vs p
        ∂m∂p = ForwardDiff.gradient(p -> measp(p, sol), params)

        dudp = state.matrix \ vec(state.dudp[1])
        dmdp = -∂m∂u * dudp + ∂m∂p
        push!(DH, dmdp[1])
    end

    vis = GridVisualizer(; Plotter, legend = :lt)
    scalarplot!(vis, P, H; color = :red, label = "h")
    scalarplot!(vis, P, DH; color = :blue, label = "dh", clear = false, show = true)
    return sum(DH)
end

using Test
function runtests()
    testval = 489.3432830184927
    @test runf() ≈ testval
    @test rung() ≈ testval
    @test runh() ≈ testval
    return nothing
end

end
