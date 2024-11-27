#=
 # 002: check edge reaction
=#

module Example002_EdgeReaction

using Printf
using VoronoiFVM
using ExtendableGrids
using ExtendableSparse
using GridVisualize
using LinearAlgebra
using SimplexGridFactory
using Triangulate

function main(; nref = 0, dim = 2, Plotter = nothing, verbose = "and", case = :compare_max, assembly = :edgewise)
    X = 0:(0.25 * 2.0^-nref):1
    i0::Int = 0
    i1::Int = 0
    if dim == 1
        grid = simplexgrid(X)
        i0 = 1
        i1 = 2
    elseif dim == 2
        b = SimplexGridBuilder(; Generator = Triangulate)
        p00 = point!(b, 0, 0)
        p10 = point!(b, 1, 0)
        p11 = point!(b, 1, 1)
        p01 = point!(b, 0, 1)
        pxx = point!(b, 0.3, 0.3)

        facetregion!(b, 1)
        facet!(b, p00, p10)
        facetregion!(b, 2)
        facet!(b, p10, p11)
        facetregion!(b, 3)
        facet!(b, p11, p01)
        facetregion!(b, 4)
        facet!(b, p01, p00)
        grid = simplexgrid(b; maxvolume = 0.01 * 4.0^(-nref))
        i0 = 1
        i1 = 3
    elseif dim == 3
        grid = simplexgrid(X, X, X)
        i0 = 5
        i1 = 6
    end

    function storage!(y, u, node, data)
        y[1] = u[1]
        return nothing
    end

    function flux!(y, u, edge, data)
        y[1] = u[1, 1] - u[1, 2]
        return nothing
    end

    # Three ways to give a constant reaction term. As a consequence,
    # these need to yield the same solution.
    # 1: classical node reaction, multiplied  by control volume size
    function reaction!(y, u, node, data)
        y[1] = -1
        return nothing
    end

    # 2: Edge reaction. Here we give it as a constant, and wie need
    # to turn the multiplication with σ/h into a multiplication with the
    # half diamond volume.
    #
    # Half diamond volume calculation
    #         /|\
    #        / | \
    #       /  |s \
    #       -------
    #         h
    #  A=s*h/2d . Our formfactor:  σ=s/h =>  A=σ*h^2
    #  - make transfer area to volume
    #
    #  τ=1/h    v= s*h/2d = σ*h^2/2d
    #
    function edgereaction!(y, u, edge, data)
        h = meas(edge)
        y[1] = -1 * h^2 / (2 * dim)
        return nothing
    end

    #
    # 3: "Joule heat:" |∇ϕ|^2=1 after 3.17 in Bradji/Herbin
    # Here we divide twice by "h" to get the constant squared gradient.
    # The multiplication with dim in 3.17 compensates the division
    # we had before
    ϕ = grid[Coordinates][1, :]

    function edgereaction2!(y, u, edge, data)
        ϕK = ϕ[edge.node[1]]
        ϕL = ϕ[edge.node[2]]
        y[1] = -(ϕK - ϕL) * (ϕK - ϕL) / 2
        return nothing
    end

    if case == :compare_max
        function bcondition!(y, u, node, data)
            boundary_dirichlet!(y, u, node; species = 1, region = 1, value = 0)
            boundary_dirichlet!(y, u, node; species = 1, region = 2, value = 0)
            boundary_dirichlet!(y, u, node; species = 1, region = 3, value = 0)
            boundary_dirichlet!(y, u, node; species = 1, region = 4, value = 0)
            boundary_dirichlet!(y, u, node; species = 1, region = 5, value = 0)
            boundary_dirichlet!(y, u, node; species = 1, region = 6, value = 0)
            return nothing
        end

        sys_noderea = VoronoiFVM.System(
            grid; bcondition = bcondition!, flux = flux!,
            reaction = reaction!, storage = storage!,
            species = [1], is_linear = true, assembly
        )
        sys_edgerea = VoronoiFVM.System(
            grid; bcondition = bcondition!, flux = flux!,
            edgereaction = edgereaction!, storage = storage!,
            species = [1], is_linear = true, assembly
        )
        sys_edgerea2 = VoronoiFVM.System(
            grid; bcondition = bcondition!, flux = flux!,
            edgereaction = edgereaction2!, storage = storage!,
            species = [1], is_linear = true, assembly
        )

        sol_noderea = solve(sys_noderea; verbose)
        sol_edgerea = solve(sys_edgerea; verbose)
        sol_edgerea2 = solve(sys_edgerea2; verbose)

        vis = GridVisualizer(; Plotter, layout = (2, 2))
        scalarplot!(
            vis[1, 1], grid, sol_noderea[1, :]; title = "node reaction",
            colormap = :hot
        )
        scalarplot!(vis[2, 1], grid, sol_edgerea[1, :]; title = "edgerea1", colormap = :hot)
        scalarplot!(
            vis[1, 2], grid, sol_edgerea2[1, :]; title = "edgerea2",
            colormap = :hot
        )

        reveal(vis)
        return maximum.([sol_noderea, sol_edgerea, sol_edgerea2])
    end

    if case == :compare_flux
        function bcondition2!(y, u, node, data)
            boundary_dirichlet!(y, u, node; species = 1, region = i1, value = 0)
            return nothing
        end

        sys2_noderea = VoronoiFVM.System(
            grid; bcondition = bcondition2!, flux = flux!,
            reaction = reaction!, storage = storage!,
            species = [1], is_linear = true
        )
        sys2_edgerea = VoronoiFVM.System(
            grid; bcondition = bcondition2!, flux = flux!,
            edgereaction = edgereaction!, storage = storage!,
            species = [1], is_linear = true
        )
        sys2_edgerea2 = VoronoiFVM.System(
            grid; bcondition = bcondition2!, flux = flux!,
            edgereaction = edgereaction2!, storage = storage!,
            species = [1], is_linear = true
        )

        sol2_noderea = solve(sys2_noderea; verbose)
        sol2_edgerea = solve(sys2_edgerea; verbose)
        sol2_edgerea2 = solve(sys2_edgerea2; verbose)

        tfac2_noderea = TestFunctionFactory(sys2_noderea)
        tfc2_noderea = testfunction(tfac2_noderea, [i0], [i1])

        tfac2_edgerea = TestFunctionFactory(sys2_edgerea)
        tfc2_edgerea = testfunction(tfac2_edgerea, [i0], [i1])

        tfac2_edgerea2 = TestFunctionFactory(sys2_edgerea2)
        tfc2_edgerea2 = testfunction(tfac2_edgerea2, [i0], [i1])

        vis = GridVisualizer(; Plotter, layout = (2, 2))
        scalarplot!(
            vis[1, 1], grid, sol2_noderea[1, :]; title = "node reaction",
            colormap = :hot
        )
        scalarplot!(
            vis[2, 1], grid, sol2_edgerea[1, :]; title = "edgerea1",
            colormap = :hot
        )
        scalarplot!(
            vis[1, 2], grid, sol2_edgerea2[1, :]; title = "edgerea2",
            colormap = :hot
        )
        reveal(vis)

        I_noderea = integrate(sys2_noderea, tfc2_noderea, sol2_noderea)
        I_edgerea = integrate(sys2_edgerea, tfc2_edgerea, sol2_edgerea)
        I_edgerea2 = integrate(sys2_edgerea2, tfc2_edgerea2, sol2_edgerea2)

        return I_noderea, I_edgerea, I_edgerea2
    end
end

using Test
function runtests()
    res = fill(false, 3)
    for dim in 1:3
        result_max = main(; case = :compare_max, assembly = :cellwise)
        result_flux = main(; case = :compare_flux, assembly = :cellwise)
        res[dim] = isapprox(result_max[1], result_max[2]; atol = 1.0e-6) &&
            isapprox(result_max[1], result_max[3]; atol = 1.0e-3) &&
            isapprox(result_flux[1], result_flux[2]; atol = 1.0e-10) &&
            isapprox(result_flux[1], result_flux[3]; atol = 1.0e-10)
    end
    res1 = all(a -> a, res)

    res = fill(false, 3)
    for dim in 1:3
        result_max = main(; case = :compare_max, assembly = :edgwise)
        result_flux = main(; case = :compare_flux, assembly = :edgwise)
        res[dim] = isapprox(result_max[1], result_max[2]; atol = 1.0e-6) &&
            isapprox(result_max[1], result_max[3]; atol = 1.0e-3) &&
            isapprox(result_flux[1], result_flux[2]; atol = 1.0e-10) &&
            isapprox(result_flux[1], result_flux[3]; atol = 1.0e-10)
    end
    res2 = all(a -> a, res)

    @test res1 && res2
    return nothing
end

end
