#=

# 420: Discontinuous Quantities
 ([source code](@__SOURCE_URL__))

Test  jumping species and quantity handling

=#
module Example420_DiscontinuousQuantities
using Printf
using VoronoiFVM
using SparseArrays
using ExtendableGrids
using GridVisualize
using LinearAlgebra

function main(; N = 5, Plotter = nothing, unknown_storage = :sparse, assembly = :edgewise)
    XX = collect(0:0.1:1)
    xcoord = XX
    for i in 1:(N - 1)
        xcoord = glue(xcoord, XX .+ i)
    end
    grid2 = simplexgrid(xcoord)
    for i in 1:N
        cellmask!(grid2, [i - 1], [i], i)
    end
    for i in 1:(N - 1)
        bfacemask!(grid2, [i], [i], i + 2)
    end

    params = zeros(2, num_cellregions(grid2))
    for i in 1:num_cellregions(grid2)
        params[1, i] = i
        params[2, i] = 10 * i
    end

    system = VoronoiFVM.System(grid2; unknown_storage = unknown_storage, assembly = assembly)

    ## First, we introduce a continuous quantity which we name "cspec". Note that the "species number" can be assigned automatically if not given explicitly.
    cspec = ContinuousQuantity(system, 1:N; ispec = 1, id = 1)

    ## A discontinuous quantity can be introduced as well. by default, each reagion gets a new species number. This can be overwritten by the user.
    dspec = DiscontinuousQuantity(system, 1:N; regionspec = [2 + i % 2 for i in 1:N], id = 2)

    # check 1D array access with quantities
    carrierList = [cspec dspec]
    numberCarriers = length(carrierList)

    params2 = zeros(1, numberCarriers)

    for icc in carrierList
        params2[icc] = 2
    end

    for i in 1:numberCarriers
        @assert params2[i] == 2
    end

    # check 2D array access with quantities
    for i in 1:num_cellregions(grid2)
        @assert params[cspec, i] == i
        @assert params[dspec, i] == 10 * i
    end

    for i in 1:num_cellregions(grid2)
        params[cspec, i] = -i
        params[dspec, i] = -10 * i
    end

    for i in 1:num_cellregions(grid2)
        @assert params[1, i] == -i
        @assert params[2, i] == -10 * i
    end

    ##For both quantities, we define simple diffusion fluxes:

    function flux(f, u, edge, data)
        f[dspec] = u[dspec, 1] - u[dspec, 2]
        f[cspec] = u[cspec, 1] - u[cspec, 2]
        return nothing
    end

    d1 = 1
    q1 = 0.2

    function breaction(f, u, bnode, data)

        # left outer boundary value for dspec
        if bnode.region == 1
            f[dspec] = u[dspec] + 0.5
        end

        ## Define a thin layer interface condition for `dspec` and an interface source for `cspec`.
        if bnode.region > 2
            react = (u[dspec, 1] - u[dspec, 2]) / d1
            f[dspec, 1] = react
            f[dspec, 2] = -react
            f[cspec] = -q1 * u[cspec]
        end
        return nothing
    end

    physics!(
        system, VoronoiFVM.Physics(;
            flux = flux,
            breaction = breaction
        )
    )

    ## Set boundary conditions
    boundary_dirichlet!(system, dspec, 2, 0.1)
    boundary_dirichlet!(system, cspec, 1, 0.1)
    boundary_dirichlet!(system, cspec, 2, 1.0)
    subgrids = VoronoiFVM.subgrids(dspec, system)

    U = solve(system)

    dvws = views(U, dspec, subgrids, system)
    cvws = views(U, cspec, subgrids, system)
    vis = GridVisualizer(; resolution = (600, 300), Plotter = Plotter)
    for i in eachindex(dvws)
        scalarplot!(
            vis, subgrids[i], dvws[i]; flimits = (-0.5, 1.5), clear = false,
            color = :red
        )
        scalarplot!(
            vis, subgrids[i], cvws[i]; flimits = (-0.5, 1.5), clear = false,
            color = :green
        )
    end
    reveal(vis)
    I = integrate(system, system.physics.storage, U)
    return sum(I[dspec, :]) + sum(I[cspec, :])
end

using Test
function runtests()
    testval = 4.2
    @test main(; unknown_storage = :sparse, assembly = :edgewise) ≈ testval &&
        main(; unknown_storage = :dense, assembly = :edgewise) ≈ testval &&
        main(; unknown_storage = :sparse, assembly = :cellwise) ≈ testval &&
        main(; unknown_storage = :dense, assembly = :cellwise) ≈ testval
    return nothing
end

end
