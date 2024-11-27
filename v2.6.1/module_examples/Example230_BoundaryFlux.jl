#=

# 103: Boundary flux
([source code](@__SOURCE_URL__))

We consider two test problems.

Testproblem A: Consider in $\Omega_1=(0,1)$
```math
- d_1 \Delta u_1  + k_1 u_1 = c_1
```
in  with homogeneous Neumann boundary conditions.

Testproblem B: Consider in $\Omega_2=(0,1) x (0, 1) $ 
```math
- d_2 \Delta u_2  + k_2 u_2 = c_2
```
in  with homogeneous Neumann boundary conditions
and at the right boundary, i.e. $ {1} x (0, 1) $
```math
- d_b \Delta v  + k_b v = c_b.
```
If d_1 = d_b, k_1 = k_b and c_1 = c_b, then u and v should coincide.
=#

module Example230_BoundaryFlux

using VoronoiFVM
using ExtendableGrids
using GridVisualize

function main(;
        n = 2 * 10, # n musst be an even number
        d1 = 5.0, db = 5.0, # prefactors (before diffusive part)
        kmax = 2.0, cmax = 3.0,
        Plotter = nothing,
        unknown_storage = :sparse, assembly = :edgewise
    )

    ###########################################################################
    ######################          1D problem           ######################
    ###########################################################################

    ispec_1D = 1
    bulk_1D = 1

    X = range(0.0; stop = 1.0, length = n)
    length_x = length(X)
    length_x_half = Int(length_x / 2)

    grid_1D = simplexgrid(X)

    k1 = zeros(length_x)
    c1 = zeros(length_x)

    k1[1:length_x_half] .= kmax
    k1[(length_x_half + 1):length_x] .= 0.0  # prefactor before reactive part
    c1[1:length_x_half] .= 0.0
    c1[(length_x_half + 1):length_x] .= cmax # source term

    #### discretization functions ####

    function flux!(f, u, edge, data)
        f[1] = d1 * (u[1, 1] - u[1, 2])
        return nothing
    end

    function reaction!(f, u, node, data)
        f[1] = k1[node.index] * u[1]
        return nothing
    end

    function source!(f, node::VoronoiFVM.Node, data)
        f[1] = c1[node.index]
        return nothing
    end

    sys_1D = VoronoiFVM.System(
        grid_1D,
        VoronoiFVM.Physics(;
            flux = flux!, reaction = reaction!,
            source = source!
        )
    )

    # enable species in only region
    enable_species!(sys_1D, ispec_1D, [bulk_1D])

    ## Stationary solution of both problems
    sol_1D = solve(sys_1D; inival = 0)

    p = GridVisualizer(;
        Plotter = Plotter, layout = (2, 1), clear = true,
        resolution = (800, 500)
    )

    scalarplot!(
        p[1, 1], grid_1D, sol_1D[1, :]; show = true,
        title = "1D calculation (d1 = $d1, kmax = $kmax, cmax = $cmax)"
    )

    ###########################################################################
    ######################          2D problem           ######################
    ###########################################################################

    grid_2D = simplexgrid(X, X)

    ispec_2D = 1
    ispec_boundary = 2
    bulk_2D = 1
    active_boundary = 2

    # parameters for the bulk problem
    d2 = 1.0
    k2 = 1.0
    c2 = 1.0

    #### discretization functions for bulk species ####
    function flux2D!(f, u, edge, data)
        f[ispec_2D] = d2 * (u[ispec_2D, 1] - u[ispec_2D, 2])
        return nothing
    end

    function reaction2D!(f, u, node, data)
        f[ispec_2D] = k2 * u[ispec_2D]
        return nothing
    end

    function source2D!(f, node, data)
        f[ispec_2D] = c2
        return nothing
    end

    #### discretization functions for boundary species at active boundary ####
    function bflux!(f, u, bedge, data)
        if bedge.region == active_boundary
            f[ispec_boundary] = db * (u[ispec_boundary, 1] - u[ispec_boundary, 2])
        end
        return nothing
    end

    function breaction!(f, u, bnode, data)
        if bnode.region == active_boundary
            if bnode.coord[2, bnode.index] <= 0.5
                kb = kmax
            else
                kb = 0.0
            end

            f[ispec_boundary] = kb * u[ispec_boundary]
        end
        return nothing
    end

    function bsource!(f, bnode, data)
        if bnode.region == active_boundary
            if bnode.coord[2, bnode.index] <= 0.5
                cb = 0.0
            else
                cb = cmax
            end

            f[ispec_boundary] = cb
        end
        return nothing
    end

    sys_2D = VoronoiFVM.System(
        grid_2D,
        VoronoiFVM.Physics(;
            flux = flux2D!, reaction = reaction2D!,
            source = source2D!,
            bflux = bflux!, breaction = breaction!,
            bsource = bsource!
        );
        unknown_storage = unknown_storage, assembly = assembly
    )

    # enable species in only region
    enable_species!(sys_2D, ispec_2D, [bulk_2D])
    enable_boundary_species!(sys_2D, ispec_boundary, [active_boundary])

    sol_2D = solve(sys_2D; inival = 0)

    # this is for variable transformation, since we consider right outer boundary and want to transform to x-axis.
    function tran32!(a, b)
        a[1] = b[2]
        return nothing
    end

    # note that if adjusting active_boundary to 3 or 4, then transform needs to be deleted.
    bgrid_2D = subgrid(grid_2D, [active_boundary]; boundary = true, transform = tran32!)
    sol_bound = view(sol_2D[ispec_boundary, :], bgrid_2D)

    scalarplot!(
        p[2, 1], bgrid_2D, sol_bound; show = true, cellwise = true,
        title = "Active boundary in 2D (db = $db, kb = $kmax, cb = $cmax)"
    )

    errorsol = VoronoiFVM.norm(sys_1D, sol_bound - sol_1D', 2)

    return errorsol
end # main

using Test
function runtests()
    @test main(; unknown_storage = :dense, assembly = :edgewise) < 1.0e-14 &&
        main(; unknown_storage = :sparse, assembly = :edgewise) < 1.0e-14 &&
        main(; unknown_storage = :dense, assembly = :cellwise) < 1.0e-14 &&
        main(; unknown_storage = :sparse, assembly = :cellwise) < 1.0e-14
    return nothing
end

end # module
