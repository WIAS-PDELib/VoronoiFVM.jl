#=

# 422: Drift-Diffusion with Discontinuous and Interface Potentials
([source code](@__SOURCE_URL__))

Nondimensionalized semiconductor device equations (with artificial doping)
with Discontinuousquantities and additional Interfacequantities.

=#

module Example422_InterfaceQuantities

using VoronoiFVM
using ExtendableGrids
using GridVisualize

function main(;
        n = 5, Plotter = nothing, tend = 20.0, unknown_storage = :sparse,
        reactionN = 5.0e0, reactionP = 5.0e0, assembly = :edgewise
    )

    ################################################################################
    #### grid
    ################################################################################
    h1 = 1.0
    h2 = 1.0
    h_total = h1 + h2

    # region numbers
    region1 = 1
    region2 = 2
    regions = [region1, region2]
    numberOfRegions = length(regions)

    # boundary region numbers
    bregion1 = 1
    bregion2 = 2
    bjunction = 3

    coord_1 = collect(range(0.0; stop = h1, length = n))
    coord_2 = collect(range(h1; stop = h1 + h2, length = n))
    coord = glue(coord_1, coord_2)

    grid = simplexgrid(coord)

    # specify inner regions
    cellmask!(grid, [0.0], [h1], region1)
    cellmask!(grid, [h1], [h1 + h2], region2)

    # specify outer regions
    bfacemask!(grid, [0.0], [0.0], bregion1)
    bfacemask!(grid, [h_total], [h_total], bregion2)

    # inner interfaces
    bfacemask!(grid, [h1], [h1], bjunction)

    #gridplot(grid, Plotter = nothing, legend=:rt)

    ################################################################################
    #########  system
    ################################################################################

    sys = VoronoiFVM.System(grid; unknown_storage = unknown_storage, assembly = assembly)
    iphin = DiscontinuousQuantity(sys, 1:numberOfRegions; id = 1)
    iphip = DiscontinuousQuantity(sys, 1:numberOfRegions; id = 2)
    iphinb = InterfaceQuantity(sys, [bjunction]; id = 3)
    iphipb = InterfaceQuantity(sys, [bjunction]; id = 4)
    ipsi = ContinuousQuantity(sys, 1:numberOfRegions; id = 5)

    NA = [10.0, 0.0]
    ND = [0.0, 10.0]

    function storage!(f, u, node, data)
        etan = -((u[iphin] - u[ipsi]))
        etap = ((u[iphip] - u[ipsi]))

        f[iphin] = -exp(etan)
        f[iphip] = exp(etap)

        f[ipsi] = 0.0
        return nothing
    end

    function reaction!(f, u, node, data)
        etan = -((u[iphin] - u[ipsi]))
        etap = ((u[iphip] - u[ipsi]))

        f[ipsi] = -(ND[node.region] - exp(etan) + exp(etap) - NA[node.region])
        ########################
        r0 = 1.0e-4
        recomb = r0 * exp(etan) * exp(etap)

        f[iphin] = -recomb
        f[iphip] = recomb
        return nothing
    end

    function flux!(f, u, node, data)
        f[ipsi] = -(u[ipsi, 2] - u[ipsi, 1])

        ########################
        bp, bm = fbernoulli_pm(-(u[ipsi, 2] - u[ipsi, 1]))

        etan1 = -((u[iphin, 1] - u[ipsi, 1]))
        etap1 = ((u[iphip, 1] - u[ipsi, 1]))

        etan2 = -((u[iphin, 2] - u[ipsi, 2]))
        etap2 = ((u[iphip, 2] - u[ipsi, 2]))

        f[iphin] = (bm * exp(etan2) - bp * exp(etan1))
        f[iphip] = -(bp * exp(etap2) - bm * exp(etap1))
        return nothing
    end

    function breaction!(f, u, bnode, data)
        if bnode.region == bjunction
            # left values
            nleft = exp(-((u[iphin, 1] - u[ipsi])))
            pleft = exp(((u[iphip, 1] - u[ipsi])))

            # interface species
            n_interf = exp(-((u[iphinb] - u[ipsi])))
            p_interf = exp(((u[iphipb] - u[ipsi])))

            # right values
            nright = exp(-((u[iphin, 2] - u[ipsi])))
            pright = exp(((u[iphip, 2] - u[ipsi])))
            ################

            # left and right reaction for n
            f[iphin, 1] = reactionN * (nleft - n_interf)
            f[iphin, 2] = reactionN * (nright - n_interf)

            # left and right reaction for p
            f[iphip, 1] = reactionP * (pleft - p_interf)
            f[iphip, 2] = reactionP * (pright - p_interf)

            # interface species reaction
            f[iphinb] = -(f[iphin, 1] + f[iphin, 2])
            f[iphipb] = -(f[iphip, 1] + f[iphip, 2])
        end
        return nothing
    end

    function bstorage!(f, u, bnode, data)
        f[ipsi] = 0.0

        if bnode.region == bjunction
            etan = -((u[iphinb] - u[ipsi]))
            etap = ((u[iphipb] - u[ipsi]))

            f[iphinb] = -exp(etan)
            f[iphipb] = exp(etap)
        end
        return nothing
    end

    physics!(
        sys,
        VoronoiFVM.Physics(;
            flux = flux!,
            storage = storage!,
            reaction = reaction!,
            breaction = breaction!,
            bstorage = bstorage!
        )
    )

    boundary_dirichlet!(sys, iphin, bregion1, 4.0)
    boundary_dirichlet!(sys, iphip, bregion1, 4.0)
    boundary_dirichlet!(sys, ipsi, bregion1, 0.0)
    boundary_dirichlet!(sys, iphin, bregion2, 0.0)
    boundary_dirichlet!(sys, iphip, bregion2, 0.0)
    boundary_dirichlet!(sys, ipsi, bregion2, 5.0)

    ################################################################################
    #########  time loop
    ################################################################################

    ## Create a solution array
    sol = unknowns(sys)
    sol .= 0.0
    t0 = 1.0e-6
    ntsteps = 10
    tvalues = range(t0; stop = tend, length = ntsteps)

    for istep in 2:ntsteps
        t = tvalues[istep]       # Actual time
        Δt = t - tvalues[istep - 1] # Time step size

        #println("Δt = ", t)

        sol = solve(sys; inival = sol, tstep = Δt)
    end # time loop

    ################################################################################
    #########  Bias Loop
    ################################################################################
    biasval = range(0; stop = 2.0, length = 10)
    Idspec = zeros(0)

    for Δu in biasval
        boundary_dirichlet!(sys, iphin, bregion1, 4.0 + Δu)
        boundary_dirichlet!(sys, iphip, bregion1, 4.0 + Δu)
        boundary_dirichlet!(sys, ipsi, bregion1, 0.0 + Δu)

        #println("Δu = ", Δu)

        sol = solve(sys; inival = sol)

        ## get current
        factory = TestFunctionFactory(sys)
        tf = testfunction(factory, [1], [2])
        I = integrate(sys, tf, sol)

        val = 0.0
        for ii in 1:(length(I) - 1)
            val = val + I[ii]
        end

        push!(Idspec, val)
    end # bias loop

    ################################################################################
    #########  Plotting
    ################################################################################

    vis = GridVisualizer(; Plotter = Plotter, layout = (2, 1))

    subgridBulk = subgrids(iphin, sys)
    phin_sol = views(sol, iphin, subgridBulk, sys)
    phip_sol = views(sol, iphip, subgridBulk, sys)
    psi_sol = views(sol, ipsi, subgridBulk, sys)

    for i in eachindex(phin_sol)
        scalarplot!(vis[1, 1], subgridBulk[i], phin_sol[i]; clear = false, color = :green)
        scalarplot!(vis[1, 1], subgridBulk[i], phip_sol[i]; clear = false, color = :red)
        scalarplot!(vis[1, 1], subgridBulk[i], psi_sol[i]; clear = false, color = :blue)
    end

    scalarplot!(vis[2, 1], biasval, Idspec; clear = false, color = :red)

    bgrid = subgrids(iphinb, sys)
    sol_bound = views(sol, iphinb, bgrid, sys)

    return sol_bound[1]
end # main

using Test
function runtests()
    testval = 0.35545473758267826
    @test main(; unknown_storage = :sparse, assembly = :edgewise) ≈ testval &&
        main(; unknown_storage = :dense, assembly = :edgewise) ≈ testval &&
        main(; unknown_storage = :sparse, assembly = :cellwise) ≈ testval &&
        main(; unknown_storage = :dense, assembly = :cellwise) ≈ testval
    return nothing
end

end # module
