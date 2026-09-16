#=

# 421: Current Calculation for AbstractQuantities
 ([source code](@__SOURCE_URL__))

Test current calculation for jumping species. Here, we have three cases:
    a. Problem initialized as usual
    b. Problem initialized with Continuousquantity
    c. Problem initialized with Discontinuousquantity with adjusted reaction rate
We see that the resulting current coincides for all three cases when adjusting the
reaction rate.

=#
module Example421_AbstractQuantities_TestFunctions

using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize
using LinearAlgebra

mutable struct Data
    rate::Float64 # rate which is within DiscontinuousQuantities
    Data() = new()
end

function main(; N = 3, Plotter = nothing, unknown_storage = :sparse, assembly = :edgewise)
    XX = collect(0:0.1:1)
    xcoord = XX
    for i in 1:(N - 1)
        xcoord = glue(xcoord, XX .+ i)
    end
    grid = simplexgrid(xcoord)
    for i in 1:N
        cellmask!(grid, [i - 1], [i], i)
    end
    for i in 1:(N - 1)
        bfacemask!(grid, [i], [i], i + 2)
    end

    sysQ = VoronoiFVM.System(grid; unknown_storage = unknown_storage)
    cspec = ContinuousQuantity(sysQ, 1:N; id = 1)                    # continuous quantity
    dspec = DiscontinuousQuantity(sysQ, 1:N; id = 2)                 # discontinuous quantity

    data = Data()
    rate = 0.0
    data.rate = rate

    function fluxQ(f, u, edge, data) # For both quantities, we define simple diffusion fluxes
        f[dspec] = u[dspec, 1] - u[dspec, 2]
        f[cspec] = u[cspec, 1] - u[cspec, 2]
        return nothing
    end

    function breactionQ(f, u, bnode, data)
        ## Define a thin layer interface condition for `dspec`.
        if bnode.region > 2
            react = (u[dspec, 1] - u[dspec, 2]) / data.rate
            f[dspec, 1] = react
            f[dspec, 2] = -react
        end
        return nothing
    end

    physics!(
        sysQ, VoronoiFVM.Physics(;
            data = data,
            flux = fluxQ,
            breaction = breactionQ
        )
    )

    ##########################################################
    icc = 1 # for system without AbstractQuantities

    function flux!(f, u, edge, data) # analogous as for other system
        f[icc] = u[icc, 1] - u[icc, 2]
        return nothing
    end

    # other system to which we compare current calculation
    sys = VoronoiFVM.System(
        grid; flux = flux!, species = icc,
        unknown_storage = unknown_storage
    )

    ## Set left boundary conditions
    boundary_dirichlet!(sysQ, dspec, 1, 0.0)
    boundary_dirichlet!(sysQ, cspec, 1, 0.0)
    boundary_dirichlet!(sys, icc, 1, 0.0)

    subgrids = VoronoiFVM.subgrids(dspec, sysQ)

    # solve
    UQ = unknowns(sysQ)
    U = unknowns(sys)
    UQ .= 0.0
    U .= 0.0
    biasval = range(0; stop = 2.0, length = 5)

    Icspec = zeros(length(biasval))
    Idspec = zeros(length(biasval))
    Iicc = zeros(length(biasval))

    for data.rate in [1.0e2, 1.0e0, 1.0e-2, 1.0e-4, 1.0e-6]
        count = 1
        for Δu in biasval
            # first problem
            boundary_dirichlet!(sysQ, dspec, 2, Δu)
            boundary_dirichlet!(sysQ, cspec, 2, Δu)

            UQ = solve(sysQ; inival = UQ)

            ## get current
            factoryQ = TestFunctionFactory(sysQ)
            tfQ = testfunction(factoryQ, [1], [2])
            IQ = integrate(sysQ, tfQ, UQ)

            val = 0.0
            for ii in dspec.regionspec # current is calculated regionwise
                val = val + IQ[ii]
            end
            Icspec[count] = IQ[cspec]
            Idspec[count] = val

            # second problem
            boundary_dirichlet!(sys, icc, 2, Δu)

            U = solve(sys; inival = U)

            factory = TestFunctionFactory(sys)
            tf = testfunction(factory, [1], [2])
            I = integrate(sys, tf, U)

            Iicc[count] = I[icc]

            count = count + 1
        end # bias loop

        # plot
        dvws = views(UQ, dspec, subgrids, sysQ)
        cvws = views(UQ, cspec, subgrids, sysQ)

        vis = GridVisualizer(; layout = (2, 1), resolution = (600, 300), Plotter = Plotter)

        for i in eachindex(dvws)
            scalarplot!(
                vis[1, 1], subgrids[i], dvws[i]; flimits = (-0.5, 1.5),
                title = @sprintf("Solution with rate=%.2f", data.rate),
                label = "discont quantity", clear = false, color = :red
            )
            scalarplot!(
                vis[1, 1], subgrids[i], cvws[i]; label = "cont quantity",
                clear = false, color = :green
            )
        end
        scalarplot!(
            vis[1, 1], grid, U[icc, :]; label = "without quantity", clear = false,
            linestyle = :dot, color = :blue
        )

        scalarplot!(
            vis[2, 1], biasval, Idspec; clear = false,
            title = @sprintf("IV with rate=%.2f", data.rate),
            label = "discont quantity", color = :red
        )
        scalarplot!(
            vis[2, 1], biasval, Icspec; clear = false, title = "Current",
            label = "cont quantity", color = :green
        )
        scalarplot!(
            vis[2, 1], biasval, Iicc; clear = false, label = "discont quantity",
            linestyle = :dot, color = :blue, show = true
        )

        reveal(vis)
        sleep(0.2)
    end # rate loop

    errorIV = norm(Idspec - Icspec, 2)

    return errorIV
end

using Test
function runtests()
    testval = 6.085802139465579e-7
    @test main(; unknown_storage = :sparse, assembly = :edgewise) ≈ testval &&
        main(; unknown_storage = :dense, assembly = :edgewise) ≈ testval &&
        main(; unknown_storage = :sparse, assembly = :cellwise) ≈ testval &&
        main(; unknown_storage = :dense, assembly = :cellwise) ≈ testval
    return nothing
end

end
