#=

# 160: Bipolar drift-diffusion with different definition of current
([source code](@__SOURCE_URL__))

The problem consists of a Poisson equation for the electrostatic potential $\psi$:

```math
-\nabla \lambda^2 \nabla \psi = (p-n + C)
```
and drift-diffusion equations of the transport of electrons and holes:

```math
z_n \partial_t n  + \nabla\cdot j_n = z_n R(n, p),\\
z_p \partial_t p  + \nabla\cdot j_p = z_p R(n, p).
```

In particular, this example explores different definitions of the carrier currents and also includes the displacement flux.
=#

module Example161_BipolarDriftDiffusionCurrent

using VoronoiFVM
using ExtendableGrids
using GridVisualize

function main(;
        n = 20, # number of nodes
        Plotter = nothing,
        plotting = false,
        verbose = false
    )

    ################################################################################
    #### grid
    ################################################################################
    h1 = 1.0; h2 = 4.0; h3 = 2.0
    h_total = h1 + h2 + h3

    # region numbers
    region1 = 1
    region2 = 2
    region3 = 3
    regions = [region1, region2, region3]

    # boundary region numbers
    bregion1 = 1
    bregion2 = 2

    # 558 nodes; spacing ≈ 0.0126;
    coord1 = collect(range(0.0; stop = h1, length = n))
    coord2 = collect(range(h1; stop = h1 + h2, length = 4 * n))
    coord3 = collect(range(h1 + h2; stop = h_total, length = 2 * n))
    coord = glue(coord1, coord2)
    coord = glue(coord, coord3)

    grid = simplexgrid(coord)

    cellmask!(grid, [0.0], [h1], region1)
    cellmask!(grid, [h1], [h1 + h2], region2)
    cellmask!(grid, [h1 + h2], [h_total], region3)

    # specify outer regions
    bfacemask!(grid, [0.0], [0.0], bregion1)
    bfacemask!(grid, [h_total], [h_total], bregion2)

    ################################################################################
    #########  system
    ################################################################################

    tPrecond = 5.0
    tExt = 10.0
    tRamp = 1.0e-5
    tEnd = tPrecond + tRamp + tExt
    Vprecond = 2.5
    VExt = 0.0

    ## Define scan protocol function
    function scanProtocol(t)
        if 0.0 <= t && t <= tPrecond
            biasVal = Vprecond
        elseif tPrecond <= t  && t <= tPrecond + tRamp # ramping down
            biasVal = Vprecond - (Vprecond - VExt) / tRamp * (t - tPrecond)
        elseif tPrecond + tRamp <= t  && t <= tEnd
            biasVal = VExt
        else
            biasVal = 0.0
        end

        return biasVal
    end

    Cn = 100        # n-doped
    Cp = 100        # p-doped
    Ca = 0.0        # intrinsic
    En = 2.0        # conduction band edge
    Ep = 0.0        # valence band edge
    μn = 1.0e-1     # electron mobility
    μp = 1.0e-1     # hole mobility
    lambda = 1.0 # Debye length
    DirichletVal = 0.0

    sys = VoronoiFVM.System(grid; unknown_storage = :sparse)
    iphin = 1; zn = -1; enable_species!(sys, iphin, regions)
    iphip = 2; zp = 1;  enable_species!(sys, iphip, regions)
    ipsi = 3; enable_species!(sys, ipsi, regions)

    function reaction!(f, u, node, data)

        nn = exp(zn * (u[iphin] - u[ipsi] + En))
        np = exp(zp * (u[iphip] - u[ipsi] + Ep))

        if node.region == region1
            C = Cn
        elseif node.region == region2
            C = Ca
        elseif node.region == region3
            C = -Cp
        end

        f[ipsi] = -(C + zn * nn + zp * np)
        ########################
        r0 = 1.0

        recomb = (r0 + 1 / (nn + np)) * (nn * np * (1.0 - exp(u[iphin] - u[iphip])))

        f[iphin] = zn * recomb
        f[iphip] = zp * recomb

        return nothing
    end

    function flux!(f, u, node, data)

        f[ipsi] = - lambda^2 * (u[ipsi, 2] - u[ipsi, 1])

        ########################
        bp, bm = fbernoulli_pm(-(u[ipsi, 2] - u[ipsi, 1]))

        nn1 = exp(zn * (u[iphin, 1] - u[ipsi, 1] + En))
        np1 = exp(zp * (u[iphip, 1] - u[ipsi, 1] + Ep))

        nn2 = exp(zn * (u[iphin, 2] - u[ipsi, 2] + En))
        np2 = exp(zp * (u[iphip, 2] - u[ipsi, 2] + Ep))

        f[iphin] = - zn * μn * (bm * nn2 - bp * nn1)
        f[iphip] = - zp * μp * (bp * np2 - bm * np1)
        return nothing
    end

    function breaction!(f, u, bnode, data)

        boundary_dirichlet!(f, u, bnode, iphin, bregion1, 0.0)
        boundary_dirichlet!(f, u, bnode, iphip, bregion1, 0.0)
        boundary_dirichlet!(f, u, bnode, ipsi, bregion1, 0.5 * (En + Ep) + asinh(Cn / (2 * sqrt(exp(- (En - Ep))))))

        boundary_dirichlet!(f, u, bnode, iphin, bregion2, DirichletVal + scanProtocol(bnode.time))
        boundary_dirichlet!(f, u, bnode, iphip, bregion2, DirichletVal + scanProtocol(bnode.time))
        boundary_dirichlet!(f, u, bnode, ipsi, bregion2, 0.5 * (En + Ep) + asinh(-Cp / (2 * sqrt(exp(- (En - Ep))))) + DirichletVal + scanProtocol(bnode.time))
        return nothing
    end

    function storage!(f, u, node, data)
        nn = exp(zn * (u[iphin] - u[ipsi] + En))
        np = exp(zp * (u[iphip] - u[ipsi] + Ep))

        f[iphin] = zn * nn
        f[iphip] = zp * np
        return nothing
    end

    physics!(
        sys,
        VoronoiFVM.Physics(;
            flux = flux!,
            reaction = reaction!,
            breaction = breaction!,
            storage = storage!
        )
    )

    ################################################################################
    #########  time loop
    ################################################################################

    control = SolverControl()
    control.Δu_opt = Inf
    control.max_round = 3
    control.damp_initial = 0.9
    control.damp_growth = 1.61 # >= 1
    control.verbose = verbose
    control.abstol = 1.0e-5
    control.reltol = 1.0e-5
    control.tol_round = 1.0e-5

    ## Create a solution array
    inival2 = unknowns(sys)

    inival2[iphin, :] = 0.0 .+ DirichletVal ./ h_total .* coord
    inival2[iphip, :] = 0.0 .+ DirichletVal ./ h_total .* coord
    inival2[ipsi, :] = asinh(Cn / 2) .+ ((asinh(-Cp / 2) + DirichletVal) - asinh(Cn / 2)) ./ h_total .* coord

    solPrecond = solve(sys, inival = inival2, times = (0.0, tPrecond), control = control)
    control.Δt_min = 1.0e-3
    control.Δt = 1.0e-3
    solRamp = solve(sys, inival = solPrecond.u[end], times = (tPrecond, tPrecond + tRamp), control = control)
    #####
    control.Δt_min = 1.0e-3
    control.Δt = 1.0e-3
    solExt = solve(sys, inival = solRamp.u[end], times = (tPrecond + tRamp, tEnd), control = control)

    ################################################################################
    ## Current calculation
    ################################################################################

    # for saving Current data
    I1 = zeros(0)
    In1 = zeros(0); Ip1 = zeros(0); Iψ1 = zeros(0)

    I2 = zeros(0)
    In2 = zeros(0); Ip2 = zeros(0); Iψ2 = zeros(0)

    I3 = zeros(0)

    tvalues = solExt.t
    number_tsteps = length(tvalues)

    factory = TestFunctionFactory(sys)
    tf = testfunction(factory, [bregion1], [bregion2])

    for istep in 2:number_tsteps

        Δt = tvalues[istep] - tvalues[istep - 1] # Time step size
        inival = solExt.u[istep - 1]
        solution = solExt.u[istep]

        ## Variant 1
        II = integrate(sys, tf, solution, inival, Δt)
        IIDisp = VoronoiFVM.integrate_displacement(sys, tf, solution, inival, Δt)

        push!(In1, II[iphin]); push!(Ip1, II[iphip]); push!(Iψ1, IIDisp[ipsi])
        push!(I1, In1[istep - 1] + Ip1[istep - 1] + Iψ1[istep - 1])

        push!(I3, In1[istep - 1] + Ip1[istep - 1])

        ## Variant 2
        IIEdge = VoronoiFVM.integrate_edgebatch(sys, tf, solution, inival, Δt)
        IIDispEdge = VoronoiFVM.integrate_displacement_edgebatch(sys, tf, solution, inival, Δt)

        push!(In2, IIEdge[iphin]); push!(Ip2, IIEdge[iphip]); push!(Iψ2, IIDispEdge[ipsi])
        push!(I2, In2[istep - 1] + Ip2[istep - 1] + Iψ2[istep - 1])

    end

    ################################################################################
    #########  Plotting
    ################################################################################

    tvalues_shift = tvalues .- (tPrecond + tRamp)

    if plotting

        vis1 = GridVisualizer(; layout = (3, 1), xlabel = "space", ylabel = "potential", legend = :lt, Plotter = Plotter, fignumber = 1)

        sol1 = solExt.u[1]
        sol2 = solExt.u[end]

        scalarplot!(vis1[1, 1], coord, sol1[iphin, :]; color = :darkgreen, label = "phin (start)", linewidth = 5, clear = false)
        scalarplot!(vis1[1, 1], coord, sol2[iphin, :]; color = :lightgreen, label = "phin (end)", linewidth = 5, clear = false)
        ####
        scalarplot!(vis1[2, 1], coord, sol1[iphip, :]; color = :darkred, label = "phip (start)", linewidth = 5, clear = false)
        scalarplot!(vis1[2, 1], coord, sol2[iphip, :]; color = :red, label = "phip (end)", linewidth = 5, clear = false)
        ####
        scalarplot!(vis1[3, 1], coord, sol1[ipsi, :]; color = :darkblue, label = "psi (start)", linewidth = 5, clear = false)
        scalarplot!(vis1[3, 1], coord, sol2[ipsi, :]; color = :lightblue, label = "psi (end)", linewidth = 5, clear = false)

        ###################
        vis2 = GridVisualizer(; layout = (3, 1), xlabel = "time", ylabel = "current", legend = :lt, xscale = :log, yscale = :log, Plotter = Plotter, fignumber = 2)

        scalarplot!(vis2[1, 1], tvalues_shift[2:end], abs.(In1); color = :darkgreen, label = "Jn", linewidth = 5, clear = false)
        scalarplot!(vis2[1, 1], tvalues_shift[2:end], abs.(Ip1); color = :darkred, label = "Jp", clear = false)
        scalarplot!(vis2[1, 1], tvalues_shift[2:end], abs.(Iψ1); color = :darkblue, label = "Jdisp", clear = false)
        scalarplot!(vis2[1, 1], tvalues_shift[2:end], abs.(I1); color = :black, linestyle = :dash, label = "Jtot", clear = false)
        ####
        scalarplot!(vis2[2, 1], tvalues_shift[2:end], abs.(In2); color = :darkgreen, label = "Jn", linewidth = 5, clear = false)
        scalarplot!(vis2[2, 1], tvalues_shift[2:end], abs.(Ip2); color = :darkred, label = "Jp", clear = false)
        scalarplot!(vis2[2, 1], tvalues_shift[2:end], abs.(Iψ2); color = :darkblue, label = "Jdisp", clear = false)
        scalarplot!(vis2[2, 1], tvalues_shift[2:end], abs.(I2); color = :black, linestyle = :dash, label = "Jtot", clear = false)
        #####
        scalarplot!(vis2[3, 1], tvalues_shift[2:end], abs.(I1); color = :red, linestyle = :solid, linewidth = 5, label = "Jtot (#1)", clear = false)
        scalarplot!(vis2[3, 1], tvalues_shift[2:end], abs.(I2); color = :green, linestyle = :dash, label = "Jtot (#2)", clear = false)
        scalarplot!(vis2[3, 1], tvalues_shift[2:end], abs.(I3); linestyle = :solid, color = :gray, label = "Jtot without disp", clear = false)

    end

    ###################################
    if sum(abs.(I1 - I2)) < 1.0e-11
        return true
    else
        return false
    end
end # main

using Test
function runtests()
    @test main() == true
    return nothing
end

end # module
