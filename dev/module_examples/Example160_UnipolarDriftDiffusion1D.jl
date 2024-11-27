#=

# 160: Unipolar degenerate drift-diffusion
([source code](@__SOURCE_URL__))

See: C. Cancès, C. Chainais-Hillairet, J. Fuhrmann, and B. Gaudeul, "A numerical-analysis-focused comparison of several finite volume schemes for a unipolar degenerate drift-diffusion model" IMA Journal of Numerical Analysis, vol. 41, no. 1, pp. 271–314, 2021.  

Available from [https://doi.org/10.1093/imanum/draa002](https://doi.org/10.1093/imanum/draa002),
the preprint is on [arxiv1907.11126](https://arxiv.org/abs/1907.11126).

The problem consists of a Poisson equation for the electrostatic potential $\phi$:

```math
-\nabla \varepsilon \nabla \phi = z(2c-1)
```
and a degenerate drift-diffusion equation of the transport of a charged species $c$:

```math
\partial_t u  - \nabla\cdot \left(\nabla c  + c \nabla (\phi - \log (1-c) )\right)
```

In particular, the paper, among others, investigates the "sedan" flux discretization which is able to handle the degeneracy coming from the $\log (1-c)$ term. The earliest reference to this scheme we found in the source code of the [SEDAN III](http://www-tcad.stanford.edu/tcad/programs/sedan3.html) semiconductor device simulator.
=#

module Example160_UnipolarDriftDiffusion1D

using Printf

using VoronoiFVM
using ExtendableGrids
using GridVisualize
using LinearSolve

mutable struct Data
    eps::Float64
    z::Float64
    ic::Int32
    iphi::Int32
    V::Float64
    Data() = new()
end

function classflux!(f, u, edge, data)
    ic = data.ic
    iphi = data.iphi
    f[iphi] = data.eps * (u[iphi, 1] - u[iphi, 2])
    bp, bm = fbernoulli_pm(u[iphi, 1] - u[iphi, 2])
    f[ic] = bm * u[ic, 1] - bp * u[ic, 2]
    return nothing
end

function storage!(f, u, node, data)
    ic = data.ic
    iphi = data.iphi
    f[iphi] = 0
    f[ic] = u[ic]
    return nothing
end

function reaction!(f, u, node, data)
    ic = data.ic
    iphi = data.iphi
    f[iphi] = data.z * (1 - 2 * u[ic])
    f[ic] = 0
    return nothing
end
const eps_reg = 1.0e-10
function sedanflux!(f, u, edge, data)
    ic = data.ic
    iphi = data.iphi
    f[iphi] = data.eps * (u[iphi, 1] - u[iphi, 2])
    mu1 = -log1p(max(-1 + eps_reg, -u[ic, 1]))
    mu2 = -log1p(max(-1 + eps_reg, -u[ic, 2]))
    bp, bm = fbernoulli_pm(data.z * 2 * (u[iphi, 1] - u[iphi, 2]) + (mu1 - mu2))
    f[ic] = bm * u[ic, 1] - bp * u[ic, 2]
    return nothing
end

function bcondition!(f, u, bnode, data)
    V = ramp(bnode.time; dt = (0, 1.0e-2), du = (0, data.V))
    boundary_dirichlet!(f, u, bnode; species = data.iphi, region = 1, value = V)
    boundary_dirichlet!(f, u, bnode; species = data.iphi, region = 2, value = 0)
    boundary_dirichlet!(f, u, bnode; species = data.ic, region = 2, value = 0.5)
    return nothing
end

function main(;
        n = 20,
        Plotter = nothing,
        dlcap = false,
        verbose = false,
        phimax = 1,
        dphi = 1.0e-1,
        unknown_storage = :sparse,
        assembly = :edgewise,
    )
    h = 1.0 / convert(Float64, n)
    grid = simplexgrid(collect(0:h:1))

    data = Data()
    data.eps = 1.0e-3
    data.z = -1
    data.iphi = 1
    data.ic = 2
    data.V = 5
    ic = data.ic
    iphi = data.iphi

    physics = VoronoiFVM.Physics(;
        data = data,
        flux = sedanflux!,
        reaction = reaction!,
        breaction = bcondition!,
        storage = storage!,
    )

    sys = VoronoiFVM.System(
        grid,
        physics;
        unknown_storage = unknown_storage,
        species = [1, 2],
        assembly = assembly,
    )

    inival = unknowns(sys)
    @views inival[iphi, :] .= 0
    @views inival[ic, :] .= 0.5

    if !dlcap
        ## Create solver control info for constant time step size
        tstep = 1.0e-5
        control = VoronoiFVM.NewtonControl()
        control.verbose = false
        control.Δt_min = tstep
        control.Δt = tstep
        control.Δt_grow = 1.1
        control.Δt_max = 0.1
        control.Δu_opt = 0.1
        control.damp_initial = 0.5

        tsol = solve(
            sys;
            method_linear = UMFPACKFactorization(),
            inival,
            times = [0.0, 10],
            control = control,
        )

        vis = GridVisualizer(; Plotter = Plotter, layout = (1, 1), fast = true)
        for log10t in -4:0.025:0
            time = 10^(log10t)
            sol = tsol(time)
            scalarplot!(
                vis[1, 1],
                grid,
                sol[iphi, :];
                label = "ϕ",
                title = @sprintf("time=%.3g", time),
                flimits = (0, 5),
                color = :green,
            )
            scalarplot!(
                vis[1, 1],
                grid,
                sol[ic, :];
                label = "c",
                flimits = (0, 5),
                clear = false,
                color = :red,
            )
            reveal(vis)
        end
        return sum(tsol.u[end])

    else  # Calculate double layer capacitance
        U = unknowns(sys)
        control = VoronoiFVM.NewtonControl()
        control.damp_initial = 1.0e-5
        delta = 1.0e-4
        @views inival[iphi, :] .= 0
        @views inival[ic, :] .= 0.5
        sys.boundary_values[iphi, 1] = 0

        delta = 1.0e-4
        vplus = zeros(0)
        cdlplus = zeros(0)
        vminus = zeros(0)
        cdlminus = zeros(0)
        cdl = 0.0
        vis = GridVisualizer(; Plotter = Plotter, layout = (2, 1), fast = true)
        for dir in [1, -1]
            phi = 0.0
            U .= inival
            while phi < phimax
                data.V = dir * phi
                U = solve(sys; inival = U, control, time = 1.0)
                Q = integrate(sys, physics.reaction, U)
                data.V = dir * phi + delta
                U = solve(sys; inival = U, control, time = 1.0)
                Qdelta = integrate(sys, physics.reaction, U)
                cdl = (Qdelta[iphi] - Q[iphi]) / delta

                if Plotter != nothing
                    scalarplot!(
                        vis[1, 1],
                        grid,
                        U[iphi, :];
                        label = "ϕ",
                        title = @sprintf("Δϕ=%.3g", phi),
                        flimits = (-5, 5),
                        clear = true,
                        color = :green,
                    )
                    scalarplot!(
                        vis[1, 1],
                        grid,
                        U[ic, :];
                        label = "c",
                        flimits = (0, 5),
                        clear = false,
                        color = :red,
                    )
                end
                if dir == 1
                    push!(vplus, dir * phi)
                    push!(cdlplus, cdl)
                else
                    push!(vminus, dir * phi)
                    push!(cdlminus, cdl)
                end

                if Plotter != nothing
                    scalarplot!(vis[2, 1], [0, 1.0e-1], [0, 0.05]; color = :white, clear = true)
                end
                v = vcat(reverse(vminus), vplus)
                c = vcat(reverse(cdlminus), cdlplus)
                if length(v) >= 2
                    scalarplot!(
                        vis[2, 1],
                        v,
                        c;
                        color = :green,
                        clear = false,
                        title = "C_dl",
                    )
                end

                phi += dphi
                reveal(vis)
            end
        end

        return cdl
    end
end

using Test
function runtests()


    evolval = 18.721369939565655
    dlcapval = 0.025657355479449806
    rtol = 1.0e-5
    @test isapprox(
        main(; unknown_storage = :sparse, dlcap = false, assembly = :edgewise),
        evolval;
        rtol = rtol,
    )
    @test isapprox(
        main(; unknown_storage = :sparse, dlcap = true, assembly = :edgewise),
        dlcapval;
        rtol = rtol,
    )
    @test isapprox(
        main(; unknown_storage = :dense, dlcap = false, assembly = :edgewise),
        evolval;
        rtol = rtol,
    )
    @test isapprox(
        main(; unknown_storage = :dense, dlcap = true, assembly = :edgewise),
        dlcapval;
        rtol = rtol,
    )
    @test isapprox(
        main(; unknown_storage = :sparse, dlcap = false, assembly = :cellwise),
        evolval;
        rtol = rtol,
    )
    @test isapprox(
        main(; unknown_storage = :sparse, dlcap = true, assembly = :cellwise),
        dlcapval;
        rtol = rtol,
    )
    @test isapprox(
        main(; unknown_storage = :dense, dlcap = false, assembly = :cellwise),
        evolval;
        rtol = rtol,
    )
    @test isapprox(
        main(; unknown_storage = :dense, dlcap = true, assembly = :cellwise),
        dlcapval;
        rtol = rtol,
    )
    return nothing
end
end
