#=

# 103: 1D Convection-diffusion equation
([source code](@__SOURCE_URL__))

Solve the equation

```math
\partial_t u -\nabla ( D \nabla u - v u) = 0
```
in $\Omega=(0,1)$ with homogeneous Neumann boundary condition
at $x=0$ and outflow boundary condition at $x=1$.
=#

module Example103_ConvectionDiffusion1D
using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize

## Bernoulli function used in the exponential fitting discretization
function bernoulli(x)
    if abs(x) < nextfloat(eps(typeof(x)))
        return 1
    end
    return x / (exp(x) - 1)
end

function exponential_flux!(f, u, edge, data)
    vh = project(edge, data.v)
    Bplus = data.D * bernoulli(vh / data.D)
    Bminus = data.D * bernoulli(-vh / data.D)
    f[1] = Bminus * u[1, 1] - Bplus * u[1, 2]
    return nothing
end

function outflow!(f, u, node, data)
    if node.region == 2
        f[1] = data.v[1] * u[1]
    end
    return nothing
end

function main(; n = 10, Plotter = nothing, D = 0.01, v = 1.0, tend = 100)

    ## Create a one-dimensional discretization
    h = 1.0 / n
    grid = simplexgrid(0:h:1)

    data = (v = [v], D = D)

    sys = VoronoiFVM.System(
        grid,
        VoronoiFVM.Physics(;
            flux = exponential_flux!, data = data,
            breaction = outflow!
        )
    )

    ## Add species 1 to region 1
    enable_species!(sys, 1, [1])

    ## Set boundary conditions
    boundary_neumann!(sys, 1, 1, 0.0)

    ## Create a solution array
    inival = unknowns(sys)
    inival[1, :] .= map(x -> 1 - 2x, grid)

    ## Transient solution of the problem
    control = VoronoiFVM.NewtonControl()
    control.Δt = 0.01 * h
    control.Δt_min = 0.01 * h
    control.Δt_max = 0.1 * tend
    tsol = solve(sys; inival, times = [0, tend], control)

    vis = GridVisualizer(; Plotter = Plotter)
    for i in 1:length(tsol.t)
        scalarplot!(
            vis[1, 1], grid, tsol[1, :, i]; flimits = (0, 1),
            title = "t=$(tsol.t[i])", show = true
        )
        sleep(0.01)
    end
    return tsol
end

using Test
function runtests()
    tsol = main()
    @test maximum(tsol) <= 1.0 && maximum(tsol.u[end]) < 1.0e-20
    return nothing
end

end
