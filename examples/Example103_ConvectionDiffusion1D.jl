#=

# 103: 1D Transient convection-diffusion equation
([source code](@__SOURCE_URL__))

Solve the time-dependent convection-diffusion equation

```math
\partial_t u -\nabla ( D \nabla u - v u) = 0
```
in $\Omega=(0,1)$ with homogeneous Neumann boundary condition
at $x=0$ and outflow boundary condition at $x=1$.

This is the time-dependent version of the convection-diffusion problem from Example102.
The equation models the evolution of a scalar quantity $u$ (e.g., concentration, temperature)
under the combined effects of diffusion (with coefficient $D$) and advection (with velocity $v$).

## Physical Interpretation

- **Diffusion term**: $D \nabla u$ represents spreading due to random motion
- **Convection term**: $v u$ represents transport by bulk motion of the medium
- **Time evolution**: $\partial_t u$ describes how the quantity changes over time

## Boundary Conditions

- **Left boundary** ($x=0$): Homogeneous Neumann condition $\partial_n u = 0$ (no flux)
- **Right boundary** ($x=1$): Outflow condition $D\partial_n u = 0$, equivalent to $  ( D \nabla u - v u)\cdot \vec n =  v u \cdot \vec n $
allowing material transported by convection towards the boundary  to leave the domain.

## Discretization

The spatial discretization uses the **exponential fitting scheme** (Scharfetter-Gummel method) which maintains monotonicity properties and provides accurate solutions even for convection-dominated transport (high Peclet numbers). This method uses the Bernoulli function:

```math
B(x) = \frac{x}{e^x - 1}
```

to construct fluxes that exactly solve the local two-point boundary value problem on each edge.

## Initial Condition

The simulation starts with a linear profile $u(x,0) = 1 - 2x$, which evolves under the
combined effects of diffusion and convection until it reaches a steady state or is
advected out of the domain.

=#

module Example103_ConvectionDiffusion1D
using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize

## Mutable struct to hold problem parameters
## This encapsulates all physical and numerical parameters
mutable struct ProblemData
    D::Float64      ## Diffusion coefficient
    v::Vector{Float64}  ## Velocity vector
end

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

    ## Initialize problem parameters in data structure
    problem_data = ProblemData(D, [v])

    sys = VoronoiFVM.System(
        grid,
        VoronoiFVM.Physics(;
            flux = exponential_flux!, 
            breaction = outflow!,
            data = problem_data
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
    control = VoronoiFVM.SolverControl()
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
