#=

# 102: 1D Stationary convection-diffusion equation
([source code](@__SOURCE_URL__))

Solve the equation

```math
-\nabla ( D \nabla u - v u) = 0
```
in $\Omega=(0,1)$ with boundary condition $u(0)=0$ and $u(1)=1$.
$v$ could be e.g. the velocity
of a moving medium or the gradient of an electric field.

This is a convection dominant second order boundary value problem which obeys
a local and a global maximum principle:
the solution which is bounded by the values at the boundary and has no local extrema in the
interior.
If $v$ is large compared to $D$, a boundary layer is observed.

The maximum principle of the solution can only be guaranteed it the discretization is
performed accordingly: the flux function must monotonically increase in the first argument
and monotonically decrease in the second argument.

The example describes three possible ways to define the flux function and demonstrates
the impact on the qualitative properties of the solution.

=#

module Example102_StationaryConvectionDiffusion1D
using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize

## Central difference flux. The velocity term is discretized using the
## average of the solution in the endpoints of the grid. If the local Peclet
## number v*h/D>1, the monotonicity property is lost.  Grid refinement
## can fix this situation by decreasing $h$.

function central_flux!(f, u, edge, data)
    f_diff = data.D * (u[1, 1] - u[1, 2])
    vh = project(edge, data.v)
    f[1] = f_diff + vh * (u[1, 1] + u[1, 2]) / 2
    return nothing
end

## The simple upwind flux corrects the monotonicity properties essentially
## via brute force and loses one order of convergence for small $h$ compared
## to the central flux.

function upwind_flux!(f, u, edge, data)
    fdiff = data.D * (u[1] - u[1, 2])
    vh = project(edge, data.v)
    if vh > 0
        f[1] = fdiff + vh * u[1, 1]
    else
        f[1] = fdiff + vh * u[1, 2]
    end
    return nothing
end

## The exponential fitting flux has the proper monotonicity properties and
## kind of interpolates in a clever way between central
## and upwind flux. It can be derived by solving the two-point boundary value problem
## at the grid interval analytically.

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

function calculate(grid, data, flux, verbose)
    sys = VoronoiFVM.System(grid, VoronoiFVM.Physics(; flux = flux, data = data))

    ## Add species 1 to region 1
    enable_species!(sys, 1, [1])

    ## Set boundary conditions
    boundary_dirichlet!(sys, 1, 1, 0.0)
    boundary_dirichlet!(sys, 1, 2, 1.0)

    ## Create a solution array
    inival = unknowns(sys; inival = 0.5)
    solution = unknowns(sys)

    ## Create solver control info
    control = VoronoiFVM.NewtonControl()
    control.verbose = verbose

    ## Stationary solution of the problem
    solution = solve(sys; inival, verbose)
    return solution
end

function main(; n = 10, Plotter = nothing, verbose = false, D = 0.01, v = 1.0)

    ## Create a one-dimensional discretization
    h = 1.0 / convert(Float64, n)
    grid = simplexgrid(collect(0:h:1))

    data = (v = [v], D = D)

    # Calculate three stationary solutions with different ways to calculate flux
    solution_exponential = calculate(grid, data, exponential_flux!, verbose)
    solution_upwind = calculate(grid, data, upwind_flux!, verbose)
    solution_central = calculate(grid, data, central_flux!, verbose)

    # Visualize solutions using GridVisualize.jl
    p = GridVisualizer(; Plotter = Plotter, layout = (3, 1))
    scalarplot!(p[1, 1], grid, solution_exponential[1, :]; title = "exponential")
    scalarplot!(p[2, 1], grid, solution_upwind[1, :]; title = "upwind")
    scalarplot!(p[3, 1], grid, solution_central[1, :]; title = "centered", show = true)

    ## Return test value
    return sum(solution_exponential) + sum(solution_upwind) + sum(solution_central)
end

using Test
function runtests()
    testval = 2.523569744561089
    @test main() ≈ testval
    return nothing
end

end
