#=

# 121: 1D Poisson Equation with Point Charge
([source code](@__SOURCE_URL__))

This example demonstrates how to handle **point charges (Dirac delta sources)** in the Poisson equation using VoronoiFVM. The problem showcases different methods for implementing singular source terms at specific locations.

## Mathematical Problem

Solve the Poisson equation with a point charge:

```math
-\Delta u = Q \delta(x) \quad \text{in } \Omega = (-1,1)
```

with boundary conditions:
- $u(-1) = 1$ (left boundary)
- $u(1) = 0$ (right boundary)

where $Q$ is the point charge strength located at $x = 0$, and $\delta(x)$ is the Dirac delta function.

## Physical Interpretation

This problem models:
- **Electrostatic potential** $u(x)$ in a 1D conductor with a point charge
- **Heat conduction** with a localized heat source
- **Pressure distribution** in a pipe with a point pressure source

The point charge creates a **discontinuity in the gradient** $\nabla u$ at $x = 0$, while the potential $u$ itself remains continuous.

## Analytical Solution

For this 1D problem, the exact solution is piecewise linear:

```math
u(x) = \begin{cases}
1 + \frac{Q}{2}(x+1) & \text{for } x \in [-1,0] \\
\frac{Q}{2}(1-x) & \text{for } x \in [0,1]
\end{cases}
```

The jump in the derivative at $x = 0$ is: $[\partial_x u]_0 = Q$

## Implementation Methods

The example demonstrates clean parameter management using a data structure:

### Data Structure Approach (Recommended)
Use a mutable struct to encapsulate problem parameters:
```julia
mutable struct ProblemData
    Q::Float64  # Point charge strength
end
```

The point charge is implemented through a unified boundary reaction that accesses the data:
```math
\int_{\partial\omega_0} Q \, ds = Q \quad \text{(at x=0)}
```
with Dirichlet conditions applied directly in the same function.

### Alternative: Direct Boundary Value Assignment
Set the charge directly as a boundary value at the special boundary region using the system's boundary value array.

## Modern VoronoiFVM Features

This example showcases modern VoronoiFVM syntax:
- **Parameter encapsulation**: Using mutable structs for clean parameter management
- **Unified boundary treatment**: Both point sources and boundary conditions in `breaction!`
- **Direct species specification**: Using `species = [1]` in system creation
- **Data passing**: Clean parameter access through the `data` argument

## Grid Design

The computational grid uses **geometric refinement** around $x = 0$ to accurately capture the singular behavior:
- Fine mesh near the point charge ($h_{min}$)
- Coarser mesh away from the singularity ($h_{max}$)
- Special boundary region at $x = 0$ for the point charge

This adaptive mesh design ensures both accuracy and efficiency.

=#

module Example121_PoissonPointCharge1D

using Printf

using VoronoiFVM
using ExtendableGrids
using GridVisualize

## Mutable struct to hold problem parameters
## This allows clean parameter passing and modification during simulation
mutable struct ProblemData
    Q::Float64  ## Point charge strength
end

function main(;
        nref = 0, Plotter = nothing, verbose = false, unknown_storage = :sparse,
        assembly = :edgewise
    )

    ## Create geometrically graded grid in (-1,1) with refinement around x=0
    ## This captures the singular behavior of the point charge accurately
    hmax = 0.2 / 2.0^nref   ## Coarse mesh size away from singularity
    hmin = 0.05 / 2.0^nref  ## Fine mesh size near the point charge

    ## Generate geometric spacing: fine near 0, coarse at boundaries
    X1 = geomspace(-1.0, 0.0, hmax, hmin)  ## Left half: coarse to fine
    X2 = geomspace(0.0, 1.0, hmin, hmax)   ## Right half: fine to coarse
    X = glue(X1, X2)  ## Combine the two halves
    grid = simplexgrid(X)

    ## Configure grid regions for multi-physics setup
    ## Create special boundary region 3 at x=0 for the point charge
    bfacemask!(grid, [0.0], [0.0], 3)

    ## Define material regions (could have different properties)
    cellmask!(grid, [-1.0], [0.0], 1)  ## Left material region
    cellmask!(grid, [0.0], [1.0], 2)   ## Right material region

    ## Initialize problem data structure with zero charge
    problem_data = ProblemData(0.0)

    ## Flux function: implements -∇u (Laplacian discretization)
    ## For Poisson equation: flux = gradient of potential
    function flux!(f, u, edge, data)
        f[1] = u[1, 1] - u[1, 2]  ## Discrete gradient: (u_left - u_right)
        return nothing
    end

    ## Storage function: for time-dependent problems (not used here)
    ## Included for completeness and potential extensions
    function storage!(f, u, node, data)
        f[1] = u[1]  ## Simple storage: ∂u/∂t
        return nothing
    end

    ## Boundary reaction: implements both the point charge and boundary conditions
    ## This unified approach handles all boundary terms in one function
    ## Note: negative sign for point charge due to VoronoiFVM convention (LHS formulation)
    function breaction!(f, u, node, data)
        if node.region == 3  ## Apply point charge at the special boundary region (x=0)
            f[1] = -data.Q  ## Point charge contribution: -Q⋅δ(x), now from data struct
        end
        ## Apply Dirichlet boundary conditions directly in the boundary reaction
        boundary_dirichlet!(f, u, node; region = 1, value = 1.0)  ## u(-1) = 1 (left boundary)
        boundary_dirichlet!(f, u, node; region = 2, value = 0.0)  ## u(1) = 0 (right boundary)
        return nothing
    end

    ## Assemble the physics: Laplacian + point source + boundary conditions
    ## Pass the problem data structure to make parameters accessible to all physics functions
    physics = VoronoiFVM.Physics(;
        flux = flux!,        ## Diffusion/conduction term
        storage = storage!,  ## Time derivative (unused here)
        breaction = breaction!,  ## Point charge source and boundary conditions
        data = problem_data  ## Problem parameters
    )

    ## Create the finite volume system with species 1 enabled
    ## This modern approach directly specifies species in system creation
    sys = VoronoiFVM.System(grid, physics; unknown_storage = :dense, assembly = assembly, species = [1])


    ## Initialize solution array
    U = unknowns(sys)
    U .= 0  ## Start with zero potential everywhere

    ## Configure Newton solver
    control = VoronoiFVM.NewtonControl()
    control.verbose = verbose

    ## Set up visualization
    vis = GridVisualizer(; Plotter = Plotter)

    ## Parameter study: solve for increasing point charge strengths
    ## This demonstrates how the solution changes with source strength
    for q in [0.1, 0.2, 0.4, 0.8, 1.6]
        ## Method 1: Modify the charge in the problem data structure (recommended approach)
        ## This cleanly updates the parameter without touching global variables
        problem_data.Q = q

        ## Solve the linear system (Newton converges in 1 iteration for linear problems)
        U = solve(sys; inival = U, control)

        ## Visualize the solution
        ## Should show piecewise linear profile with slope discontinuity at x=0
        scalarplot!(
            vis, grid, U[1, :]; title = @sprintf("Q=%.2f", q), clear = true,
            show = true
        )
    end

    ## Return sum of solution values for testing purposes
    return sum(U)
end

using Test
function runtests()
    testval = 20.254591679579015
    @test main(; assembly = :edgewise) ≈ testval &&
        main(; assembly = :cellwise) ≈ testval
    return nothing
end
end
