#=

# 125: Terminal Flux Calculation via Test Functions, 1D
([source code](@__SOURCE_URL__))

This example demonstrates **test function-based flux calculation** for computing boundary fluxes in finite volume problems. The method uses auxiliary test functions to convert difficult-to-compute boundary integrals into volume integrals that can be accurately evaluated using the finite volume discretization.

## Problem Setup

Consider a 1D reaction-diffusion system on $\Omega = [0,1]$ with two competing species:

```math
\begin{aligned}
\partial_t u_1 - \nabla \cdot (D_1 \nabla u_1) + R &= 0\\
\partial_t u_2 - \nabla \cdot (D_2 \nabla u_2) - R &= 0
\end{aligned}
```

where the reaction term $R = 10(u_1 - u_2)$ represents conversion between species:
- Species 1 converts to species 2 when $u_1 > u_2$
- Species 2 converts to species 1 when $u_2 > u_1$

## Boundary Conditions

- **Left boundary** ($x=0$): Neumann flux condition for species 1: $\partial_n u_1 = 0.01$
- **Right boundary** ($x=1$): Dirichlet condition for species 2: $u_2 = 0$

## Test Function Method

Computing the exact flux through boundaries in finite volume methods can be challenging. The **test function approach** provides an elegant solution:

### Theory

To compute the flux $\int_{\Gamma} \mathbf{j} \cdot \mathbf{n} \, ds$ through boundary $\Gamma$, we construct a test function $T(x)$ that solves:

```math
-\nabla^2 T = 0 \quad \text{in } \Omega
```

with boundary conditions:
- $T = 1$ on the boundary where we want to measure flux
- $T = 0$ on other boundaries  
- $\partial_n T = 0$ on remaining boundaries

### Implementation

The example constructs two test functions:
- `tf1`: Measures flux from boundary 2 to boundary 1
- `tf2`: Measures flux from boundary 1 to boundary 2

Using Green's theorem, the boundary flux integral becomes:

```math
\int_{\Gamma} \mathbf{j} \cdot \mathbf{n} \, ds = \int_\Omega \nabla T \cdot \mathbf{j} + T(\nabla \cdot \mathbf{j}) \, d\omega
```

### Complete Derivation with Reaction Terms

For our reaction-diffusion system, we have the flux $\mathbf{j} = -D\nabla u$ and the PDE:
```math
\partial_t u - \nabla \cdot \mathbf{j} + R = 0
```

Therefore: $\nabla \cdot \mathbf{j} = \partial_t u + R$

For the stationary case ($\partial_t u = 0$), we get: $\nabla \cdot \mathbf{j} = R$

The complete flux calculation becomes:
```math
\begin{aligned}
\int_{\Gamma} \mathbf{j} \cdot \mathbf{n} \, ds &= \int_{\Gamma} T\mathbf{j} \cdot \mathbf{n} \, ds \quad \text{(since } T=1 \text{ on } \Gamma\text{)}\\
&= \int_{\partial\Omega} T\mathbf{j} \cdot \mathbf{n} \, ds \quad \text{(since } T=0 \text{ on other boundaries)}\\
&= \int_\Omega \nabla \cdot (T \mathbf{j}) \, d\omega \quad \text{(Gauss theorem)}\\
&= \int_\Omega \nabla T \cdot \mathbf{j} \, d\omega + \int_\Omega T \nabla\cdot \mathbf{j} \, d\omega\\
&= \int_\Omega \nabla T \cdot \mathbf{j} \, d\omega + \int_\Omega T \cdot R \, d\omega
\end{aligned}
```

### Finite Volume Approximation

The VoronoiFVM implementation approximates these integrals using the finite volume discretization:

```math
\int_\Omega \nabla T \cdot \mathbf{j} \, d\omega \approx \sum_{k,l} \frac{|\omega_k \cap \omega_l|}{h_{k,l}} g(u_k, u_l) (T_k - T_l)
```

where:
- The sum runs over pairs of neighboring control volumes
- $g(u_k, u_l) = D(u_k - u_l)$ is the finite volume flux
- $|\omega_k \cap \omega_l|/h_{k,l}$ is the interface area divided by distance

The reaction term integral is approximated as:
```math
\int_\Omega T \cdot R \, d\omega \approx \sum_k |\omega_k| T_k R_k
```

This volume integral can be accurately computed using the finite volume discretization.

## Parameter Study

The example varies the diffusion coefficients $D_1 = D_2 = \varepsilon$ from 1.0 to 0.01, demonstrating how:
- **High diffusion** ($\varepsilon = 1.0$): Species profiles are smooth
- **Low diffusion** ($\varepsilon = 0.01$): Sharp gradients form, requiring careful flux calculation

## Physical Interpretation

This models competing transport processes:
- **Species injection**: Species 1 enters at the left boundary
- **Reaction zone**: Conversion between species occurs throughout the domain
- **Species removal**: Species 2 exits at the right boundary

The test function method accurately quantifies the net transport rates between boundaries, which is crucial for:
- Mass balance verification
- Process optimization
- Validation of numerical schemes

For a more comprehensive explanation, see [Example225: Terminal flux calculation via test functions, nD](@ref).

=#

module Example125_TestFunctions1D
using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize

## Mutable struct to hold problem parameters
## This encapsulates diffusion coefficients and other physical parameters
mutable struct ProblemData
    eps::Vector{Float64}  ## Diffusion coefficients for both species
end

function main(; n = 100, Plotter = nothing, verbose = false, unknown_storage = :sparse, assembly = :edgewise)
    ## Create a 1D grid with n intervals on [0,1]
    h = 1 / n
    grid = simplexgrid(collect(0:h:1))

    ## Initialize problem parameters in data structure
    ## Diffusion coefficients for both species (will be varied in the loop)
    problem_data = ProblemData([1.0, 1.0e-1])

    ## Define the physics of the problem
    physics = VoronoiFVM.Physics(
        ## Reaction terms: Species 1 converts to species 2 and vice versa
        ## Rate depends on concentration difference (Le Chatelier principle)
        reaction = function (f, u, node, data)
            f[1] = 10 * (u[1] - u[2])  ## Species 1 decreases when u1 > u2
            f[2] = 10 * (u[2] - u[1])  ## Species 2 increases when u1 > u2
            return nothing
        end, 
        ## Flux terms: Simple diffusion for both species
        ## f[i] = D_i * (u[i,1] - u[i,2]) represents diffusive flux
        flux = function (f, u, edge, data)
            f[1] = data.eps[1] * (u[1, 1] - u[1, 2])  ## Diffusion flux for species 1
            f[2] = data.eps[2] * (u[2, 1] - u[2, 2])  ## Diffusion flux for species 2
            return nothing
        end, 
        ## Storage terms: Simple time derivative terms
        storage = function (f, u, node, data)
            f[1] = u[1]  ## ∂u1/∂t
            f[2] = u[2]  ## ∂u2/∂t
            return nothing
        end,
        data = problem_data  ## Pass problem parameters
    )
    
    ## Create the finite volume system
    sys = VoronoiFVM.System(grid, physics; unknown_storage = unknown_storage, assembly = assembly)

    ## Enable both species in the single region
    enable_species!(sys, 1, [1])  ## Species 1 active in region 1
    enable_species!(sys, 2, [1])  ## Species 2 active in region 1

    ## Set boundary conditions
    boundary_neumann!(sys, 1, 1, 0.01)    ## Constant flux of species 1 at left boundary (x=0)
    boundary_dirichlet!(sys, 2, 2, 0.0)   ## Species 2 concentration = 0 at right boundary (x=1)

    ## Create test function factory for flux calculations
    factory = TestFunctionFactory(sys)
    ## tf1: Test function = 1 at boundary 1, = 0 at boundary 2 (measures flux from right to left)
    tf1 = testfunction(factory, [2], [1])
    ## tf2: Test function = 1 at boundary 2, = 0 at boundary 1 (measures flux from left to right)  
    tf2 = testfunction(factory, [1], [2])

    ## Initialize solution arrays
    inival = unknowns(sys)
    inival[2, :] .= 0.1  ## Initial concentration of species 2
    inival[1, :] .= 0.1  ## Initial concentration of species 1

    ## Configure Newton solver
    control = VoronoiFVM.SolverControl()
    control.verbose = verbose
    control.damp_initial = 0.1  ## Use damping to help convergence
    
    ## Initialize flux integral result
    I1 = 0
    
    ## Set up visualization
    p = GridVisualizer(; Plotter = Plotter, layout = (2, 1))
    
    ## Parameter study: vary diffusion coefficients to test method robustness
    for xeps in [1.0, 0.1, 0.01]
        ## Update diffusion coefficients in the problem data structure
        problem_data.eps = [xeps, xeps]  ## Set both species to same diffusion coefficient
        
        ## Solve the stationary problem
        U = solve(sys; inival, control)
        
        ## Calculate flux integral using test function tf1
        ## This gives the net flux from boundary 2 to boundary 1
        I1 = integrate(sys, tf1, U)
        
        ## Get grid coordinates for visualization
        coord = coordinates(grid)
        
        ## Use current solution as initial guess for next iteration
        inival .= U
        
        ## Visualize the solutions
        scalarplot!(p[1, 1], grid, U[1, :])  ## Plot species 1 concentration
        scalarplot!(p[2, 1], grid, U[2, :])  ## Plot species 2 concentration
        reveal(p)
        
        ## Store a test value (concentration at grid point 5)
        u5 = U[5]
    end
    
    ## Return the flux integral for species 1
    ## This should equal the input flux (0.01) for mass balance
    return I1[1]
end

using Test
function runtests()
    testval = 0.01
    @test main(; unknown_storage = :sparse, assembly = :edgewise) ≈ testval
    @test main(; unknown_storage = :sparse, assembly = :cellwise) ≈ testval
    @test main(; unknown_storage = :dense, assembly = :cellwise) ≈ testval
    @test main(; unknown_storage = :dense, assembly = :edgewise) ≈ testval
    return nothing
end
end
