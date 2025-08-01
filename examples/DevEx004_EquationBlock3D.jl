#=

# 004: 3D Nonlinear Reaction-Diffusion System with Coupled Species
([source code](@__SOURCE_URL__))

This example implements a **3D nonlinear reaction-diffusion system** with three coupled species on a unit cube domain $Ω = [0,1]^3$. The system demonstrates various linear solver techniques for handling the resulting large, coupled nonlinear equations.

## Mathematical Formulation

The system solves the following coupled partial differential equations for three species $u_1$, $u_2$, and $u_3$:

```math
\frac{\partial u_i}{\partial t} + \nabla \cdot \mathbf{J}_i + R_i = S_i \quad \text{in } Ω, \quad i = 1, 2, 3
```

### Reaction Terms
The nonlinear reaction terms create a cyclic coupling between the three species:

```math
\begin{aligned}
R_1 &= 100 u_1 - u_2^2\\
R_2 &= 100 u_2 - u_3^2\\  
R_3 &= 100 u_3 - u_1^2
\end{aligned}
```

### Flux Terms
The flux $J_i$ for each species involves nonlinear, cross-diffusion effects. For each edge connecting two nodes, the flux is computed as:

```math
\begin{aligned}
J_1 &= -\varepsilon (u_1+u_2)\nabla u_1^2 \\
J_2 &= -\varepsilon (u_2+u_3)\nabla u_2^2 \\
J_3 &= -\varepsilon (u_3+u_1)\nabla u_3^2
\end{aligned}
```

where $\varepsilon = 1.0$ is the diffusion parameter, and $u_i^{(1)}$, $u_i^{(2)}$ denote the values of species $i$ at the two nodes of an edge. The flux combines cross-diffusion (species coupling through averaged concentrations) with nonlinear concentration-dependent transport.

### Source Terms
All three species have identical source terms with a Gaussian profile centered at $(0.5, 0.5)$:

```math
S_i = \exp\left(-20\left((x_1 - 0.5)^2 + (x_2 - 0.5)^2\right)\right), \quad i = 1, 2, 3
```

### Boundary Conditions
Dirichlet boundary conditions are applied:
- On regions 2 and 4: $u_i = 1$ for all species $i = 1, 2, 3$
- Homogeneous Neumann BC otherwise

This example serves as benchmark code for comparing the efficiency of different linear solver approaches for large-scale, coupled nonlinear PDE systems.

=#

module DevEx004_EquationBlock3D

## under development

using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize
using LinearSolve
using ExtendableSparse
using ExtendableSparse: ILUZeroPreconBuilder, JacobiPreconBuilder
using SparseArrays
import AMGCLWrap, Metis
using AlgebraicMultigrid
using LinearAlgebra
import Pardiso
using Test

function main(; nref = 0, Plotter = nothing, npart = 20, assembly = :edgewise, tol = 1.0e-7, kwargs...)
    # Create a uniform 3D grid on the unit cube [0,1]³
    # Grid resolution increases exponentially with nref parameter
    X = range(0, 1.0, length = 5 * 2^nref + 1)
    grid = simplexgrid(X, X, X)

    # Enable parallel processing if multiple threads are available
    # Partition the grid using METIS for load balancing across threads
    if Threads.nthreads() > 1
        grid = partition(grid, PlainMetisPartitioning(; npart); nodes = true, edges = true)
    end
    nn = num_nodes(grid)

    # Diffusion parameter ε appearing in the flux terms
    eps = 1.0

    # Define the nonlinear reaction terms creating cyclic coupling between species
    # R₁ = 100u₁ - u₂², R₂ = 100u₂ - u₃², R₃ = 100u₃ - u₁²
    function reaction(f, u, node, data)
        f[1] = 100 * u[1] - u[2]^2
        f[2] = 100 * u[2] - u[3]^2
        f[3] = 100 * u[3] - u[1]^2
        return nothing
    end

    # Define nonlinear cross-diffusion flux terms
    # J₁ = -ε(u₁+u₂)∇u₁², J₂ = -ε(u₂+u₃)∇u₂², J₃ = -ε(u₃+u₁)∇u₃²
    # Cross-diffusion coefficients are averaged between edge nodes
    function flux(f, u, edge, data)
        d1 = (u[2, 1] + u[2, 2]) / 2  # Average u₂ concentration for u₁ flux
        d2 = (u[3, 1] + u[3, 2]) / 2  # Average u₃ concentration for u₂ flux
        d3 = (u[1, 1] + u[1, 2]) / 2  # Average u₁ concentration for u₃ flux
        f[1] = eps * (d1 + d2) * (u[1, 1]^2 - u[1, 2]^2)
        f[2] = eps * (d2 + d3) * (u[2, 1]^2 - u[2, 2]^2)
        f[3] = eps * (d3 + d1) * (u[3, 1]^2 - u[3, 2]^2)
        return nothing
    end

    # Define Gaussian source terms centered at (0.5, 0.5) for all species
    # S = exp(-20((x₁-0.5)² + (x₂-0.5)²))
    function source(f, node, data)
        x1 = node[1] - 0.5
        x2 = node[2] - 0.5
        f[1] = exp(-20.0 * (x1^2 + x2^2))
        f[2] = f[1]
        f[3] = f[1]
        return nothing
    end

    # Define storage term for time-dependent problems (identity operator)
    function storage(f, u, node, data)
        f .= u
        return nothing
    end

    # Set Dirichlet boundary conditions: u = 1 on regions 2 and 4 for all species
    # Homogeneous Neumann conditions elsewhere
    function bcondition(f, u, node, data)
        for species in (1, 2, 3)
            boundary_dirichlet!(
                f,
                u,
                node;
                species,
                region = 2,
                value = 1,
            )
            boundary_dirichlet!(
                f,
                u,
                node;
                species,
                region = 4,
                value = 1,
            )
        end
        return nothing
    end

    # Assemble the VoronoiFVM system with all physics components
    # Three coupled species with nonlinear cross-diffusion and reaction terms
    sys = VoronoiFVM.System(
        grid; reaction, flux, source, storage, bcondition, assembly,
        species = [1, 2, 3],
    )

    # Choose platform-appropriate direct solver for reference solution
    # UMFPACK on macOS, MKL Pardiso on other platforms
    if Sys.isapple()
        direct_solver = LinearSolve.UMFPACKFactorization()
        println("Direct solver: UMFPACK")
    else
        direct_solver = LinearSolve.MKLPardisoFactorize()
        println("Direct solver: MKL Pardiso")
    end

    # Solve with direct method to obtain reference solution
    @time "Direct solution" direct_sol = solve(sys; inival = 0.5, method_linear = direct_solver, kwargs...)
    println()

    ok = true

    # Benchmark 1: AMGCL algebraic multigrid solver with block structure
    # Exploits the 3-species block structure for efficient multigrid
    @time "AMGCL" amgcl_sol = solve(sys; inival = 0.5, method_linear = AMGCLWrap.AMGSolverAlgorithm(blocksize = 3), kwargs...)
    @show n = norm(amgcl_sol - direct_sol, Inf)
    ok = ok && n < tol
    println()

    # Benchmark 2: ILU(0) preconditioned GMRES
    # Block ILU(0) factorization respects 3-species coupling structure
    @time "ilu0-gmres" ilu0_sol = solve(
        sys; inival = 0.5,
        method_linear = KrylovJL_GMRES(precs = ILUZeroPreconBuilder(blocksize = 3)),
        keepcurrent_linear = false,
        kwargs...
    )
    @show n = norm(ilu0_sol - direct_sol, Inf)
    ok = ok && n < tol
    println()

    # Benchmark 3: Direct solver as preconditioner for GMRES
    # Uses exact factorization from first Newton step as preconditioner in the subsequent steps
    @time "Direct-gmres" pgmres_sol = solve(
        sys; inival = 0.5,
        method_linear = KrylovJL_GMRES(precs = LinearSolvePreconBuilder(direct_solver)),
        keepcurrent_linear = false,
        kwargs...
    )
    @show n = norm(pgmres_sol - direct_sol, Inf)
    ok = ok && n < tol
    println()

    # Benchmark 4: Block GMRES with direct solver preconditioner
    # Applies block structure to GMRES iteration with direct preconditioning
    # Partitioning separates the three species into distinct blocks
    @time "Direct-blockgmres" pbgmres_sol = solve(
        sys; inival = 0.5,
        method_linear = KrylovJL_GMRES(
            precs = BlockPreconBuilder(
                precs = LinearSolvePreconBuilder(direct_solver),
                partitioning = A -> [1:3:size(A, 1), 2:3:size(A, 1), 3:3:size(A, 1)]
            )
        ),
        keepcurrent_linear = false,
        kwargs...
    )
    @show n = norm(pbgmres_sol - direct_sol, Inf)
    ok = ok && n < tol
    println()

    # Benchmark 5: Block GMRES with AMGCL preconditioner
    # Combines block structure with algebraic multigrid preconditioning
    # More scalable than direct methods for large 3D systems
    if !Sys.isapple() # getting some unclear OMP problems here...
        @time "AMGCL-blockgmres" amgclbgmres_sol = solve(
            sys; inival = 0.5,
            method_linear = KrylovJL_GMRES(
                precs = BlockPreconBuilder(
                    precs = AMGCLWrap.AMGPreconBuilder(),
                    partitioning = A -> [1:3:size(A, 1), 2:3:size(A, 1), 3:3:size(A, 1)]
                )
            ),
            keepcurrent_linear = false,
            kwargs...
        )
        @show n = norm(amgclbgmres_sol - direct_sol, Inf)
        ok = ok && n < tol
        println()
    end
    # Benchmark 6: Block GMRES with native Julia AMG preconditioner
    # Uses AlgebraicMultigrid.jl instead of AMGCL for comparison
    # Smoothed aggregation AMG with block structure
    @time "AMG-blockgmres" amgbgmres_sol = solve(
        sys; inival = 0.5,
        method_linear = KrylovJL_GMRES(
            precs = BlockPreconBuilder(
                precs = AlgebraicMultigrid.SmoothedAggregationPreconBuilder(),
                partitioning = A -> [1:3:size(A, 1), 2:3:size(A, 1), 3:3:size(A, 1)]
            )
        ),
        keepcurrent_linear = false,
        kwargs...
    )
    @show n = norm(amgbgmres_sol - direct_sol, Inf)
    ok = ok && n < tol
    println()

    return ok
end


using Test
function runtests()
    @test main() == true
    return nothing
end

end
