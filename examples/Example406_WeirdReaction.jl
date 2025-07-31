#=

# 406: 1D Weird Surface Reaction
 ([source code](@__SOURCE_URL__))

Species $A$ and $B$ exist in the interior of the do    function generic_operator!(f, u, sys)
        f .= 0
        idx = unknown_indices(unknowns(sys))
        f[idx[problem_data.iC, 1]] = u[idx[problem_data.iC, 1]] +
            0.1 * (u[idx[problem_data.iA, 1]] - u[idx[problem_data.iA, 2]]) / (problem_data.X[2] - problem_data.X[1])
        return nothing
    end

    # If we know the sparsity pattern, we can here create a
    # sparse matrix with values set to 1 in the nonzero
    # slots. This allows to circumvent the
    # autodetection which may takes some time.
    function generic_operator_sparsity(sys)
        idx = unknown_indices(unknowns(sys))
        sparsity = spzeros(num_dof(sys), num_dof(sys))
        sparsity[idx[problem_data.iC, 1], idx[problem_data.iC, 1]] = 1
        sparsity[idx[problem_data.iA, 1], idx[problem_data.iA, 1]] = 1
        sparsity[idx[problem_data.iC, 1], idx[problem_data.iA, 2]] = 1
        return sparsity
    lives at the boundary $\Gamma_1$.  We assume a heterogeneous reaction scheme
where $A$ reacts to $B$ with a rate depending on $\nabla A$ near the surface

```math
\begin{aligned}
      A &\leftrightarrow B\\
\end{aligned}
```

In $\Omega$, both $A$ and $B$ are transported through diffusion:

```math
\begin{aligned}
\partial_t u_B - \nabla\cdot D_A \nabla u_A & = f_A\\
\partial_t u_B - \nabla\cdot D_B \nabla u_B & = 0\\
\end{aligned}
```
Here, $f(x)$ is a source term creating $A$.
On $\Gamma_2$, we set boundary conditions
```math
\begin{aligned}
D_A \nabla u_A & = 0\\
u_B&=0
\end{aligned}
```
describing no normal flux for $A$ and zero concentration of $B$.
On $\Gamma_1$, we use the mass action law to describe the boundary reaction and
the evolution of the boundary concentration $C$. We assume that there is a limited
amount of surface sites $S$ for species C, so in fact A has to react with a free
surface site in order to become $C$ which reflected by the factor $1-u_C$. The same
is true for $B$.
```math
\begin{aligned}
R_{AB}(u_A, u_B)&=k_{AB}^+exp(u_A'(0))u_A - k_{AB}^-exp(-u_A'(0))u_B\\
- D_A \nabla u_A  +  R_{AB}(u_A, u_B)& =0 \\
- D_B \nabla u_B  -  R_{AB}(u_A, u_B)& =0 \\
\end{aligned}
```

=#
module Example406_WeirdReaction
using Printf
using VoronoiFVM
using SparseArrays
using ExtendableGrids
using GridVisualize

## Problem data structure to avoid global variables
mutable struct ProblemData
    D_A::Float64      # Diffusion coefficient for species A
    D_B::Float64      # Diffusion coefficient for species B
    kp_AB::Float64    # Forward reaction constant A->B
    km_AB::Float64    # Backward reaction constant B->A
    iA::Int           # Species index for A
    iB::Int           # Species index for B
    iC::Int           # Species index for C
    X::Vector{Float64} # Grid coordinates (needed for gradient calculation)
end

function main(;
        n = 10,
        Plotter = nothing,
        verbose = false,
        tend = 1,
        unknown_storage = :sparse,
        autodetect_sparsity = true
    )
    h = 1.0 / convert(Float64, n)
    X = collect(0.0:h:1.0)
    N = length(X)

    grid = simplexgrid(X)
    ## By default, \Gamma_1 at X[1] and \Gamma_2 is at X[end]

    ## Species numbers
    iA = 1
    iB = 2
    iC = 3

    ## Create problem data structure
    problem_data = ProblemData(1.0, 1.0e-2, 1.0, 0.1, iA, iB, iC, X)

    ## Diffusion flux for species A and B
    function flux!(f, u, edge, data)
        f[data.iA] = data.D_A * (u[data.iA, 1] - u[data.iA, 2])
        f[data.iB] = data.D_B * (u[data.iB, 1] - u[data.iB, 2])
        return nothing
    end

    ## Storage term of species A and B
    function storage!(f, u, node, data)
        f[data.iA] = u[data.iA]
        f[data.iB] = u[data.iB]
        return nothing
    end

    ## Source term for species a around 0.5
    function source!(f, node, data)
        x1 = node[1] - 0.5
        f[data.iA] = exp(-100 * x1^2)
        return nothing
    end

    function breaction!(f, u, node, data)
        if node.region == 1
            R = data.kp_AB * exp(u[data.iC]) * u[data.iA] - exp(-u[data.iC]) * data.km_AB * u[data.iB]
            f[data.iA] += R
            f[data.iB] -= R
        end
        return nothing
    end

    ## This generic operator works on the full solution seen as linear vector, and indexing
    ## into it needs to be performed with the help of idx (defined below for a solution vector)
    ## Its sparsity is detected automatically using SparsityDetection.jl
    ## Here, we calculate the gradient of u_A at the boundary and store the value in u_C which
    ## is then used as a parameter in the boundary reaction
    function generic_operator!(f, u, sys)
        f .= 0
        idx = unknown_indices(unknowns(sys))
        f[idx[problem_data.iC, 1]] = u[idx[problem_data.iC, 1]] +
            0.1 * (u[idx[problem_data.iA, 1]] - u[idx[problem_data.iA, 2]]) / (problem_data.X[2] - problem_data.X[1])
        return nothing
    end

    # If we know the sparsity pattern, we can here create a
    # sparse matrix with values set to 1 in the nonzero
    # slots. This allows to circumvent the
    # autodetection which may takes some time.
    function generic_operator_sparsity(sys)
        idx = unknown_indices(unknowns(sys))
        sparsity = spzeros(num_dof(sys), num_dof(sys))
        sparsity[idx[problem_data.iC, 1], idx[problem_data.iC, 1]] = 1
        sparsity[idx[problem_data.iC, 1], idx[problem_data.iA, 1]] = 1
        sparsity[idx[problem_data.iC, 1], idx[problem_data.iA, 2]] = 1
        return sparsity
    end

    if autodetect_sparsity
        physics = VoronoiFVM.Physics(;
            breaction = breaction!,
            generic = generic_operator!,
            flux = flux!,
            storage = storage!,
            source = source!,
            data = problem_data
        )
    else
        physics = VoronoiFVM.Physics(;
            breaction = breaction!,
            generic = generic_operator!,
            generic_sparsity = generic_operator_sparsity,
            flux = flux!,
            storage = storage!,
            source = source!,
            data = problem_data
        )
    end
    sys = VoronoiFVM.System(grid, physics; unknown_storage = unknown_storage)

    ## Enable species in bulk resp
    enable_species!(sys, problem_data.iA, [1])
    enable_species!(sys, problem_data.iB, [1])

    ## Enable surface species
    enable_boundary_species!(sys, problem_data.iC, [1])

    ## Set Dirichlet bc for species B on \Gamma_2
    boundary_dirichlet!(sys, problem_data.iB, 2, 0.0)

    ## Initial values
    U = unknowns(sys)
    U .= 0.0
    idx = unknown_indices(U)

    tstep = 0.01
    time = 0.0
    T = Float64[]
    u_C = Float64[]

    control = VoronoiFVM.SolverControl()
    control.verbose = verbose
    p = GridVisualizer(; Plotter = Plotter, layout = (2, 1))
    while time < tend
        time = time + tstep
        U = solve(sys; inival = U, time, tstep, control)
        if verbose
            @printf("time=%g\n", time)
        end
        ## Record  boundary pecies
        push!(T, time)
        push!(u_C, U[problem_data.iC, 1])

        scalarplot!(
            p[1, 1], grid, U[problem_data.iA, :]; label = "[A]",
            title = @sprintf(
                "max_A=%.5f max_B=%.5f u_C=%.5f", maximum(U[problem_data.iA, :]),
                maximum(U[problem_data.iB, :]), u_C[end]
            ), color = :red
        )
        scalarplot!(p[1, 1], grid, U[problem_data.iB, :]; label = "[B]", clear = false, color = :blue)
        scalarplot!(p[2, 1], copy(T), copy(u_C); label = "[C]", clear = true, show = true)
    end
    return U[problem_data.iC, 1]
end

using Test
function runtests()
    testval = 0.007027597470502758
    @test main(; unknown_storage = :sparse) ≈ testval &&
        main(; unknown_storage = :dense) ≈ testval &&
        main(; unknown_storage = :sparse, autodetect_sparsity = false) ≈ testval &&
        main(; unknown_storage = :dense, autodetect_sparsity = false) ≈ testval
    return nothing
end

end
