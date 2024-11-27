#=

# 115: 1D heterogeneous catalysis
 ([source code](@__SOURCE_URL__))

Let $\Omega=(0,1)$, $\Gamma_1=\{0\}$, $\Gamma_2=\{1\}$
Regard a system of three species: $A,B,C$ and let 
$u_A=[A]$, $u_B=[B]$ and $u_C=[C]$ be their corresponding concentrations.

Species $A$ and $B$ exist in the interior of the domain, species $C$
lives a the boundary $\Gamma_1$.  We assume a heterogeneous reaction scheme
where $A$ reacts to $C$ and $C$ reacts to $B$:

```math
\begin{aligned}
      A &\leftrightarrow C\\
      C &\leftrightarrow B 
\end{aligned}
```
with reaction constants $k_{AC}^\pm$ and k_{BC}^\pm$.

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
R_{AC}(u_A, u_C)&=k_{AC}^+ u_A(1-u_C) - k_{AC}^-u_C\\
R_{BC}(u_C, u_B)&=k_{BC}^+ u_B(1-u_C) - k_{BC}^-u_C\\
- D_A \nabla u_A  + S R_{AC}(u_A, u_C)& =0 \\
- D_B \nabla u_B  + S R_{BC}(u_B, u_C)& =0 \\
\partial_t C  - R_{AC}(u_A, u_C) - R_{BC}(u_B, u_C) &=0
\end{aligned}
```

=#

module Example115_HeterogeneousCatalysis1D
using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize
using LinearAlgebra
using OrdinaryDiffEqRosenbrock
using SciMLBase: NoInit

function main(;
        n = 10, Plotter = nothing, verbose = false, tend = 1,
        unknown_storage = :sparse, assembly = :edgewise,
        diffeq = false
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

    ## Diffusion flux for species A and B
    D_A = 1.0
    D_B = 1.0e-2
    function flux!(f, u, edge, data)
        f[iA] = D_A * (u[iA, 1] - u[iA, 2])
        f[iB] = D_B * (u[iB, 1] - u[iB, 2])
        return nothing
    end

    ## Storage term of species A and B
    function storage!(f, u, node, data)
        f[iA] = u[iA]
        f[iB] = u[iB]
        return nothing
    end

    ## Source term for species a around 0.5
    function source!(f, node, data)
        x1 = node[1] - 0.5
        f[iA] = exp(-100 * x1^2)
        return nothing
    end

    ## Reaction constants (p = + , m = -)
    ## Chosen to prefer path A-> C -> B
    ## More over, A reacts faster than to C than C to B
    ## leading to "catalyst poisoning", i.e. C taking up most of the
    ## available catalyst sites
    kp_AC = 100.0
    km_AC = 1.0

    kp_BC = 0.1
    km_BC = 1.0

    S = 0.01

    R_AC(u_A, u_C) = kp_AC * u_A * (1 - u_C) - km_AC * u_C
    R_BC(u_B, u_C) = kp_BC * u_B * (1 - u_C) - km_BC * u_C

    function breaction!(f, u, node, data)
        if node.region == 1
            f[iA] = S * R_AC(u[iA], u[iC])
            f[iB] = S * R_BC(u[iB], u[iC])
            f[iC] = -R_BC(u[iB], u[iC]) - R_AC(u[iA], u[iC])
        end
        return nothing
    end

    ## This is for the term \partial_t u_C at the boundary
    function bstorage!(f, u, node, data)
        if node.region == 1
            f[iC] = u[iC]
        end
        return nothing
    end

    physics = VoronoiFVM.Physics(;
        breaction = breaction!,
        bstorage = bstorage!,
        flux = flux!,
        storage = storage!,
        source = source!
    )

    sys = VoronoiFVM.System(grid, physics; unknown_storage = unknown_storage)

    ## Enable species in bulk resp
    enable_species!(sys, iA, [1])
    enable_species!(sys, iB, [1])

    ## Enable surface species
    enable_boundary_species!(sys, iC, [1])

    ## Set Dirichlet bc for species B on \Gamma_2
    boundary_dirichlet!(sys, iB, 2, 0.0)

    ## Initial values
    inival = unknowns(sys)
    inival .= 0.0
    U = unknowns(sys)

    tstep = 0.01
    time = 0.0

    ## Data to store surface concentration vs time

    p = GridVisualizer(; Plotter = Plotter, layout = (3, 1))
    if diffeq
        inival = unknowns(sys, inival = 0)
        problem = ODEProblem(sys, inival, (0, tend))
        ## use fixed timesteps just for the purpose of CI
        odesol = solve(problem, Rosenbrock23(); initializealg = NoInit(), dt = tstep, adaptive = false)
        tsol = reshape(odesol, sys)
    else
        control = fixed_timesteps!(VoronoiFVM.NewtonControl(), tstep)
        tsol = solve(sys; inival, times = [0, tend], control, verbose = verbose)
    end

    p = GridVisualizer(; Plotter = Plotter, layout = (3, 1), fast = true)
    for it in 1:length(tsol)
        time = tsol.t[it]
        scalarplot!(
            p[1, 1], grid, tsol[iA, :, it]; clear = true,
            title = @sprintf("[A]: (%.3f,%.3f)", extrema(tsol[iA, :, it])...)
        )
        scalarplot!(
            p[2, 1], grid, tsol[iB, :, it]; clear = true,
            title = @sprintf("[B]: (%.3f,%.3f)", extrema(tsol[iB, :, it])...)
        )
        scalarplot!(
            p[3, 1], tsol.t[1:it], tsol[iC, 1, 1:it]; title = @sprintf("[C]"),
            clear = true, show = true
        )
    end

    return tsol[iC, 1, end]
end

using Test
function runtests()
    testval = 0.87544440641274
    testvaldiffeq = 0.8891082547874963
    @test isapprox(main(; unknown_storage = :sparse, assembly = :edgewise), testval; rtol = 1.0e-12)
    @test isapprox(main(; unknown_storage = :dense, assembly = :edgewise), testval; rtol = 1.0e-12)
    @test isapprox(main(; unknown_storage = :sparse, assembly = :cellwise), testval; rtol = 1.0e-12)
    @test isapprox(main(; unknown_storage = :dense, assembly = :cellwise), testval; rtol = 1.0e-12)

    @test isapprox(main(; diffeq = true, unknown_storage = :sparse, assembly = :edgewise), testvaldiffeq; rtol = 1.0e-12)
    @test isapprox(main(; diffeq = true, unknown_storage = :dense, assembly = :edgewise), testvaldiffeq; rtol = 1.0e-12)
    @test isapprox(main(; diffeq = true, unknown_storage = :sparse, assembly = :cellwise), testvaldiffeq; rtol = 1.0e-12)
    @test isapprox(main(; diffeq = true, unknown_storage = :dense, assembly = :cellwise), testvaldiffeq; rtol = 1.0e-12)
    return nothing
end
end
