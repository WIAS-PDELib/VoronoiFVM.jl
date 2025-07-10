#=

# 225: Terminal flux calculation via test functions, nD
([source code](@__SOURCE_URL__))

After calculating solutions based on the finite volume method, it
may be interesting to obtain information about the solution besides of the graphical representation.

Here, we focus on the following data:
- integrals of the solution
- flux through parts of the boundary

Let us define the following reaction - diffusion system in a domain $\Omega$:

```math
\begin{aligned}
\partial_t u_1 - \nabla \cdot \nabla u_1 + r(u_1, u_2) &= f=1.0\\
\partial_t u_2 - \nabla \cdot \nabla u_1 - r(u_1, u_2) &= 0
\end{aligned}
```

with boundary conditions $u_2=0$ on $\Gamma_2\subset\partial\Omega$ and
$r(u_1,u_2)=u_1 + 0.1 u_2$

The source $f$ creates species $u_1$ which reacts to $u_2$, $u_2$ then leaves
the domain at boundary $\Gamma_2$.

### Stationary problem

For the stationary problem, we have the following flux balances derived from the equations and from Gauss theorem:

```math
\begin{aligned}
\int_\Omega r(u_1,u_2) d\omega &= \int_\Omega f d\omega \\
\int_\Omega -r(u_1,u_2) d\omega &= \int_{\Gamma_2} \nabla u \cdot \vec n ds \\
\end{aligned}
```

The volume integrals can be approximated based on the finite volume subdivision
$\Omega=\cup_{i\in \mathcal N} \omega_i$:

```math
\begin{aligned}
\int_\Omega r(u_1,u_2) d\omega &\approx \sum_{i\in \mathcal N} |\omega_i| r(u_{1,i},u_{2,i})\\
\int_\Omega f d\omega &\approx \sum_{i\in \mathcal N} |\omega_i| f_i
\end{aligned}
```

But what about  the boundary integral ?  Here, we use a  trick to cast
the surface  integral to the  integral to  a volume integral  with the
help of a test function:

Let $T(x)$ be the solution of the Laplace problem $-\nabla^2 T =0$ in $\Omega$ and the boundary conditions

```math
\begin{aligned}
T &=0 \quad \text{at}\; \Gamma_4\\
T &=1 \quad \text{at}\; \Gamma_2\\
\partial_n T &=0\quad \text{at}\;  \Gamma_1,\Gamma_3
\end{aligned}
```

Write ``\vec j=-\nabla u``. and assume $\nabla\cdot \vec j + r =f$.

```math
\begin{aligned}
\int_{\Gamma_2} \vec j \cdot \vec n ds&=\int_{\Gamma_2} T\vec j \cdot \vec n ds \quad \text{due to}\; T=1\; \text{on}\; \Gamma_2\\
 &=\int_{\partial\Omega}  T\vec j \cdot \vec n ds\quad \text{due to}\; T=0\; \text{on}\; \Gamma_4, \quad\vec j\cdot \vec n=0\; \text{on}\; \Gamma_1, \Gamma_3\\
&= \int_\Omega \nabla \cdot (T \vec j) d\omega \quad \text{(Gauss)}\\
&= \int_\Omega \nabla T \cdot \vec j d\omega + \int_\Omega T \nabla\cdot j d\omega\\
&=  \int_\Omega \nabla T \cdot \vec j d\omega + \int_\Omega T(f-r)dω\\
\end{aligned}
```

and we approximate

```math
\begin{aligned}
\int_\Omega \nabla T \cdot \vec j d\omega \approx \sum_{k,l}
\frac{|\omega_k\cap\omega_l|}{h_{k,l}}g(u_k, u_l) (T_k-T_l)
\end{aligned}
```

where the sum runs over pairs of neighboring control volumes.

The `integrate` method with a  test function parameter returns a value
for each species, the sign convention assumes that species leaving the
domain lead to negative values.

### Transient problem

The amount  of species created via  the source term (measured  in `F`)
integrated  over time  should be  equal to  the sum  of the  amount of
species left in  the domain at the  very end of the  evolution and the
amount of species which left the domain:

$\int_{t_0}^{t_{end}} \int_\Omega f d\omega dt= \int_\Omega (u_1+u_2)dω + \int_{t_0}^{t_{end}} \int_{\Gamma_2} \nabla u_2 \cdot \vec n ds$

Literature references:

- H. Gajewski "Analysis und Numerik von Ladungstransport in Halbleitern", WIAS Berlin, Report No.6
- Yoder, P. D., K. Gärtner, and W. Fichtner. "A generalized Ramo–Shockley theorem for classical to quantum transport at arbitrary frequencies." Journal of Applied Physics 79.4 (1996): 1951-1954.
- P. Farrell, N. Rotundo, D. H. Doan, M. Kantner, J. Fuhrmann, and T. Koprucki, "Numerical methods for drift-diffusion models", in Handbook of optoelectronic device modeling and simulation: Lasers, modulators, photodetectors, solar cells, and numerical methods, vol. 2, J. Piprek, Ed. Boca Raton: CRC Press, 2017, pp. 733–771.

=#

module Example225_TestFunctions2D

using VoronoiFVM, GridVisualize, ExtendableGrids

function main(;
        n = 10, Plotter = nothing, verbose = false, unknown_storage = :sparse, assembly = :edgewise,
        dim = 2, tend = 5, dt = 0.2
    )
    n = [101, 21, 5]
    X = collect(range(0.0, 1; length = n[dim]))
    if dim == 1
        grid = simplexgrid(X)
        Γ_where_T_equal_1 = [2]
        Γ_where_T_equal_0 = [1]
    elseif dim == 2
        grid = simplexgrid(X, X)
        Γ_where_T_equal_1 = [2]
        Γ_where_T_equal_0 = [4]
    elseif dim == 3
        grid = simplexgrid(X, X, X)
        Γ_where_T_equal_1 = [2]
        Γ_where_T_equal_0 = [4]
    end

    function storage(f, u, node, data)
        f .= u
        return nothing
    end

    function flux(f, u, edge, data)
        f[1] = u[1, 1] - u[1, 2]
        f[2] = u[2, 1] - u[2, 2]
        return nothing
    end

    r(u1, u2) = u1 - 0.1 * u2

    function reaction(f, u, node, data)
        f[1] = r(u[1], u[2])
        f[2] = -r(u[1], u[2])
        return nothing
    end

    function source(f, node, data)
        f[1] = 1.0
        return nothing
    end

    physics = VoronoiFVM.Physics(;
        flux = flux,
        storage = storage,
        reaction = reaction,
        source = source
    )

    system = VoronoiFVM.System(grid, physics; assembly = assembly)

    enable_species!(system, 1, [1])
    enable_species!(system, 2, [1])
    boundary_dirichlet!(system, 2, 2, 0.0)

    sol = solve(system; inival = 0.0)

    vis = GridVisualizer(;
        Plotter = Plotter, layout = (1, 2), resolution = (600, 300),
        fignumber = 1
    )
    scalarplot!(vis[1, 1], grid, sol[1, :]; flimits = (0, 1.5), title = "u_1")
    scalarplot!(vis[1, 2], grid, sol[2, :]; flimits = (0, 1.5), title = "u_2", show = true)

    """
        The `integrate` method of `VoronoiFVM`  provides a possibility to calculate
        the volume integral of a function of a solution as described above.
        It returns a `num_specie` x `num_regions` matrix of the integrals
        of the function of the unknowns over the different subdomains (here, we have only one):
    """

    """
        Amount of u_1 and u_2 in the domain aka integral over identity storage function:
    """
    U = integrate(system, storage, sol)

    """
    Amount of species created by source term per unit time:
    """
    F = integrate(system, (f, u, node, data) -> source(f, node, data), sol)

    """
    Amount of  reaction per unit time:
    """
    R = integrate(system, reaction, sol)

    tf = TestFunctionFactory(system)
    T = testfunction(tf, Γ_where_T_equal_0, Γ_where_T_equal_1)

    I = integrate(system, T, sol)

    t0 = 0.0

    control = fixed_timesteps!(VoronoiFVM.NewtonControl(), dt)

    tsol = solve(system; inival = 0.0, times = [t0, tend], control)

    vis1 = GridVisualizer(;
        Plotter = Plotter, layout = (1, 2), resolution = (600, 300),
        fignumber = 4
    )

    for i in 1:length(tsol)
        sol = tsol.u[i]
        scalarplot!(vis1[1, 1], grid, sol[1, :]; flimits = (0, 1.5), clear = true)
        scalarplot!(vis1[1, 2], grid, sol[2, :]; flimits = (0, 1.5), show = true)
    end

    outflow_rate = integrate(system, T, tsol, rate = false)

    vis2 = GridVisualizer(;
        Plotter = Plotter, layout = (1, 1), resolution = (600, 300), legend = :lc,
        fignumber = 2
    )
    scalarplot!(vis2[1, 1], [0, tend], -[I[2], I[2]]; label = "stationary", clear = true, color = :red)
    scalarplot!(vis2[1, 1], tsol.t[2:end], -outflow_rate[2, :]; label = "transient", clear = false, show = true, color = :blue)

    all_outflow = 0.0
    for i in 1:(length(tsol) - 1)
        all_outflow -= outflow_rate[2, i] * (tsol.t[i + 1] - tsol.t[i])
    end

    Uend = integrate(system, storage, tsol.u[end])
    isapprox(F[1], R[1]; rtol = 1.0e-12) ? true : return false
    isapprox(I[1], 0.0; atol = 1.0e-12) ? true : return false
    isapprox(R[2], I[2]; rtol = 1.0e-12) ? true : return false
    return isapprox(F[1] * (tend - t0), (Uend[1] + Uend[2] + all_outflow); rtol = 1.0e-12) ? true :
        return false
end

using Test
function runtests()
    @test main(; dim = 1, unknown_storage = :sparse, assembly = :edgewise)
    @test main(; dim = 1, unknown_storage = :dense, assembly = :edgewise)
    @test main(; dim = 2, unknown_storage = :sparse, assembly = :edgewise)
    @test main(; dim = 2, unknown_storage = :dense, assembly = :edgewise)
    @test main(; dim = 3, unknown_storage = :sparse, assembly = :edgewise)
    @test main(; dim = 3, unknown_storage = :dense, assembly = :edgewise)

    @test main(; dim = 1, unknown_storage = :sparse, assembly = :cellwise)
    @test main(; dim = 1, unknown_storage = :dense, assembly = :cellwise)
    @test main(; dim = 2, unknown_storage = :sparse, assembly = :cellwise)
    @test main(; dim = 2, unknown_storage = :dense, assembly = :cellwise)
    @test main(; dim = 3, unknown_storage = :sparse, assembly = :cellwise)
    @test main(; dim = 3, unknown_storage = :dense, assembly = :cellwise)
    return nothing
end

end
