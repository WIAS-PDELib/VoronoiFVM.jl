#=
# 105: 1D Nonlinear Poisson equation 
([source code](@__SOURCE_URL__))

Solve the nonlinear Poisson equation

```math
-\nabla \varepsilon \nabla u + e^{u}-e^{-u} = f
```
in $\Omega=(0,1)$ with boundary condition $u(0)=0$ and $u(1)=1$ with 
```math
f(x)=
    \begin{cases}
    1&,x>0.5\\
    -1&, x<0.5
    \end{cases}.
```

This stationary problem is an example of a nonlinear Poisson equation or Poisson-Boltzmann equation.
Such equation occur e.g. in simulations of electrochemical systems and semiconductor devices.

=#

module Example105_NonlinearPoisson1D
using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize

function main(; n = 10, Plotter = nothing, verbose = false, unknown_storage = :sparse, assembly = :edgewise)

    ## Create a one-dimensional discretization
    h = 1.0 / convert(Float64, n)
    grid = simplexgrid(collect(0:h:1))

    ## A parameter which is "passed" to the flux function via scope
    ϵ = 1.0e-3

    ## Flux function which describes the flux
    ## between neighboring control volumes
    function flux!(f, u, edge, data)
        f[1] = ϵ * (u[1, 1] - u[1, 2])
        return nothing
    end

    ## Source term
    function source!(f, node, data)
        if node[1] <= 0.5
            f[1] = 1
        else
            f[1] = -1
        end
        return nothing
    end

    ## Reaction term
    function reaction!(f, u, node, data)
        f[1] = exp(u[1]) - exp(-u[1])
        return nothing
    end

    ## Create a physics structure
    physics = VoronoiFVM.Physics(;
        flux = flux!,
        source = source!,
        reaction = reaction!
    )

    ## Create a finite volume system - either
    ## in the dense or  the sparse version.
    ## The difference is in the way the solution object
    ## is stored - as dense or as sparse matrix

    sys = VoronoiFVM.System(grid, physics; unknown_storage = unknown_storage, assembly = assembly)

    ## Add species 1 to region 1
    enable_species!(sys, 1, [1])

    ## Set boundary conditions
    boundary_dirichlet!(sys, 1, 1, 0.0)
    boundary_dirichlet!(sys, 1, 2, 1.0)

    ## Create a solution array
    inival = unknowns(sys; inival = 0.5)

    ## Stationary solution of the problem
    solution = solve(sys; inival, verbose)

    scalarplot(grid, solution[1, :]; title = "Nonlinear Poisson", Plotter = Plotter)

    return sum(solution)
end

using Test
function runtests()
    testval = 1.5247901344230088
    @test main(; unknown_storage = :sparse, assembly = :edgewise) ≈ testval &&
        main(; unknown_storage = :dense, assembly = :edgewise) ≈ testval &&
        main(; unknown_storage = :sparse, assembly = :cellwise) ≈ testval &&
        main(; unknown_storage = :dense, assembly = :cellwise) ≈ testval
    return nothing
end

end
