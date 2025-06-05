#=
# 108: 1D Nonlinear Diffusion equation with ODE
([source code](@__SOURCE_URL__))

Solve the nonlinear diffusion equation

```math
\partial_t u -\Delta u^m = 0
```
in $\Omega=(-1,1)$ with homogeneous Neumann boundary conditions using the implicit Euler method.

This equation is also called  "porous medium equation". 
The Barenblatt solution is an exact solution of this problem which for m>1 has a finite support.
We initialize this problem with the exact solution for $t=t_0=0.001$.

(see Barenblatt, G. I. "On nonsteady motions of gas and fluid in porous medium." Appl. Math. and Mech.(PMM) 16.1 (1952): 67-78.)
=#

module Example108_OrdinaryDiffEq1D
using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize
using OrdinaryDiffEqRosenbrock

function barenblatt(x, t, m)
    tx = t^(-1.0 / (m + 1.0))
    xx = x * tx
    xx = xx * xx
    xx = 1 - xx * (m - 1) / (2.0 * m * (m + 1))
    if xx < 0.0
        xx = 0.0
    end
    return tx * xx^(1.0 / (m - 1.0))
end

function main(;
        n = 20, m = 2, Plotter = nothing, verbose = false,
        unknown_storage = :sparse, tend = 0.01, assembly = :edgewise, solver = Rosenbrock23()
    )

    ## Create a one-dimensional discretization
    h = 1.0 / convert(Float64, n / 2)
    X = collect(-1:h:1)
    grid = simplexgrid(X)

    ## Flux function which describes the flux
    ## between neighboring control volumes
    function flux!(f, u, edg, data)
        f[1] = u[1, 1]^m - u[1, 2]^m
        return nothing
    end

    ## Storage term
    function storage!(f, u, node, data)
        f[1] = u[1]
        return nothing
    end

    ## Create a physics structure
    physics = VoronoiFVM.Physics(;
        flux = flux!,
        storage = storage!
    )

    ## Create a finite volume system - either
    ## in the dense or  the sparse version.
    ## The difference is in the way the solution object
    ## is stored - as dense or as sparse matrix
    sys = VoronoiFVM.System(grid, physics; unknown_storage = unknown_storage, assembly = assembly)

    ## Add species 1 to region 1
    enable_species!(sys, 1, [1])

    ## Create a solution array
    inival = unknowns(sys)
    t0 = 0.001

    ## Broadcast the initial value
    inival[1, :] .= map(x -> barenblatt(x, t0, m), X)

    problem = ODEProblem(sys, inival, (t0, tend))
    odesol = solve(problem, solver)
    tsol = reshape(odesol, sys)


    p = GridVisualizer(; Plotter = Plotter, layout = (1, 1), fast = true)
    for i in 1:length(tsol)
        time = tsol.t[i]
        scalarplot!(
            p[1, 1], grid, tsol[1, :, i]; title = @sprintf("t=%.3g", time),
            color = :red, label = "numerical",
            markershape = :circle, markevery = 1
        )
        scalarplot!(
            p[1, 1], grid, map(x -> barenblatt(x, time, m), grid); clear = false,
            color = :green,
            label = "exact", markershape = :none
        )
        reveal(p)
        sleep(1.0e-2)
    end
    return sum(tsol.u[end])
end

using Test
function runtests()
    testval = 46.66666666671521
    @test main(; unknown_storage = :sparse, assembly = :edgewise) ≈ testval
    @test main(; unknown_storage = :dense, assembly = :edgewise) ≈ testval
    @test main(; unknown_storage = :sparse, assembly = :cellwise) ≈ testval
    @test main(; unknown_storage = :dense, assembly = :cellwise) ≈ testval
    return nothing
end

end
