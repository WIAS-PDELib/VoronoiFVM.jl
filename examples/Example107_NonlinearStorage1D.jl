#=

# 107: 1D Nonlinear Storage 
([source code](@__SOURCE_URL__))

This equation comes from the transformation of the nonlinear diffuision equation.
```math
\partial_t u^\frac{1}{m} -\Delta u = 0
```
in $\Omega=(-1,1)$ with homogeneous Neumann boundary conditions.
We can derive an exact solution from the Barenblatt solution of the previous
example.

=#

module Example107_NonlinearStorage1D
using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize

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
        n = 20, m = 2.0, Plotter = nothing, verbose = false,
        unknown_storage = :sparse, tend = 0.01, tstep = 0.0001, assembly = :edgewise
    )

    ## Create a one-dimensional discretization
    h = 1.0 / convert(Float64, n / 2)
    X = collect(-1:h:1)
    grid = simplexgrid(X)

    ## Flux function which describes the flux
    ## between neighboring control volumes
    function flux!(f, u, edge, data)
        f[1] = u[1, 1] - u[1, 2]
        return nothing
    end

    ϵ = 1.0e-10

    ## Storage term
    ## This needs to be regularized as its derivative
    ## at 0 is infinity
    function storage!(f, u, node, data)
        f[1] = (ϵ + u[1])^(1.0 / m)
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
    solution = unknowns(sys)
    t0 = 0.001

    ## Broadcast the initial value
    inival[1, :] .= map(x -> barenblatt(x, t0, m)^m, X)

    ## Create solver control info
    control = VoronoiFVM.NewtonControl()
    control.verbose = verbose
    control.Δu_opt = 0.1
    control.force_first_step = true
    tsol = solve(sys; inival, times = [t0, tend], control)

    if Plotter != nothing
        p = GridVisualizer(; Plotter = Plotter, layout = (1, 1), fast = true)
        for i in 1:length(tsol)
            time = tsol.t[i]
            scalarplot!(
                p[1, 1], grid, tsol[1, :, i]; title = @sprintf("t=%.3g", time),
                color = :red, label = "numerical"
            )
            scalarplot!(
                p[1, 1], grid, map(x -> barenblatt(x, time, m)^m, grid); clear = false,
                color = :green, label = "exact"
            )
            reveal(p)
            sleep(1.0e-2)
        end
    end
    return sum(tsol.u[end])
end

using Test
function runtests()
    testval = 174.72418935404414
    @test main(; unknown_storage = :sparse, assembly = :edgewise) ≈ testval rtol = 1.0e-5
    @test main(; unknown_storage = :dense, assembly = :edgewise) ≈ testval rtol = 1.0e-5
    @test main(; unknown_storage = :sparse, assembly = :cellwise) ≈ testval rtol = 1.0e-5
    @test main(; unknown_storage = :dense, assembly = :cellwise) ≈ testval rtol = 1.0e-5

    return nothing
end

end
