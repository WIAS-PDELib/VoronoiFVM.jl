# # 206: 2D Joule heating
# ([source code](@__SOURCE_URL__))
#=

```math
\begin{aligned}
-\nabla \left\cdot (\kappa(T) \nabla \phi\right) &= 0\\
\partial_t (cT) - \nabla\cdot \left(\lambda \nabla T\right) &= \kappa(T) |\nabla \phi|^2\\
\kappa(T)&= \kappa_0 exp(\alpha(T-T0))
\end{aligned}
```
The discretization uses the approach developed in
A. Bradji, R. Herbin, [DOI 10.1093/imanum/drm030](https://doi.org/10.1093/imanum/drm030).
=#

module Example206_JouleHeat

using Printf
using VoronoiFVM
using ExtendableGrids
using ExtendableSparse
using GridVisualize
using LinearAlgebra
using SimplexGridFactory
using LinearSolve
import Triangulate
import Metis

function main(;
        nref = 0, Plotter = nothing, verbose = "and", unknown_storage = :sparse, assembly = :edgewise,
        ythin = 0.25
    )

    ## Create grid
    b = SimplexGridBuilder(; Generator = Triangulate)
    p00 = point!(b, 0, 0)
    p30 = point!(b, 3, 0)
    p32 = point!(b, 3, 1)
    p21 = point!(b, 2, ythin)
    p11 = point!(b, 1, ythin)
    p02 = point!(b, 0, 1)

    facetregion!(b, 4)
    facet!(b, p00, p30)
    facetregion!(b, 2)
    facet!(b, p30, p32)
    facetregion!(b, 3)
    facet!(b, p32, p21)
    facet!(b, p21, p11)
    facet!(b, p11, p02)
    facetregion!(b, 1)
    facet!(b, p02, p00)

    grid = simplexgrid(b; maxvolume = 0.01 * 4.0^(-nref))
    grid = partition(grid, PlainMetisPartitioning(npart = 20); nodes = true, edges = true)
    @show grid
    ## Describe problem
    iϕ::Int = 1
    iT::Int = 2
    κ0::Float64 = 1
    α::Float64 = 1
    T0::Float64 = 0.5
    λ::Float64 = 1
    c::Float64 = 1

    function storage!(y, u, node, data)
        y[iT] = c * u[iT]
        return nothing
    end

    κ(T) = κ0 * exp(α * (T - T0))

    function flux!(y, u, edge, data)
        y[iϕ] = κ(y[iT]) * (u[iϕ, 1] - u[iϕ, 2])
        y[iT] = λ * (u[iT, 1] - u[iT, 2])
        return nothing
    end

    ## The convention in VoronoiFVM.jl is to have all terms depending on the solution
    ## on the left hand side of the equation. That is why we have the minus sign here.
    function jouleheat!(y, u, edge, data)
        y[iT] = -κ(y[iT]) * (u[iϕ, 1] - u[iϕ, 2]) * (u[iϕ, 1] - u[iϕ, 2])
        return nothing
    end

    function bcondition!(y, u, node, data)
        boundary_dirichlet!(y, u, node; species = iϕ, region = 1, value = -10)
        boundary_dirichlet!(y, u, node; species = iϕ, region = 2, value = 10)

        boundary_robin!(y, u, node; species = iT, region = 1, value = T0, factor = 0.5)
        boundary_robin!(y, u, node; species = iT, region = 2, value = T0, factor = 0.5)
        boundary_robin!(y, u, node; species = iT, region = 3, value = T0, factor = 0.5)
        boundary_robin!(y, u, node; species = iT, region = 4, value = T0, factor = 0.5)
        return nothing
    end

    sys = VoronoiFVM.System(
        grid; bcondition = bcondition!, flux = flux!,
        edgereaction = jouleheat!, storage = storage!,
        species = [iϕ, iT], assembly = assembly
    )

    sol = solve(
        sys; verbose,
        method_linear = KrylovJL_BICGSTAB(precs = LinearSolvePreconBuilder(UMFPACKFactorization())),
        keepcurrent_linear = false
    )

    vis = GridVisualizer(; Plotter, layout = (2, 1))
    scalarplot!(vis[1, 1], grid, sol[iϕ, :]; title = "ϕ", colormap = :bwr)
    scalarplot!(vis[2, 1], grid, sol[iT, :]; title = "T", colormap = :hot)
    reveal(vis)
    return norm(sol, Inf)
end

using Test
function runtests()
    testval = 24.639120035942938
    @test main(; assembly = :edgewise) ≈ testval &&
        main(; assembly = :cellwise) ≈ testval
    return nothing
end
end
