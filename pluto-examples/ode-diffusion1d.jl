### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 6d29fa9a-54b9-490f-bb96-e44244522d27
# Convert this cell to markdown in order to enable Pluto's inbuilt package manager
if isdefined(Main, :PlutoRunner)
    using Pkg
    docsdir = joinpath(@__DIR__, "..", "docs")
    if isdir(docsdir)
        Pkg.activate(docsdir)
    end
end


# ╔═╡ 6d467640-b19c-4f77-845d-f9b4aca62104
begin
    using Test
    using Revise
    using Printf
    using ExtendableGrids
    using VoronoiFVM
    using OrdinaryDiffEqBDF
    using OrdinaryDiffEqRosenbrock
    using OrdinaryDiffEqSDIRK
    using LinearAlgebra
    using PlutoUI
    using DataStructures
    using GridVisualize
    if isdefined(Main, :PlutoRunner)
        using CairoMakie
        default_plotter!(CairoMakie)
    end
end

# ╔═╡ 02424193-41e8-4cec-8f52-6fd66173ace8
md"""
Solve the nonlinear diffusion equation

```math
\partial_t u -\Delta u^m = 0
```
in $\Omega=(-1,1)$ with homogeneous Neumann boundary conditions using the implicit Euler method.

This equation is also called  "porous medium equation". 
The Barenblatt solution 
```math
b(x,t)=\max\left(0,t^{-\alpha}\left(1-\frac{\alpha(m-1)r^2}{2dmt^{\frac{2\alpha}{d}}}\right)^{\frac{1}{m-1}}\right)
```

is an exact solution of this problem which for m>1 has a finite support.
We initialize this problem with the exact solution for $t=t_0=0.001$.

(see Barenblatt, G. I. "On nonsteady motions of gas and fluid in porous medium." Appl. Math. and Mech.(PMM) 16.1 (1952): 67-78.)

Here, we compare the implicit Euler approach in VoronoiFVM with the ODE solvers in [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/) and demonstrate the possibility to use VoronoiFVM to define differential
operators compatible with its [ODEFunction](https://diffeq.sciml.ai/stable/features/performance_overloads/#ODEFunction)
interface.

"""

# ╔═╡ 870b8b91-cd74-463e-b258-c092cd0af200
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

# ╔═╡ 78208d7b-9d79-415c-ab3a-85948251e635
function create_porous_medium_problem(n, m)
    h = 1.0 / convert(Float64, n / 2)
    X = collect(-1:h:1)
    grid = simplexgrid(X)

    function flux!(f, u, edge, data)
        f[1] = u[1, 1]^m - u[1, 2]^m
        return nothing
    end

    storage!(f, u, node, data) = f[1] = u[1]

    sys = VoronoiFVM.System(grid, flux = flux!, storage = storage!, species = 1)
    return sys, X
end

# ╔═╡ 4ef024a4-cb1d-443d-97fb-ab3a32a78ffd
begin
    function run_vfvm(; n = 20, m = 2, t0 = 0.001, tend = 0.01, tstep = 1.0e-6)
        sys, X = create_porous_medium_problem(n, m)
        inival = unknowns(sys)
        inival[1, :] .= map(x -> barenblatt(x, t0, m), X)
        sol = VoronoiFVM.solve(sys; inival, times = (t0, tend), Δt = tstep, Δu_opt = 0.01, Δt_min = tstep, store_all = true, log = true, reltol = 1.0e-3)
        err = norm(sol[1, :, end] - map(x -> barenblatt(x, tend, m), X))
        return sol, sys, err
    end
    run_vfvm(m = 2, n = 10) # "Precompile"
end;

# ╔═╡ 2a8ac57e-486d-4825-95ab-f0402b910dbd
diffeqmethods = OrderedDict(
    "Rosenbrock23 (Rosenbrock)" => Rosenbrock23,
    "QNDF2 (Like matlab's ode15s)" => QNDF2,
    "FBDF" => FBDF,
    "Implicit Euler" => ImplicitEuler
)

# ╔═╡ 9239409b-6de0-4157-8a35-412c909efa96
begin
    function run_diffeq(; n = 20, m = 2, t0 = 0.001, tend = 0.01, solver = nothing)
        sys, X = create_porous_medium_problem(n, m)
        inival = unknowns(sys)
        inival[1, :] .= map(x -> barenblatt(x, t0, m), X)
        state = VoronoiFVM.SystemState(sys)
        problem = ODEProblem(state, inival, (t0, tend))
        odesol = solve(problem, solver)
        sol = reshape(odesol, sys; state)
        err = norm(sol[1, :, end] - map(x -> barenblatt(x, tend, m), X))
        return sol, sys, err
    end
    for method in diffeqmethods
        run_diffeq(m = 2, n = 10, solver = method.second()) # "Precompile"
    end
end;

# ╔═╡ 3a004ab9-2705-4f5c-8e6e-10d508cc9a1b
md"""
method: $(@bind method Select([keys(diffeqmethods)...]))
"""

# ╔═╡ 3e1e62ec-c50a-499e-b516-8478904429c5
m = 2; n = 50;

# ╔═╡ 12ab322c-60ae-419f-9334-82f2f7ee7b59
t1 = @elapsed sol1, sys1, err1 = run_vfvm(m = m, n = n);history_summary(sol1)

# ╔═╡ 604898ba-1e8f-4c7c-9711-9958a8351854
t2 = @elapsed sol2, sys2, err2 = run_diffeq(m = m, n = n, solver = diffeqmethods[method]());history_summary(sol2)

# ╔═╡ 0676e28e-4e4e-4976-ab57-fb2d2e062625
let
    aspect = 600
    vis = GridVisualizer(layout = (1, 2), resolution = (650, 300))
    scalarplot!(vis[1, 1], sys1, sol1; aspect)
    scalarplot!(vis[1, 2], sys2, sol2; aspect)
    reveal(vis)
end

# ╔═╡ 2535425f-da87-4488-a834-f523a260bfde
md"""
Left: $(@sprintf("VoronoiFVM implicit Euler: %.0f ms e=%.2e",t1*1000,err1))

Right: $(@sprintf("    %s: %.0f ms, e=%.2e",method,t2*1000,err2))
"""

# ╔═╡ 84a29eb7-d936-4bd3-b15f-2886a4ca4985
@test err2 < err1


# ╔═╡ Cell order:
# ╠═6d29fa9a-54b9-490f-bb96-e44244522d27
# ╠═6d467640-b19c-4f77-845d-f9b4aca62104
# ╟─02424193-41e8-4cec-8f52-6fd66173ace8
# ╠═870b8b91-cd74-463e-b258-c092cd0af200
# ╠═78208d7b-9d79-415c-ab3a-85948251e635
# ╠═4ef024a4-cb1d-443d-97fb-ab3a32a78ffd
# ╠═9239409b-6de0-4157-8a35-412c909efa96
# ╟─2a8ac57e-486d-4825-95ab-f0402b910dbd
# ╠═12ab322c-60ae-419f-9334-82f2f7ee7b59
# ╟─3a004ab9-2705-4f5c-8e6e-10d508cc9a1b
# ╠═3e1e62ec-c50a-499e-b516-8478904429c5
# ╠═604898ba-1e8f-4c7c-9711-9958a8351854
# ╟─0676e28e-4e4e-4976-ab57-fb2d2e062625
# ╟─2535425f-da87-4488-a834-f523a260bfde
# ╠═84a29eb7-d936-4bd3-b15f-2886a4ca4985
