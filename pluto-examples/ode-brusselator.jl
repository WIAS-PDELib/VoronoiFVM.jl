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

# ╔═╡ efb2a070-80b6-4c97-b5a3-6aa5aebc1706
# Convert this cell to markdown in order to enable Pluto's inbuilt package manager
if isdefined(Main, :PlutoRunner)
    using Pkg
    docsdir = joinpath(@__DIR__, "..", "docs")
    if isdir(docsdir)
        Pkg.activate(docsdir)
    end
end


# ╔═╡ 7b4b635d-c600-4d6a-9cef-64bd97eb6a3e
begin
    using Test
    using Revise
    using Printf
    using VoronoiFVM
    using OrdinaryDiffEqBDF
    using OrdinaryDiffEqRosenbrock
    using OrdinaryDiffEqSDIRK
    using LinearAlgebra
    using PlutoUI
    using ExtendableGrids
    using DataStructures
    using GridVisualize
    if isdefined(Main, :PlutoRunner)
        using CairoMakie
        default_plotter!(CairoMakie)
        CairoMakie.activate!(type = "png")
    end
end

# ╔═╡ adea3b1c-854f-4f9d-9013-ab6a2a6f7fd7
md"""
##  A Brusselator problem
"""

# ╔═╡ 9cefef6c-b487-4daa-a814-efad85f6a634
md"""
Two species diffusing and interacting via a reaction
```math
\begin{aligned}
  \partial_t u_1 - \nabla \cdot (D_1 \nabla u_1) &+ (B+1)u_1-A-u_1^2u_2  =0\\
  \partial_t u_2 - \nabla \cdot (D_2 \nabla u_2) &+ u_1^2u_2 -B u_1 =0\\
\end{aligned}
```
"""


# ╔═╡ 32107aac-050d-4f50-b95c-383e1bb38652
begin
    const bruss_A = 2.25
    const bruss_B = 7.0
    const bruss_D_1 = 0.025
    const bruss_D_2 = 0.25
    const pert = 0.1
    const bruss_tend = 150
end;


# ╔═╡ 2fb1c53a-7ff8-4ac9-ae78-83bcbc57c926
function bruss_storage(f, u, node, data)
    f[1] = u[1]
    f[2] = u[2]
    return nothing
end;

# ╔═╡ 71b7e770-5cd4-4671-a76a-8e29eda04eec
function bruss_diffusion(f, u, edge, data)
    f[1] = bruss_D_1 * (u[1, 1] - u[1, 2])
    f[2] = bruss_D_2 * (u[2, 1] - u[2, 2])
    return nothing
end;

# ╔═╡ f1e7a242-d631-4624-b0c0-ac44f139d77c
function bruss_reaction(f, u, node, data)
    f[1] = (bruss_B + 1.0) * u[1] - bruss_A - u[1]^2 * u[2]
    f[2] = u[1]^2 * u[2] - bruss_B * u[1]
    return nothing
end;

# ╔═╡ 7e214c83-9c5c-40a9-8b00-db79dfec9a88
diffeqmethods = OrderedDict(
    "Rosenbrock23 (Rosenbrock)" => Rosenbrock23,
    "QNDF2 (Like matlab's ode15s)" => QNDF2,
    "FBDF" => FBDF,
    "Implicit Euler" => ImplicitEuler
)

# ╔═╡ 95be1da7-5f98-4a15-bd8e-7db1ee324768
begin

    function ODESolver(system, inival, solver)
        state = VoronoiFVM.SystemState(system)
        problem = ODEProblem(state, inival, (0, bruss_tend))
        odesol = solve(
            problem,
            solver,
            dt = 1.0e-5, reltol = 1.0e-4
        )
        return reshape(odesol, system; state)
    end

    sys0 = VoronoiFVM.System(simplexgrid(0:0.1:1), species = [1, 2], flux = bruss_diffusion, storage = bruss_storage, reaction = bruss_reaction)
    problem0 = ODEProblem(sys0, unknowns(sys0), (0, 0.1))

    for method in diffeqmethods
        solve(problem0, method.second()) #precompile
    end
end

# ╔═╡ 1462d783-93d3-4ad4-8701-90bde88c7553
md"""
dim:$(@bind bruss_dim Scrubbable(1:2,default=1)) ``\;``  method: $(@bind bruss_method Select([keys(diffeqmethods)...])) ``\;`` t: $(@bind t_bruss PlutoUI.Slider(0:bruss_tend/200:bruss_tend,show_value=true,default=bruss_tend))
"""

# ╔═╡ d48ad585-9d0a-4b7e-a54b-3c76d8a5ca21
if bruss_dim == 1
    bruss_X = -1:0.01:1
    bruss_grid = simplexgrid(bruss_X)
else
    bruss_X = -1:0.1:1
    bruss_grid = simplexgrid(bruss_X, bruss_X)
end;

# ╔═╡ 719a15e1-7a69-4e70-b20e-d75fa448458e
bruss_system = VoronoiFVM.System(
    bruss_grid, species = [1, 2],
    flux = bruss_diffusion, storage = bruss_storage, reaction = bruss_reaction
);

# ╔═╡ 62a1dad1-b095-4df9-b1f8-e6a97084d8f8
begin
    inival = unknowns(bruss_system, inival = 0)
    coord = bruss_grid[Coordinates]
    fpeak(x) = exp(-norm(10 * x)^2)
    for i in 1:size(inival, 2)
        @views inival[1, i] = fpeak(coord[:, i])
        @views inival[2, i] = fpeak(coord[:, i])
    end
end

# ╔═╡ e71a2ed0-5f39-473f-87a0-6f61748f2793
t_run = @elapsed bruss_tsol = ODESolver(bruss_system, inival, diffeqmethods[bruss_method]());

# ╔═╡ c1da7d8e-2921-4366-91f0-dc8e1834595b
(t_run = t_run, VoronoiFVM.history_details(bruss_tsol)...)

# ╔═╡ e7a8aae1-8e7a-4b7d-8ce6-701ea586b89a
let
    bruss_sol = bruss_tsol(t_bruss)

    vis = GridVisualizer(; layout = (1, 2), size = (600, 300))
    scalarplot!(vis[1, 1], bruss_grid, bruss_sol[1, :], limits = (0, 10), show = true, colormap = :summer)
    scalarplot!(vis[1, 2], bruss_grid, bruss_sol[2, :], limits = (0.5, 3), show = true, colormap = :summer)
end


# ╔═╡ 4a4a2f78-4d4c-4b9b-883e-1e57de41c9a9
scalarplot(bruss_tsol.t[1:(end - 1)], bruss_tsol.t[2:end] - bruss_tsol.t[1:(end - 1)], yscale = :log, resolution = (500, 200), xlabel = "t", ylabel = "Δt", title = "timesteps")


# ╔═╡ Cell order:
# ╠═efb2a070-80b6-4c97-b5a3-6aa5aebc1706
# ╠═7b4b635d-c600-4d6a-9cef-64bd97eb6a3e
# ╟─adea3b1c-854f-4f9d-9013-ab6a2a6f7fd7
# ╟─9cefef6c-b487-4daa-a814-efad85f6a634
# ╠═32107aac-050d-4f50-b95c-383e1bb38652
# ╠═2fb1c53a-7ff8-4ac9-ae78-83bcbc57c926
# ╠═71b7e770-5cd4-4671-a76a-8e29eda04eec
# ╠═f1e7a242-d631-4624-b0c0-ac44f139d77c
# ╠═95be1da7-5f98-4a15-bd8e-7db1ee324768
# ╟─7e214c83-9c5c-40a9-8b00-db79dfec9a88
# ╠═d48ad585-9d0a-4b7e-a54b-3c76d8a5ca21
# ╠═719a15e1-7a69-4e70-b20e-d75fa448458e
# ╠═62a1dad1-b095-4df9-b1f8-e6a97084d8f8
# ╠═e71a2ed0-5f39-473f-87a0-6f61748f2793
# ╠═c1da7d8e-2921-4366-91f0-dc8e1834595b
# ╟─1462d783-93d3-4ad4-8701-90bde88c7553
# ╟─e7a8aae1-8e7a-4b7d-8ce6-701ea586b89a
# ╟─4a4a2f78-4d4c-4b9b-883e-1e57de41c9a9
