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

# ╔═╡ f7c86694-ca5c-48fa-a83c-88a7fb359396
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

# ╔═╡ c9bf8371-95ec-4a4f-a6dc-5c570bccd94a
md"""
## The wave equation as system of equations
"""


# ╔═╡ 8b8b8167-feab-4a66-8396-d991e7afb2d2
md"""
This is the n-dimensional wave equation:

```math
u_{tt}- c^2 \Delta u = 0
```


We can create a system of first order in time PDEs out of this:

```math 
    \begin{aligned}
        u_t - v&=0\\
     v_t -c^2\Delta u&=0
\end{aligned}
```

This allows for a quick implementation in VoronoiFVM (which may be not the optimal way, in particular with respect to time discretization).
"""


# ╔═╡ f4002117-f4c4-4191-9198-d75a2a2adc9a
const iu = 1; const iv = 2;

# ╔═╡ 335129d6-1da0-4d1f-ad14-02719e7fd215
storage(y, u, node, data) = y .= u;

# ╔═╡ 44f570e7-7834-4fae-a222-80e7ed29eb28
reaction(y, u, node, data) = y[iu] = -u[iv];

# ╔═╡ bee98dec-af80-422f-97ee-126154612708
flux(y, u, edge, data) = y[iv] = data.c^2 * (u[iu, 1] - u[iu, 2]);

# ╔═╡ 5b233195-e614-407f-824d-ca4c3cb9287c
md"""
Implementation of transparent or mirror bc
"""

# ╔═╡ 7b195fdb-9190-4eb4-8341-5a32947278d2
function brea(y, u, node, data)
    if node.region == 2
        if data.bctype == :transparent
            y[iu] = data.c * u[iu]
        elseif data.bctype == :mirror
            boundary_dirichlet!(y, u, node, species = 1, region = 2, value = 0)
        end
    end
    return nothing
end;


# ╔═╡ d00668ed-0837-408f-b235-9fc54267b97e
md"""
Wave velocity: 
"""

# ╔═╡ f8d02094-e1ce-4099-bc35-e940bbf81033
const c = 0.1

# ╔═╡ 9478f888-3722-4a1f-8981-e717b890ddf2
md"""
Domain length:
"""

# ╔═╡ 90617940-a390-4446-96c1-5308482c4acb
L = 4

# ╔═╡ ccf238c2-ad24-4f53-99a4-dc17d3ed249d
N = 151

# ╔═╡ dfe4b7d9-47b9-45b5-b382-b609947f0a2b
dt = 1.0e-2; tend = 100.0;

# ╔═╡ 2c3f949c-4d0d-4f54-a580-94e90a79ddb5
grid = simplexgrid(range(-L, L, length = N));

# ╔═╡ 8f5eac27-de80-4a02-b776-6205b8a805f8
md"""
Perturbation in the center of the domain: 
"""

# ╔═╡ 54a04c00-db2c-43ab-b3f4-aa86119ab4cb
md"""
Boundary condition at x=L: $(@bind bc2type Select(["reflection","mirror","transparent"]))
"""

# ╔═╡ 3cf0f943-f692-4031-81e2-a05d0bc5448b
sys = VoronoiFVM.System(grid, storage = storage, flux = flux, breaction = brea, reaction = reaction, data = (c = c, bctype = Symbol(bc2type)), species = [1, 2])

# ╔═╡ 5b33c8e6-f270-40bc-aa14-611804e1b265
if bc2type == "reflection"
    md"""Reflection (Neumann) bc ``\partial_x u|_{x=L}=0``"""
elseif bc2type == "mirror"
    md"""Mirror (Dirichlet) bc ``u|_{x=L}=0``"""
elseif bc2type == "transparent"
    md"""Transparent  bc ``(\partial_t u - c\partial_xu)|_{x=0}=0`` """
end

# ╔═╡ fc50e3fc-4dc9-4721-97a6-47d8f5ed8851
diffeqmethods = OrderedDict(
    "QNDF2" => QNDF2,
    "FBDF" => FBDF,
    "Rosenbrock23 (Rosenbrock)" => Rosenbrock23,
    "Implicit Euler" => ImplicitEuler,
    "Implicit Midpoint" => ImplicitMidpoint,
)

# ╔═╡ 46c6c09e-e610-4823-b0dc-fc6eaa557fc8
md"""
Package wave number κ: $(@bind κ Scrubbable(0:5,default=3)) method: $(@bind method Select([keys(diffeqmethods)...]))

t=$(@bind t  PlutoUI.Slider(range(0,tend,length=10001),default=tend/2;show_value=true)) 

"""

# ╔═╡ e65b01d0-d0bd-4ac3-b94a-d9d36b7c42c3
begin
    inival = unknowns(sys, inival = 0)
    inival[1, :] .= map(x -> cos(κ * π * x) * exp(-x^2 / 0.1), grid)
end;

# ╔═╡ 05bac3fc-92de-494f-8ae1-098155ea506f
problem = ODEProblem(sys, inival, (0.0, tend));

# ╔═╡ 18f636ab-645f-446b-aef0-66dbdcc50384
tsol = solve(
    problem, diffeqmethods[method]();
    force_dtmin = true,
    adaptive = true,
    reltol = 1.0e-4,
    abstol = 1.0e-5,
    dtmin = dt,
    progress = true,
    progress_steps = 1,
    dt = dt
);

# ╔═╡ cbcd3633-a376-4c34-9eea-ed94534ae8b3
tsol1 = reshape(tsol, sys);

# ╔═╡ 41ce301c-3e71-4f3f-80b8-f6e115b3b5e0
let
    vis = GridVisualizer()
    scalarplot!(vis, sys, tsol1, colormap = :bwr, limits = (-1, 1), levels = (-0.9:0.2:0.9))
    reveal(vis)
end

# ╔═╡ 7e5140d7-3304-4eac-be4c-981324dc346b
let
    u = tsol1(t)
    scalarplot(grid, u[1, :], flimits = (-1, 1), clear = true, show = true, title = "t=$(t)", resolution = (600, 150))
end


# ╔═╡ Cell order:
# ╠═f7c86694-ca5c-48fa-a83c-88a7fb359396
# ╠═6d467640-b19c-4f77-845d-f9b4aca62104
# ╟─c9bf8371-95ec-4a4f-a6dc-5c570bccd94a
# ╟─8b8b8167-feab-4a66-8396-d991e7afb2d2
# ╠═f4002117-f4c4-4191-9198-d75a2a2adc9a
# ╠═335129d6-1da0-4d1f-ad14-02719e7fd215
# ╠═44f570e7-7834-4fae-a222-80e7ed29eb28
# ╠═bee98dec-af80-422f-97ee-126154612708
# ╟─5b233195-e614-407f-824d-ca4c3cb9287c
# ╠═7b195fdb-9190-4eb4-8341-5a32947278d2
# ╟─d00668ed-0837-408f-b235-9fc54267b97e
# ╠═f8d02094-e1ce-4099-bc35-e940bbf81033
# ╟─9478f888-3722-4a1f-8981-e717b890ddf2
# ╠═90617940-a390-4446-96c1-5308482c4acb
# ╠═ccf238c2-ad24-4f53-99a4-dc17d3ed249d
# ╠═dfe4b7d9-47b9-45b5-b382-b609947f0a2b
# ╠═2c3f949c-4d0d-4f54-a580-94e90a79ddb5
# ╠═3cf0f943-f692-4031-81e2-a05d0bc5448b
# ╟─8f5eac27-de80-4a02-b776-6205b8a805f8
# ╠═e65b01d0-d0bd-4ac3-b94a-d9d36b7c42c3
# ╠═05bac3fc-92de-494f-8ae1-098155ea506f
# ╠═18f636ab-645f-446b-aef0-66dbdcc50384
# ╟─54a04c00-db2c-43ab-b3f4-aa86119ab4cb
# ╟─5b33c8e6-f270-40bc-aa14-611804e1b265
# ╟─46c6c09e-e610-4823-b0dc-fc6eaa557fc8
# ╠═7e5140d7-3304-4eac-be4c-981324dc346b
# ╠═41ce301c-3e71-4f3f-80b8-f6e115b3b5e0
# ╠═cbcd3633-a376-4c34-9eea-ed94534ae8b3
# ╠═fc50e3fc-4dc9-4721-97a6-47d8f5ed8851
