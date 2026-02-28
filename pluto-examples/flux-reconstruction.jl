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

# ╔═╡ f154cd22-4843-45b5-b735-16c1b5b27b19
# Convert this cell to markdown in order to enable Pluto's inbuilt package manager
if isdefined(Main, :PlutoRunner)
    using Pkg
    docsdir = joinpath(@__DIR__, "..", "docs")
    if isdir(docsdir)
        Pkg.activate(docsdir)
    end
end


# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin
    using Revise
    using Test
    using SimplexGridFactory, Triangulate, ExtendableGrids, VoronoiFVM
    using PlutoUI, GridVisualize
    if isdefined(Main, :PlutoRunner)
        using CairoMakie
        CairoMakie.activate!(; type = "png", visible = false)
        GridVisualize.default_plotter!(CairoMakie)
    end
end;

# ╔═╡ de468cb9-b34d-4d2e-b911-9b93920caca1
md"""
# Flux reconstruction and visualization for the Laplace operator
"""

# ╔═╡ 531de061-d943-4b5a-85f2-cbd48bb049ce
md"""
We demonstrate the reconstruction of the gradient vector field from the solution of the Laplace operator and its visualization.
"""

# ╔═╡ 184193b6-39ef-4d0c-92a3-157fa5809832
TableOfContents(; aside = false)

# ╔═╡ 30dc968f-44df-45ea-bdb3-c938a8026224
md"""
## Grid
"""

# ╔═╡ 4064355c-6367-4797-96f6-020bca9b854e
md"""
Define a "Swiss cheese domain" with punched-out holes, where each hole boundary corresponds to a different boundary condition.
"""

# ╔═╡ 928a70c5-4706-40a1-9387-abcb71c09443
function swiss_cheese_2d()
    function circlehole!(builder, center, radius; n = 20)
        points = [
            point!(builder, center[1] + radius * sin(t), center[2] + radius * cos(t))
                for t in range(0, 2π; length = n)
        ]
        for i in 1:(n - 1)
            facet!(builder, points[i], points[i + 1])
        end
        facet!(builder, points[end], points[1])
        return holepoint!(builder, center)
    end

    builder = SimplexGridBuilder(; Generator = Triangulate)
    cellregion!(builder, 1)
    maxvolume!(builder, 0.1)
    regionpoint!(builder, 0.1, 0.1)

    p1 = point!(builder, 0, 0)
    p2 = point!(builder, 10, 0)
    p3 = point!(builder, 10, 10)
    p4 = point!(builder, 0, 10)

    facetregion!(builder, 1)
    facet!(builder, p1, p2)
    facet!(builder, p2, p3)
    facet!(builder, p3, p4)
    facet!(builder, p4, p1)

    holes = [
        1.0 2.0
        8.0 9.0
        2.0 8.0
        8.0 4.0
        9.0 1.0
        3.0 4.0
        4.0 6.0
        7.0 9.0
        4.0 7.0
        7.0 5.0
        2.0 1.0
        4.0 1.0
        4.0 8.0
        3.0 6.0
        4.0 9.0
        6.0 9.0
        3.0 5.0
        1.0 4.0
    ]'

    radii = [
        0.15,
        0.15,
        0.1,
        0.35,
        0.2,
        0.3,
        0.1,
        0.4,
        0.1,
        0.4,
        0.2,
        0.2,
        0.2,
        0.35,
        0.15,
        0.25,
        0.15,
        0.25,
    ]

    for i in 1:length(radii)
        facetregion!(builder, i + 1)
        circlehole!(builder, holes[:, i], radii[i])
    end

    return simplexgrid(builder)
end

# ╔═╡ bc304085-69c3-4974-beb4-6f2b981ac0f1
grid = swiss_cheese_2d()

# ╔═╡ ad0367e6-7308-46bc-a7ff-9c17dcfe470d
gridplot(grid)

# ╔═╡ f27d80dc-ae41-4889-aede-daa9848488d4
md"""
## System + solution
"""

# ╔═╡ 758ee64a-55f9-495c-ad0f-79ee6cb1c922
mutable struct Params
    val11::Float64
end

# ╔═╡ c71479bb-886c-443c-b8f6-9b7ad8b34cef
params = Params(5)

# ╔═╡ 62b300e4-8476-4531-8042-fdb199355fbb
md"""
Simple flux function for Laplace operator 
"""

# ╔═╡ 236d6cd5-c190-49c7-98fe-021fae224455
flux(y, u, edge, data) = y[1] = u[1, 1] - u[1, 2];

# ╔═╡ 6bd1c695-ff0d-46c8-bbe8-9a8e4eea9a23
md"""
At hole #11, the value will be bound to a slider defined below
"""

# ╔═╡ b3f92fa7-6510-49e0-8c0d-b6ef6e897bb3
function bc(y, u, bnode, data)
    boundary_dirichlet!(y, u, bnode; region = 2, value = 10.0)
    boundary_dirichlet!(y, u, bnode; region = 3, value = 0.0)
    boundary_dirichlet!(y, u, bnode; region = 11, value = data.val11)
    return nothing
end

# ╔═╡ f4ebe6ad-4e04-4f33-9a66-6bec977adf4d
md"""
Define a finite volume system with Dirichlet boundary conditions at some of the holes
"""

# ╔═╡ f3af20fb-c3cc-41af-9782-d40adf2371f7
system = VoronoiFVM.System(grid; flux = flux, species = 1, bcondition = bc, data = params)

# ╔═╡ d86c43f3-ec2f-4fae-88d2-1068603e7044
md"""
Solve, and trigger solution upon boundary value change
"""

# ╔═╡ 90bd1da0-c371-4889-a00f-a17a27463c88
md"""
## Flux reconstruction
"""

# ╔═╡ 41bd1230-c87c-47b0-8e58-67ad55609fd3
md"""
Reconstruct the node flux. It is a ``d\times n_{spec}\times n_{nodes}`` tensor.
`nf[:,ispec,:]` then is a vector function representing the flux density of species `ispec` in each node of the domain. This readily can be fed into `GridVisualize.vectorplot`.

The reconstruction is based on the  "magic formula"
R. Eymard, T. Gallouet, R. Herbin, IMA Journal of Numerical Analysis (2006)
26, 326−353 ([Arxive version](https://arxiv.org/abs/math/0505109v1)), Lemma 2.4 .

"""

# ╔═╡ 17be52fb-f55b-4b3d-85e5-33f36134046b
vis = GridVisualizer(; dim = 2, resolution = (400, 400))

# ╔═╡ 03f582ec-4e95-4ca4-8482-9c797027810d
md"""
``v_{11}:`` $(@bind  val11 PlutoUI.Slider(0:0.1:10,default=5,show_value=true))
"""

# ╔═╡ 27efdcac-8ee4-47e4-905d-8d8c7313ddd1
begin
    params.val11 = val11
    sol = solve(system)
end;

# ╔═╡ 18d5bc33-2578-41d0-a390-c164d754b8e1
@test params.val11 != 5.0 || isapprox(sum(sol), 7842.2173682050525; rtol = 1.0e-12)

# ╔═╡ 41f427c1-b6ad-46d4-9151-1d872b4efeb6
nf = nodeflux(system, sol)

# ╔═╡ 0e34c818-021b-44c9-8ee4-1a737c3de9cb
@test params.val11 != 5.0 || isapprox(sum(nf), 978.000534849034; rtol = 1.0e-14)

# ╔═╡ 9468db0c-e924-4737-9b75-6bec753aafa9
md"""
Joint plot of solution and flux reconstruction
"""

# ╔═╡ 531edb71-6d32-4231-b117-5e36416d2fb1
begin
    scalarplot!(vis, grid, sol[1, :]; levels = 9, colormap = :summer, clear = true)
    vectorplot!(vis, grid, nf[:, 1, :]; clear = false, vscale = 1.5)
    reveal(vis)
end

# ╔═╡ 7f2081ee-67af-4d6c-971a-a383b6377fdf
md"""
# The 1D case
"""

# ╔═╡ f7c139cd-df37-4e36-bff1-115c33bb9067
src(x) = exp(-x^2 / 0.01)

# ╔═╡ aec97834-1612-49ab-a19b-7fd74c83228f
source1d(y, node, data) = y[1] = src(node[1])

# ╔═╡ 387f52d0-926e-445d-9e3e-afebbc7207f6
flux1d(y, u, edge, data) = y[1] = u[1, 1]^2 - u[1, 2]^2

# ╔═╡ 3793696c-c934-4e56-a1e7-887fc2181970
function bc1d(y, u, bnode, data)
    boundary_dirichlet!(y, u, bnode; region = 1, value = 0.01)
    boundary_dirichlet!(y, u, bnode; region = 2, value = 0.01)
    return nothing
end

# ╔═╡ 159ffdb7-a5d9-45bd-a53f-ba3751c91ae5
grid1d = simplexgrid(-1:0.01:1)

# ╔═╡ 29257fc4-d94b-4cf1-8432-30ba3fc4dc1b
sys1d = VoronoiFVM.System(
    grid1d;
    flux = flux1d,
    bcondition = bc1d,
    source = source1d,
    species = [1],
)

# ╔═╡ d8df038e-9cfc-4eb4-9845-2244ac95190b
sol1d = solve(sys1d; inival = 0.1)

# ╔═╡ 58dedb6a-ab19-44b8-90a5-a7f67700bc2f
nf1d = nodeflux(sys1d, sol1d)

# ╔═╡ 14b9e972-2538-43f1-a558-c6495543c9db
let
    vis1d = GridVisualizer(; dim = 1, resolution = (500, 250), legend = :lt)
    scalarplot!(vis1d, grid1d, map(src, grid1d); label = "rhs", color = :blue)
    scalarplot!(vis1d, grid1d, sol1d[1, :]; label = "solution", color = :red, clear = false)
    vectorplot!(vis1d, grid1d, nf1d[:, 1, :]; label = "flux", clear = false, color = :green)
    reveal(vis1d)
end

# ╔═╡ c62f06a2-70d5-40e7-b137-3dc547f5e246
@test nf1d[1, 1, 101] ≈ 0.0

# ╔═╡ cdd80b71-50dc-4b64-b0ea-37c57d65012f
@test nf1d[1, 1, 1] ≈ -nf1d[1, 1, end]

# ╔═╡ d3e64470-8e11-4ef7-af34-453088910fee
html"""<hr>"""

# ╔═╡ 2b3132a4-5e99-47f2-b85f-c55e90f7f93f
html"""<hr>"""


# ╔═╡ Cell order:
# ╟─de468cb9-b34d-4d2e-b911-9b93920caca1
# ╟─531de061-d943-4b5a-85f2-cbd48bb049ce
# ╟─184193b6-39ef-4d0c-92a3-157fa5809832
# ╠═f154cd22-4843-45b5-b735-16c1b5b27b19
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╟─30dc968f-44df-45ea-bdb3-c938a8026224
# ╟─4064355c-6367-4797-96f6-020bca9b854e
# ╠═928a70c5-4706-40a1-9387-abcb71c09443
# ╠═bc304085-69c3-4974-beb4-6f2b981ac0f1
# ╟─ad0367e6-7308-46bc-a7ff-9c17dcfe470d
# ╟─f27d80dc-ae41-4889-aede-daa9848488d4
# ╠═758ee64a-55f9-495c-ad0f-79ee6cb1c922
# ╠═c71479bb-886c-443c-b8f6-9b7ad8b34cef
# ╟─62b300e4-8476-4531-8042-fdb199355fbb
# ╠═236d6cd5-c190-49c7-98fe-021fae224455
# ╟─6bd1c695-ff0d-46c8-bbe8-9a8e4eea9a23
# ╠═b3f92fa7-6510-49e0-8c0d-b6ef6e897bb3
# ╟─f4ebe6ad-4e04-4f33-9a66-6bec977adf4d
# ╠═f3af20fb-c3cc-41af-9782-d40adf2371f7
# ╟─d86c43f3-ec2f-4fae-88d2-1068603e7044
# ╠═27efdcac-8ee4-47e4-905d-8d8c7313ddd1
# ╠═18d5bc33-2578-41d0-a390-c164d754b8e1
# ╟─90bd1da0-c371-4889-a00f-a17a27463c88
# ╟─41bd1230-c87c-47b0-8e58-67ad55609fd3
# ╠═41f427c1-b6ad-46d4-9151-1d872b4efeb6
# ╠═0e34c818-021b-44c9-8ee4-1a737c3de9cb
# ╠═17be52fb-f55b-4b3d-85e5-33f36134046b
# ╠═03f582ec-4e95-4ca4-8482-9c797027810d
# ╟─9468db0c-e924-4737-9b75-6bec753aafa9
# ╠═531edb71-6d32-4231-b117-5e36416d2fb1
# ╟─7f2081ee-67af-4d6c-971a-a383b6377fdf
# ╠═f7c139cd-df37-4e36-bff1-115c33bb9067
# ╠═aec97834-1612-49ab-a19b-7fd74c83228f
# ╠═387f52d0-926e-445d-9e3e-afebbc7207f6
# ╠═3793696c-c934-4e56-a1e7-887fc2181970
# ╠═159ffdb7-a5d9-45bd-a53f-ba3751c91ae5
# ╠═29257fc4-d94b-4cf1-8432-30ba3fc4dc1b
# ╠═d8df038e-9cfc-4eb4-9845-2244ac95190b
# ╠═58dedb6a-ab19-44b8-90a5-a7f67700bc2f
# ╠═14b9e972-2538-43f1-a558-c6495543c9db
# ╠═c62f06a2-70d5-40e7-b137-3dc547f5e246
# ╠═cdd80b71-50dc-4b64-b0ea-37c57d65012f
# ╟─d3e64470-8e11-4ef7-af34-453088910fee
# ╟─2b3132a4-5e99-47f2-b85f-c55e90f7f93f
