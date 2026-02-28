### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 629854c8-966c-4ef5-aa8b-3ab8d85c753f
# Convert this cell to markdown in order to enable Pluto's inbuilt package manager
if isdefined(Main, :PlutoRunner)
    using Pkg
    docsdir = joinpath(@__DIR__, "..", "docs")
    if isdir(docsdir)
        Pkg.activate(docsdir)
    end
end

# ╔═╡ b285aca3-dee5-4b77-9276-537563e8643b
begin
    using Revise
    using VoronoiFVM
    using ExtendableGrids
    using Test
    using PlutoUI
    using LinearAlgebra
    using GridVisualize
    if isdefined(Main, :PlutoRunner)
        using CairoMakie
        CairoMakie.activate!(; type = "png")
        GridVisualize.default_plotter!(CairoMakie)
    end
end;

# ╔═╡ 4ed0c302-26e4-468a-a40d-0e6406f802d0
md"""
# Nonlinear solver control
[Source](https://github.com/WIAS-PDELib/VoronoiFVM.jl/blob/master/pluto-examples/nonlinear-solvers.jl)

"""

# ╔═╡ eb9ea477-6122-4774-b4a4-04dd7346e2b6
md"""
Generally, nonlinear systems in this package  are solved using Newton's method.  In many cases, the default settings provided by this package work well. However, the convergence of Newton's method is only guaranteed with initial values s7ufficiently close to the exact solution. This notebook describes how change the default settings for the solution of nonlinear problems with VoronoiFVM.jl. 
"""

# ╔═╡ 7a104243-d3b9-421a-b494-5607c494b106
TableOfContents(; aside = false)

# ╔═╡ 9843b65c-6ca8-4ef8-a896-2c8cec4bff7c
md"""
Define a nonlinear Poisson equation to have an example. Let ``Ω=(0,10)`` and define
```math
\begin{aligned}
-Δ u + e^u-e^{-u} & = 0 & \text{in}\; Ω \\
	u(0)&=100\\
    u(10)&=0
\end{aligned}
```
"""

# ╔═╡ a70b8d85-66aa-4da8-8157-dd0244e3e4f6
X = 0:0.001:1

# ╔═╡ c8eda836-d719-4412-895e-c3a24fec21ec
flux(y, u, edge, data) = y[1] = u[1, 1] - u[1, 2];

# ╔═╡ c09f5dfc-fc47-4952-8051-54731ec2b00b
function reaction(y, u, node, data)
    eplus = exp(u[1])
    eminus = 1 / eplus
    y[1] = eplus - eminus
    return nothing
end

# ╔═╡ eab04557-5084-4174-b275-b4d4399238e5
function bc(y, u, node, data)
    boundary_dirichlet!(y, u, node; region = 1, value = 100)
    boundary_dirichlet!(y, u, node; region = 2, value = 0.0)
    return nothing
end;

# ╔═╡ 316112fd-6553-494a-8e4a-65b34829891d
system = VoronoiFVM.System(X; flux = flux, reaction = reaction, bcondition = bc, species = 1);

# ╔═╡ 42e54ff9-fc11-4a31-ba10-32f8d817d50c
md"""
## Solution using default settings
"""

# ╔═╡ 050ed807-1bca-4015-85f3-d7347ecb7e6b
begin
    sol = solve(system; log = true)
    hist = history(sol)
end;

# ╔═╡ b9bb8020-5470-4964-818c-7f9b3bf2a8b4
scalarplot(
    system,
    sol;
    resolution = (500, 200),
    xlabel = "x",
    ylable = "y",
    title = "solution"
)

# ╔═╡ b3124c06-1f40-46f5-abee-9c2e8e538162
md"""
With `log=true`, the `solve` method in addition to the solution records the solution
history which after finished solution can be obtatined as `history(sol)`.
"""

# ╔═╡ 973db266-eb91-46e8-a917-9beeeb2c1ea7
md"""
The history can be plotted:
"""

# ╔═╡ 20e925f3-43fa-4db1-a656-79cf9c1c3303
function plothistory(h)
    return scalarplot(
        1:length(h),
        h;
        resolution = (500, 200),
        yscale = :log,
        xlabel = "step",
        ylabel = "||δu||_∞",
        title = "Maximum norm of Newton update"
    )
end;

# ╔═╡ ebdc2c82-f72e-4e35-a63f-4ba5154e294f
plothistory(hist)

# ╔═╡ 3d49aafd-79a6-4e2f-b422-24a5be7aa87a
md"""
History can be summarized:
"""

# ╔═╡ a217a308-5569-41e4-9d9d-418217017030
summary(hist)

# ╔═╡ f951b78c-a3fb-432e-bf6e-b956049d6a0d
md"""
History can be explored in detail:
"""

# ╔═╡ fcd7beb6-51b2-4ca8-a184-43ba6b5d2c1a
VoronoiFVM.details(hist)

# ╔═╡ baed6e43-747b-4557-95c3-d4805f12b3a1
md"""
With default solver settings, for this particular problem, Newton's method needs $(length(history(sol))) iteration steps.
"""

# ╔═╡ ccef0590-d5f8-4ee2-bb7a-d48ccfbd4d99
check(sol) = isapprox(sum(sol), 2554.7106586964906; rtol = 1.0e-12)

# ╔═╡ c0432a54-85ec-4478-bd75-f5b43770a117
@test check(sol)

# ╔═╡ a4c1a2f5-dddb-4a45-bc03-e168cbd7d569
md"""
## Damping 
"""

# ╔═╡ 38539474-af65-4f0f-9aa1-2292f4f6331c
md"""
Try to use a damped version of Newton method. The damping scheme is rather simple: an initial damping value `damp_initial` is increased by a growth factor `damp_growth` in each iteration until it reaches 1.
"""

# ╔═╡ d961d026-0b55-46c2-8555-8ef0763d8016
begin
    sol1 = solve(system; log = true, inival = 1, damp_initial = 0.15, damp_growth = 1.5)
    hist1 = history(sol1)
end

# ╔═╡ e66d20f0-4987-471b-82ee-0c56160f9b01
plothistory(hist1)

# ╔═╡ 35971019-fa07-4033-aebf-7872030a0cef
VoronoiFVM.details(hist1)

# ╔═╡ 63ce84fc-e81b-4768-8122-36bfbd789727
summary(hist1)

# ╔═╡ bf3fe305-eecc-413e-956c-9737b9160f83
md"""
We see that the number of iterations decreased significantly.
"""

# ╔═╡ c8227ea2-2189-438f-bd02-e9d803031830
@test check(sol1)

# ╔═╡ f8cf5cdb-647c-4eb8-8bde-b12844f72b24
md"""
## Embedding
"""

# ╔═╡ eea852d0-937e-4af0-8688-f941e5d31697
md"""
Another possibility is the embedding (or homotopy) via a parameter: start with solving a simple problem and increase the level of complexity by increasing the parameter until the full problem is solved. This process is controlled by the parameters 
- `Δp`: initial parameter step size
- `Δp_min`: minimal parameter step size
- `Δp_max`: maximum parameter step size
- `Δp_grow`: maximum growth factor
- `Δu_opt`: optimal difference of solutions between two embedding steps

After successful solution of a parameter, the new parameter step size is calculated as
`Δp_new=min(Δp_max, Δp_grow, Δp*Δu_opt/(|u-u_old|+1.0e-14))` and adjusted to the end of the parameter interval. 

If the solution is unsuccessful, the parameter stepsize is halved and solution is retried, until the minimum step size is reached. 
"""

# ╔═╡ a71cbcd4-310e-47a8-94f9-1159995a7711
function pbc(y, u, node, data)
    boundary_dirichlet!(y, u, node; region = 1, value = 100 * embedparam(node))
    boundary_dirichlet!(y, u, node; region = 2, value = 0)
    return nothing
end;

# ╔═╡ 89435c65-0520-4430-8727-9d013df6182d
system2 = VoronoiFVM.System(
    X;
    flux = flux,
    reaction = function (y, u, node, data)
        reaction(y, u, node, data)

        y[1] = y[1] * embedparam(node)
        return nothing
    end,
    bcondition = pbc,
    species = 1,
);

# ╔═╡ cb382145-c4f1-4222-aed7-32fa1e3bd7e4
begin
    sol2 = solve(
        system2;
        inival = 0,
        log = true,
        embed = (0, 1),
        Δp = 0.1,
        Δp_grow = 1.2,
        Δu_opt = 15
    )
    history2 = history(sol2)
end

# ╔═╡ a0b2aaf5-f7b1-40eb-ac4e-9790a8bbf09d
summary(history2)

# ╔═╡ 0d2311aa-f79a-4a44-bac4-59d6c5457ca5
plothistory(vcat(history2[2:end]...))

# ╔═╡ 310675d5-d175-43fc-bb86-3fb1c2d3c24c
sol2.u[end]

# ╔═╡ 75ab714d-0251-42ef-a415-ba6bed4c688f
@test check(sol2.u[end])

# ╔═╡ fe22e2d8-fd70-49e2-8baf-d5ec3faead24
md"""
For this particular problem, embedding uses less overall Newton steps than the default settings, but the damped method is faster.
"""

# ╔═╡ 2cb7f0f4-7635-462c-89f8-93d18e7f24fb
md"""
## SolverControl
"""

# ╔═╡ 5f13831d-b73f-44fb-a870-16261f926ed5
md"""
Here we show the docsctring of `SolverControl` (formerly `NewtonControl`). This is a struct which can be passed to the `solve` method. Alternatively, as shown in this notebook, keyword arguments named like its entries can be passed directly to the `solve` method.
"""

# ╔═╡ 038d096a-1339-403c-aa9c-3112442d622d
@doc VoronoiFVM.SolverControl

# ╔═╡ fc0245fe-1bf2-45a3-aa7c-9cce8d7eef37
html"""<hr>"""

# ╔═╡ bdbe6513-70b1-4d97-a79c-71534caad2b7
html"""<hr>"""


# ╔═╡ Cell order:
# ╟─4ed0c302-26e4-468a-a40d-0e6406f802d0
# ╟─eb9ea477-6122-4774-b4a4-04dd7346e2b6
# ╟─7a104243-d3b9-421a-b494-5607c494b106
# ╠═629854c8-966c-4ef5-aa8b-3ab8d85c753f
# ╠═b285aca3-dee5-4b77-9276-537563e8643b
# ╟─9843b65c-6ca8-4ef8-a896-2c8cec4bff7c
# ╠═a70b8d85-66aa-4da8-8157-dd0244e3e4f6
# ╠═c8eda836-d719-4412-895e-c3a24fec21ec
# ╠═c09f5dfc-fc47-4952-8051-54731ec2b00b
# ╠═eab04557-5084-4174-b275-b4d4399238e5
# ╠═316112fd-6553-494a-8e4a-65b34829891d
# ╟─42e54ff9-fc11-4a31-ba10-32f8d817d50c
# ╠═050ed807-1bca-4015-85f3-d7347ecb7e6b
# ╠═b9bb8020-5470-4964-818c-7f9b3bf2a8b4
# ╟─b3124c06-1f40-46f5-abee-9c2e8e538162
# ╟─973db266-eb91-46e8-a917-9beeeb2c1ea7
# ╠═20e925f3-43fa-4db1-a656-79cf9c1c3303
# ╠═ebdc2c82-f72e-4e35-a63f-4ba5154e294f
# ╟─3d49aafd-79a6-4e2f-b422-24a5be7aa87a
# ╠═a217a308-5569-41e4-9d9d-418217017030
# ╟─f951b78c-a3fb-432e-bf6e-b956049d6a0d
# ╟─fcd7beb6-51b2-4ca8-a184-43ba6b5d2c1a
# ╟─baed6e43-747b-4557-95c3-d4805f12b3a1
# ╠═ccef0590-d5f8-4ee2-bb7a-d48ccfbd4d99
# ╠═c0432a54-85ec-4478-bd75-f5b43770a117
# ╟─a4c1a2f5-dddb-4a45-bc03-e168cbd7d569
# ╟─38539474-af65-4f0f-9aa1-2292f4f6331c
# ╠═d961d026-0b55-46c2-8555-8ef0763d8016
# ╟─e66d20f0-4987-471b-82ee-0c56160f9b01
# ╠═35971019-fa07-4033-aebf-7872030a0cef
# ╠═63ce84fc-e81b-4768-8122-36bfbd789727
# ╟─bf3fe305-eecc-413e-956c-9737b9160f83
# ╠═c8227ea2-2189-438f-bd02-e9d803031830
# ╟─f8cf5cdb-647c-4eb8-8bde-b12844f72b24
# ╟─eea852d0-937e-4af0-8688-f941e5d31697
# ╠═a71cbcd4-310e-47a8-94f9-1159995a7711
# ╠═89435c65-0520-4430-8727-9d013df6182d
# ╠═cb382145-c4f1-4222-aed7-32fa1e3bd7e4
# ╠═a0b2aaf5-f7b1-40eb-ac4e-9790a8bbf09d
# ╠═0d2311aa-f79a-4a44-bac4-59d6c5457ca5
# ╠═310675d5-d175-43fc-bb86-3fb1c2d3c24c
# ╠═75ab714d-0251-42ef-a415-ba6bed4c688f
# ╟─fe22e2d8-fd70-49e2-8baf-d5ec3faead24
# ╟─2cb7f0f4-7635-462c-89f8-93d18e7f24fb
# ╟─5f13831d-b73f-44fb-a870-16261f926ed5
# ╠═038d096a-1339-403c-aa9c-3112442d622d
# ╟─fc0245fe-1bf2-45a3-aa7c-9cce8d7eef37
# ╟─bdbe6513-70b1-4d97-a79c-71534caad2b7
