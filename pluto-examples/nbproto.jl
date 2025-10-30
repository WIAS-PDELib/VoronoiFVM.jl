### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# ╔═╡ b285aca3-dee5-4b77-9276-537563e8643b
begin
    import Pkg as _Pkg
    haskey(ENV, "PLUTO_PROJECT") && _Pkg.activate(ENV["PLUTO_PROJECT"])
    using Revise
    using VoronoiFVM
    using ExtendableGrids
    using Test
    using PlutoUI
    using GridVisualize
    if isdefined(Main, :PlutoRunner)
        using CairoMakie
        CairoMakie.activate!(type = "png", visible = false)
        GridVisualize.default_plotter!(CairoMakie)
    end
end;

# ╔═╡ 4ed0c302-26e4-468a-a40d-0e6406f802d0
md"""
# Intro
"""

# ╔═╡ 7a104243-d3b9-421a-b494-5607c494b106
TableOfContents(; aside = false)

# ╔═╡ c8eda836-d719-4412-895e-c3a24fec21ec
scalarplot(sin.(0:0.1:10), size = (500, 200))

# ╔═╡ 3eef08af-f6ba-4874-82c0-65ff53e7f7da
@test 1 == 1

# ╔═╡ e7bb8e62-228b-4b80-824b-31ea22543fba
if isdefined(Main, :PlutoRunner)
    let figure = Figure()
        axis = Axis(figure[1, 1], aspect = DataAspect())
        X = 0:0.1:1
        Y = 0:0.1:1
        pts = [Point2f(x, y) for y in X, x in X]
        li = LinearIndices(pts)
        for j in 1:(length(Y) - 1)
            for i in 1:(length(X) - 1)
                t1 = [li[i, j], li[i + 1, j], li[i + 1, j + 1]]
                @views poly!(axis, pts[t1], strokecolor = :black, color = :red, strokewidth = 1)
                t2 = [li[i, j], li[i, j + 1], li[i + 1, j + 1]]
                @views poly!(axis, pts[t2], strokecolor = :black, color = :red, strokewidth = 1)
            end
        end
        figure
    end
end

# ╔═╡ fc0245fe-1bf2-45a3-aa7c-9cce8d7eef37
html"""<hr><hr><hr>"""

# ╔═╡ f50c6497-cba3-491a-bedd-5f94f88f76fb
md"""
# Appendix: Tests & Development
"""

# ╔═╡ ad899a81-baab-4433-8b7f-1e5c3b18dae6


# ╔═╡ bdbe6513-70b1-4d97-a79c-71534caad2b7


# ╔═╡ Cell order:
# ╠═b285aca3-dee5-4b77-9276-537563e8643b
# ╟─4ed0c302-26e4-468a-a40d-0e6406f802d0
# ╟─7a104243-d3b9-421a-b494-5607c494b106
# ╠═c8eda836-d719-4412-895e-c3a24fec21ec
# ╠═3eef08af-f6ba-4874-82c0-65ff53e7f7da
# ╠═e7bb8e62-228b-4b80-824b-31ea22543fba
# ╟─fc0245fe-1bf2-45a3-aa7c-9cce8d7eef37
# ╟─f50c6497-cba3-491a-bedd-5f94f88f76fb
# ╟─ad899a81-baab-4433-8b7f-1e5c3b18dae6
# ╟─bdbe6513-70b1-4d97-a79c-71534caad2b7
