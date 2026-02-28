### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 8e7c0b82-e038-4837-958f-fffeca16401a
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
    using CairoMakie
    using Revise
    using VoronoiFVM
    using ForwardDiff: derivative
end

# ╔═╡ 5e13b3db-570c-4159-939a-7e2268f0a102
md"""
# Bernoulli function test

We test the implementation of the Bernoulli function in VoronoiFVM against the evaluation
with BigFloat. This allows to optimize thresholds for switching between evaluation expressions.
"""

# ╔═╡ 3c0166c1-cad7-46b3-9330-2d8b39ce6774
md"""
### Reference with BigFLoat
"""

# ╔═╡ b47b781b-ec11-485a-9de6-061ad0957f46
function B_Big(x)
    bx = BigFloat(x)
    return Float64(bx / expm1(bx))
end

# ╔═╡ 33f6b010-e0c4-4c7b-9cbd-460b8ba0f606
function DB_Big(x)
    bx = BigFloat(x)
    bone = one(BigFloat)
    bex = exp(bx)
    b = -(x * bex - bex + bone) / ((bex - bone) * (bex - bone))
    return Float64(b)
end

# ╔═╡ a0e3bf41-e5ca-4395-a26e-941a8b897705
md"""
## Implementation using `expm1`
"""

# ╔═╡ c95c09fe-ac77-4032-a455-2a13b0d7edc2
B(x) = x / expm1(x)

# ╔═╡ c3fd0ff2-7111-4165-ad93-d6d7257301fa
md"""
## Approximation for small x

For small values of x, a Taylor approximation implemented using a Horner scheme is utilized, as the exponential expression runs into errors in the vicinity of zero and fails to evaluate at zero.. As as long as its error is large than that of the Taylor approximation calculated with the Taylor scheme, we should use the later one. 
"""

# ╔═╡ d8ba486e-44f3-414c-ba3b-e8f0f5e9a099
B(0.0)

# ╔═╡ 827fbfb1-2e3f-402a-a0c4-a4799df5838e
fbernoulli(0.0)

# ╔═╡ f99ff517-3afc-469b-af34-8fd59233a1df
B(nextfloat(0.0))

# ╔═╡ e78571b4-6980-447b-a71d-a9609fda57b0
fbernoulli(nextfloat(0.0))

# ╔═╡ cb10e915-b753-4cd5-9355-71558c83369c
derivative(B, 0.0)

# ╔═╡ a63959f9-ed7a-47dc-b253-78c763dd359f
derivative(fbernoulli, 0.0)

# ╔═╡ 16db7f90-7d23-4864-a059-38d09f3e5d3b
derivative(B, nextfloat(0.0))

# ╔═╡ 7678ad35-f29b-448c-bd85-247140d42456
derivative(fbernoulli, nextfloat(0.0))

# ╔═╡ 56ff3f5c-6fe9-4d44-a5ae-449c42efca62
smallX = collect(-0.5:(1.0e-4 + 1.0e-8):0.5);

# ╔═╡ 6e7c197b-8ad2-4b9d-a0bc-40a48db32387
let
    p = Figure(
        ; size = (600, 400),
    )
    ax = Axis(
        p[1, 1];
        title = "|B_Big(x)-B(x)|",
        xminorticksvisible = true,
        xminorgridvisible = true,
        xminorticks = IntervalsBetween(10),
        xticks = -0.5:0.1:0.5,
        yscale = log10,
        xlabel = "x",
        ylabel = "error"
    )

    lines!(
        ax,
        smallX,
        abs.(B_Big.(smallX) .- B.(smallX)) .+ 1.0e-20
    )
    ax2 = Axis(
        p[2, 1];
        title = "|B_Big(x)-VoronoiFVM.bernoulli_horner(x)|",
        xminorticksvisible = true,
        xminorgridvisible = true,
        xminorticks = IntervalsBetween(10),
        xticks = -0.5:0.1:0.5,
        yscale = log10,
        xlabel = "x",
        ylabel = "error"
    )
    lines!(
        ax2,
        smallX,
        abs.(B_Big.(smallX) .- VoronoiFVM.bernoulli_horner.(smallX)) .+ 1.0e-20
    )
    p
end

# ╔═╡ c035b09f-92ca-4af9-9bcc-c754650e7bb1
let
    p = Figure(
        ; size = (600, 400),
    )
    ax = Axis(
        p[1, 1];
        title = "|DB_Big(x)-derivative(B,x)|",
        xminorticksvisible = true,
        xminorgridvisible = true,
        xminorticks = IntervalsBetween(10),
        xticks = -0.5:0.1:0.5,
        yscale = log10,
        xlabel = "x",
        ylabel = "error"
    )

    lines!(
        ax,
        smallX,
        abs.(DB_Big.(smallX) .- derivative.(B, smallX)) .+ 1.0e-20
    )
    ax2 = Axis(
        p[2, 1];
        title = "|DB_Big(x)-derivative(bernoulli_horner,x)|",
        xminorticksvisible = true,
        xminorgridvisible = true,
        xminorticks = IntervalsBetween(10),
        xticks = -0.5:0.1:0.5,
        yscale = log10,
        xlabel = "x",
        ylabel = "error"
    )
    lines!(
        ax2,
        smallX,
        abs.(DB_Big.(smallX) .- derivative.(VoronoiFVM.bernoulli_horner, smallX)) .+ 1.0e-20
    )
    p
end


# ╔═╡ 5a293797-beb9-493e-af12-d978c50d6148
md"""
## Error comparison for  VoronoiFVM implementation
"""

# ╔═╡ 26cdb920-291a-4b54-963f-fd9bd610662f
largeX = -100:1.00001e-3:100;

# ╔═╡ feb21ce6-0ddc-45fb-90f4-1e46261a9110
let
    p = Figure(
        ; size = (650, 400)
    )
    Label(p[0, 1], "Positive Argument"; fontsize = 15, font = :bold)
    bf = B_Big.(smallX)
    vf = first.(fbernoulli_pm.(smallX))
    err = abs.(bf - vf) ./ bf
    maxerr = maximum(err)
    ax1 = Axis(
        p[1, 0:2];
        title = "Maximum relative error for small x: $maxerr",
        xminorticksvisible = true,
        xminorgridvisible = true,
        xminorticks = IntervalsBetween(10),
        xticks = -0.5:0.1:0.5,
        xlabel = "x",
        ylabel = "error"
    )

    lines!(
        ax1,
        smallX,
        err;
        color = :red
    )

    bf = B_Big.(largeX)
    vf = first.(fbernoulli_pm.(largeX))
    err = abs.(bf - vf)
    maxerr = maximum(err)

    ax2 = Axis(
        p[2, 0:2];
        title = "Maximum absolute error for large x: $maxerr",
        xminorticksvisible = true,
        xminorgridvisible = true,
        xminorticks = IntervalsBetween(4),
        xticks = -100:20:100,
        xlabel = "x",
        ylabel = "error"
    )

    lines!(
        ax2,
        largeX,
        err;
        color = :red
    )
    save("bernoulli_posarg.png", p)
    p
end

# ╔═╡ 3fda3dd7-1603-463c-9be4-03a217ada56f
let
    p = Figure(
        ; size = (650, 400)
    )
    Label(p[0, 1], "Negative Argument"; fontsize = 15, font = :bold)

    bf = B_Big.(-smallX)
    vf = last.(fbernoulli_pm.(smallX))
    err = abs.(bf - vf) ./ bf
    maxerr = maximum(err)

    ax1 = Axis(
        p[1, 0:2];
        title = "Maximum relative error for small x: $(maxerr)",
        xminorticksvisible = true,
        xminorgridvisible = true,
        xminorticks = IntervalsBetween(10),
        xticks = -0.5:0.1:0.5,
        xlabel = "x",
        ylabel = "error"
    )
    lines!(
        ax1,
        smallX,
        err;
        color = :red
    )

    bf = B_Big.(-largeX)
    vf = last.(fbernoulli_pm.(largeX))
    err = abs.(bf - vf)
    maxerr = maximum(err)
    ax2 = Axis(
        p[2, 0:2];
        title = "Maximum absolute error for large x: $maxerr",
        xminorticksvisible = true,
        xminorgridvisible = true,
        xminorticks = IntervalsBetween(4),
        xticks = -100:20:100,
        xlabel = "x",
        ylabel = "error"
    )

    lines!(
        ax2,
        largeX,
        err;
        color = :red
    )
    save("bernoulli_negarg.png", p)
    p
end

# ╔═╡ d3717c34-9fbf-470d-99f0-dabc3cc023fc
md"""
Derivative error
"""

# ╔═╡ 68481e52-31d8-4208-b687-f8002ad27232
let
    p = Figure(
        ; size = (650, 400)
    )

    Label(p[0, 1], "Derivative"; fontsize = 15, font = :bold)
    bf = DB_Big.(smallX)
    vf = derivative.(fbernoulli, smallX)
    err = abs.(bf - vf) ./ abs.(bf)
    maxerr = maximum(err)

    ax1 = Axis(
        p[1, 0:2];
        title = "Maximum relative error for small x: $maxerr",
        xminorticksvisible = true,
        xminorgridvisible = true,
        xminorticks = IntervalsBetween(10),
        xticks = -0.5:0.1:0.5, xlabel = "x",
        ylabel = "error"
    )

    lines!(
        ax1,
        smallX,
        err;
        color = :red
    )

    bf = DB_Big.(largeX)
    vf = derivative.(fbernoulli, largeX)
    err = abs.(bf - vf)
    maxerr = maximum(err)
    ax2 = Axis(
        p[2, 0:2];
        title = "Maximum absolute error for large x: $maxerr",
        xminorticksvisible = true,
        xminorgridvisible = true,
        xminorticks = IntervalsBetween(4),
        xticks = -100:20:100,
        xlabel = "x",
        ylabel = "error"
    )

    lines!(
        ax2,
        largeX,
        err;
        color = :red
    )
    save("bernoulli_derivative.png", p)
    p
end

# ╔═╡ 12f268d2-baa6-4c1b-bed5-e9df53b469fc
html"<hr>"

# ╔═╡ Cell order:
# ╠═8e7c0b82-e038-4837-958f-fffeca16401a
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╟─5e13b3db-570c-4159-939a-7e2268f0a102
# ╟─3c0166c1-cad7-46b3-9330-2d8b39ce6774
# ╠═b47b781b-ec11-485a-9de6-061ad0957f46
# ╠═33f6b010-e0c4-4c7b-9cbd-460b8ba0f606
# ╟─a0e3bf41-e5ca-4395-a26e-941a8b897705
# ╠═c95c09fe-ac77-4032-a455-2a13b0d7edc2
# ╟─c3fd0ff2-7111-4165-ad93-d6d7257301fa
# ╠═d8ba486e-44f3-414c-ba3b-e8f0f5e9a099
# ╠═827fbfb1-2e3f-402a-a0c4-a4799df5838e
# ╠═f99ff517-3afc-469b-af34-8fd59233a1df
# ╠═e78571b4-6980-447b-a71d-a9609fda57b0
# ╠═cb10e915-b753-4cd5-9355-71558c83369c
# ╠═a63959f9-ed7a-47dc-b253-78c763dd359f
# ╠═16db7f90-7d23-4864-a059-38d09f3e5d3b
# ╠═7678ad35-f29b-448c-bd85-247140d42456
# ╠═56ff3f5c-6fe9-4d44-a5ae-449c42efca62
# ╟─6e7c197b-8ad2-4b9d-a0bc-40a48db32387
# ╟─c035b09f-92ca-4af9-9bcc-c754650e7bb1
# ╟─5a293797-beb9-493e-af12-d978c50d6148
# ╠═26cdb920-291a-4b54-963f-fd9bd610662f
# ╟─feb21ce6-0ddc-45fb-90f4-1e46261a9110
# ╟─3fda3dd7-1603-463c-9be4-03a217ada56f
# ╟─d3717c34-9fbf-470d-99f0-dabc3cc023fc
# ╟─68481e52-31d8-4208-b687-f8002ad27232
# ╟─12f268d2-baa6-4c1b-bed5-e9df53b469fc
