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

# ╔═╡ 9f57cb53-5fd0-41e5-8d8c-a68a92038767
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
    using VoronoiFVM
    using ExtendableGrids
    using PlutoUI
    using GridVisualize
    if isdefined(Main, :PlutoRunner)
        using CairoMakie
        CairoMakie.activate!(; type = "png")
        GridVisualize.default_plotter!(CairoMakie)
    end
end;

# ╔═╡ 5e13b3db-570c-4159-939a-7e2268f0a102
md"""
# Some problems with Voronoi FVM
[Source](https://github.com/WIAS-PDELib/VoronoiFVM.jl/blob/master/pluto-examples/problemcase.jl)

Draft. J. Fuhrmann, Oct. 29. 2021. Updated Dec 19, 2021.

We discuss one of the critical cases for application the Voronoi finite volume method.
We provide some practical fix and opine that the finite element method probably has the same problems.
"""

# ╔═╡ 556480e0-94f1-4e47-be9a-3e1e0e99555c
TableOfContents(; aside = false)

# ╔═╡ fae47c55-eef8-4428-bb5f-45824978753d
md"""
## Transient problem

This problem was suggested by R. Eymard.
"""

# ╔═╡ 8ba2300c-17ff-44e1-b33a-c5bdf1ce12fe
md"""
Regard the following problem coupling Darcy's equation with Fick's law and transport:
"""

# ╔═╡ 51c9517c-8797-4406-b053-301694fb0484
md"""
```math
  \begin{aligned}
    \vec v &= k \nabla p \\
    \nabla \cdot \vec v &= 0\\
    \partial_t (\phi c) - \nabla \cdot (D\nabla c + c \vec v) &= 0
  \end{aligned}
```
"""

# ╔═╡ 99341e32-9c78-4e31-bec0-d1ffbc85ec32
md"""
The domain is described by the following discretization grid:
"""

# ╔═╡ cd013964-f329-4d2c-ae4b-305093f0ac56
md"""
### Results

In the calculations, we ramp up the inlet concentration and  measure  the amount  of solute flowing  through the outlet - the breaktrough curve.
"""

# ╔═╡ afccbadb-2ca8-4c3e-8c6d-c78df59d8d7e
nref = 1

# ╔═╡ dd9f8d38-8812-40ba-88c8-f873ec7d6121
tend = 100

# ╔═╡ 5f6ac608-b1a0-450e-910e-d7d8ea2ffae0
ε_fix = 1.0e-4

# ╔═╡ 5b60c7d4-7bdb-4989-b055-6695b9fdeedc
md"""
Here, we plot the solutions for the `grid_n` case and the `grid_f` case.
"""

# ╔═╡ 98ae56dd-d42d-4a93-bb0b-5956b6e981a3
md"""
Time: $(@bind t PlutoUI.Slider(1:tend/100:tend,show_value=true,default=tend*0.1))
"""

# ╔═╡ 99c3b54b-d458-482e-8aa0-d2c2b51fdf25
md"""
## Reaction-Diffusion problem

Here we solve the following problem:

```math
    -\nabla D \nabla u + R u = 0
```
where D is large in the high permeability region and small otherwise. R is a constant.

"""

# ╔═╡ eef85cd7-eba4-4c10-9e1d-38411179314d
md"""
### Results
"""

# ╔═╡ fcd066f1-bcd8-4479-a4e4-7b8c235336c4
md"""
## Discussion

### Transient case
As there will be nearly no flow
in  y-direction, we should  get the  very same  results in  all four
cases for small permeability values in the low permeability region.

In the `grid_n` case,  the heterogeneous control volumina  ovrestimate the storage
capacity which shows itself  in a underestimation  of the
transferred solute.

With  the high  permeability contrast,  the results  for heterogeneous
domain should be essentially equal to those for 1D domain.
 However,   with  a   coarse  resolution   in
y-direction, we see large  differences in the transient behaviour of
the breaktrough curve compared to the 1D case.
The introduction of a thin  protection layer leads  to  reasonable   results.


Alternatively, the porosity of the low permeability region can be modified.
Arguably, this is the case in practice, see e.g.
[Ackerer et al, Transport in Porous Media35:345–373, 1999](https://link.springer.com/content/pdf/10.1023/A:1006564309167.pdf)
(eq. 7).

### Reaction diffusion case
In this case, we look at a homogeneous reaction in a domain divided
into a high and low diffusion region. With high contrast in the diffusion
coefficients, the reasonable assumption is that the reaction takes place only
in the high diffusion region, and the un-consumed share of species leaves at the outlet.

In this case we observe a similar related problem which can be fixed by adding a thin layer
of control volumes at the boundary. No problem occurs if the reaction rate at in the low diffusion
region is zero.


### Conclusion

Here, we indeed observe problem with the Voronoi approach: care must be taken to handle the case
of hetero interfaces in connection with transient processes and/or homogeneous reactions.
In these cases it should be analyzed if the problem occurs, and why, and it appears, that the discussion
should not be had without reference to the correct physical models. A remedy based on meshing
is available at least for straight interfaces.

### Opinion

With standard ways of using finite elements, the issue described here will occur in a similar way, so
the discussion is indeed with the alternative "cell centered" finite volume approach which places interfaces
at the boundaries of the control volumes rather than along the edges of a underlying triangulation.

#### Drawbacks of two point flux Voronoi methods based on simplicial meshes (as tested here):
- Anisotropic diffusion is only correct with aligned meshes
- Reliance on boundary conforming Delaunay property of the underlying mesh, thus narrowing the available meshing strategies
- The issue described  in the present notebook. However, in both cases discussed here, IMHO it might  "go  away"  depending on the correct physics.
  There should be more discussions with relevant application problems at hand.

#### Advantages (compared to the cell centered approach placing collocation points away from interfaces)
- Availability of P1 interpolant on simplices for visualization, interpolation, coupling etc.
- Mesh generators tend to place interfaces at triangle edges.
- Dirichlet BC can be applied exactly
- There is a straightforward way to link interface processes with bulk processes, e.g. an adsorption reaction is easily described by a reaction term at the boundary which involves interface and bulk value available at the same mesh node.


"""

# ╔═╡ c9d92201-813c-499b-b863-b138c30eb634
md"""
## Appendix
"""

# ╔═╡ a372ac90-c871-4dc0-a44b-a5bddef71823
md"""
### Domain data
"""

# ╔═╡ 124b2a0a-ef19-453e-9e3a-5b5ce7db5fac
md"""
Sizes:
"""

# ╔═╡ 1ad18670-e7cb-4f7a-be0f-3db98cdeb6a4
begin
    L = 10   # length of the high perm layer
    W = 0.5  # width of high perm layer
    Wlow = 2 # width of adjacent low perm layers
end;

# ╔═╡ cc325b2c-6174-4b8d-8e39-202ac68b5705
md"""
In the center of the domain, we assume a layer with high permeability.

As  boundary  conditions for the pressure ``p`` we choose  fixed pressure values at  the left
and right boundaries of the  domain, triggering a constant pressure gradient throughout the domain.

At the inlet of the high  permeability layer, we set ``c=1``, and at the
outlet, we set ``c=0``.

The high permeability layer has length `L`=$( L) and width `W`= $(W).

We solve the time dependent problem on three types of  rectangular grids with the same
resolution in   $x$ direction and different variants to to handle the  high permeability
layer.


- `grid_n` - a "naive" grid which just resolves the permeability layer and the surrounding material with equally spaced (in y direction) grids
- `grid_1` - a 1D grid  of the high permeability layer. With high permeability contrast, the solution of the 2D case at y=0 should coincide with the 1D solution
- `grid_f` - a "fixed" 2D grid which resolves the permeability layer and the surrounding material with equally spaced (in y direction) grids and "protection layers" of width `ε_fix`=$(ε_fix)  correcting the size of high permeability control volumes


"""

# ╔═╡ 47bc8e6a-e296-42c9-bfc5-967edfb0feb7
md"""
Boundary conditions:
"""

# ╔═╡ d1d5bad2-d282-4e7d-adb9-baf21f58155e
begin
    const Γ_top = 3
    const Γ_bot = 1
    const Γ_left = 4
    const Γ_right = 2
    const Γ_in = 5
    const Γ_out = 2
end;

# ╔═╡ 9d736062-6821-46d9-9e49-34b43b78e814
begin
    Ω_low = 1
    Ω_high = 2
end;

# ╔═╡ 83b9931f-9020-4400-8aeb-31ad391184db
function grid_2d(; nref = 0, ε_fix = 0.0)
    nx = 10 * 2^nref
    ny = 1 * 2^nref
    nylow = 3 * 2^nref
    xc = linspace(0, L, nx + 1)
    y0 = linspace(-W / 2, W / 2, ny + 1)
    if ε_fix > 0.0
        yfix = [W / 2, W / 2 + ε_fix]
        ytop = glue(yfix, linspace(yfix[end], Wlow, nylow + 1))
    else
        ytop = linspace(W / 2, Wlow, nylow + 1)
    end
    yc = glue(-reverse(ytop), glue(y0, ytop))
    grid = simplexgrid(xc, yc)
    cellmask!(grid, [0, -W / 2], [L, W / 2], Ω_high)
    bfacemask!(grid, [0, -W / 2], [0, W / 2], Γ_in)
    return bfacemask!(grid, [L, -W / 2], [L, W / 2], Γ_out)
end

# ╔═╡ 46a0f078-4165-4e37-9e69-e69af8584f6e
gridplot(grid_2d(); resolution = (400, 300))

# ╔═╡ 3f693666-4026-4c01-a7aa-8c7dcbc32372
gridplot(grid_2d(; ε_fix = 1.0e-1); resolution = (400, 300))

# ╔═╡ c402f03c-746a-45b8-aaac-902a2f196094
function grid_1d(; nref = 0)
    nx = 10 * 2^nref
    xc = linspace(0, L, nx + 1)
    grid = simplexgrid(xc)
    cellmask!(grid, [0], [L], Ω_high)
    bfacemask!(grid, [0], [0], Γ_in)
    bfacemask!(grid, [L], [L], Γ_out)
    return grid
end

# ╔═╡ d772ac1b-3cda-4a2b-b0a9-b22b63b30653
md"""
### Transient solver
"""

# ╔═╡ a63a655c-e48b-4969-9409-31cd3db3bdaa
md"""
Pressure index in solution
"""

# ╔═╡ d7009231-4b43-44bf-96ba-9a203c0b5f5a
const ip = 1;

# ╔═╡ 26965e38-91cd-4022-bdff-4c503f724bfe
md"""
Concentration index in solution
"""

# ╔═╡ c904c921-fa10-43eb-bd46-b2869fa7f431
const ic = 2;

# ╔═╡ b143c846-2294-47f7-a2d1-8a6eabe942a3
md"""
Generate breaktrough courve from transient solution
"""

# ╔═╡ 92e4e4ab-3485-4cb9-9b41-e702a211a477
function breakthrough(sys, tf, sol)
    of = similar(sol.t)
    t = sol.t
    of[1] = 0
    for i in 2:length(sol.t)
        of[i] = -integrate(sys, tf, sol.u[i], sol.u[i - 1], t[i] - t[i - 1])[ic]
    end
    return of
end

# ╔═╡ 3df8bace-b4f1-4052-84f7-dff21d3a35f0
md"""
Transient solver:
"""

# ╔═╡ e866db69-9388-4691-99f7-879cf0658418
function trsolve(
        grid;
        κ = [1.0e-3, 5],
        D = [1.0e-12, 1.0e-12],
        Δp = 1.0,
        ϕ = [1, 1],
        tend = 100,
    )
    function flux(y, u, edge, data)
        y[ip] = κ[edge.region] * (u[ip, 1] - u[ip, 2])
        bp, bm = fbernoulli_pm(y[ip] / D[edge.region])
        y[ic] = D[edge.region] * (bm * u[ic, 1] - bp * u[ic, 2])
        return nothing
    end

    function stor(y, u, node, data)
        y[ip] = 0
        y[ic] = ϕ[node.region] * u[ic]
        return nothing
    end

    dim = dim_space(grid)
    function bc(y, u, bnode, data)
        c0 = ramp(bnode.time; dt = (0, 0.001), du = (0, 1))
        boundary_dirichlet!(y, u, bnode, ic, Γ_in, c0)
        boundary_dirichlet!(y, u, bnode, ic, Γ_out, 0)

        boundary_dirichlet!(y, u, bnode, ip, Γ_in, Δp)
        boundary_dirichlet!(y, u, bnode, ip, Γ_out, 0)
        if dim > 1
            boundary_dirichlet!(y, u, bnode, ip, Γ_left, Δp)
            boundary_dirichlet!(y, u, bnode, ip, Γ_right, 0)
        end
        return nothing
    end

    sys = VoronoiFVM.System(
        grid;
        check_allocs = true,
        flux = flux,
        storage = stor,
        bcondition = bc,
        species = [ip, ic],
    )

    inival = solve(sys; inival = 0, time = 0.0)
    factory = TestFunctionFactory(sys)
    tfc = testfunction(factory, [Γ_in, Γ_left, Γ_top, Γ_bot], [Γ_out])

    sol = VoronoiFVM.solve(
        sys;
        inival = inival,
        times = [0, tend],
        Δt = 1.0e-4,
        Δt_min = 1.0e-6,
    )

    bt = breakthrough(sys, tfc, sol)
    if dim == 1
        bt = bt * W
    end

    return grid, sol, bt
end

# ╔═╡ cd88123a-b042-43e2-99b9-ec925a8794ed
grid_n, sol_n, bt_n = trsolve(grid_2d(; nref = nref); tend = tend);

# ╔═╡ 1cf0db37-42cc-4dd9-9da3-ebb94ff63b1b
sum(bt_n)

# ╔═╡ c52ed973-2250-423a-b427-e91972f7ce74
@test sum(bt_n) ≈ 18.143158169851787

# ╔═╡ 732e79fa-5b81-4401-974f-37ea3427e770
scalarplot(grid_n, sol_n(t)[ic, :]; resolution = (500, 200), show = true)

# ╔═╡ b0ad0adf-6f6c-4fb3-b58e-e05cc8c0c796
grid_1, sol_1, bt_1 = trsolve(grid_1d(; nref = nref); tend = tend);

# ╔═╡ 02330841-fdf9-4ebe-9da6-cf96529b223c
@test sum(bt_1) ≈ 20.66209910195916

# ╔═╡ e36d2aef-1b5a-45a7-9289-8d1e544bcedd
scalarplot(
    grid_1,
    sol_1(t)[ic, :];
    levels = 0:0.2:1,
    resolution = (500, 150),
    xlabel = "x",
    ylabel = "c",
    title = "1D calculation, t=$t"
)

# ╔═╡ 76b77ec0-27b0-4a02-9ae4-43d756eb09dd
grid_f, sol_f, bt_f = trsolve(grid_2d(; nref = nref, ε_fix = ε_fix); tend = tend);

# ╔═╡ d23d6634-266c-43e3-9493-b61fb390bbe7
@test sum(bt_f) ≈ 20.661131375044135

# ╔═╡ f42d4eb6-3e07-40c9-a8b3-dc772e674222
scalarplot(grid_f, sol_f(t)[ic, :]; resolution = (500, 200), show = true)

# ╔═╡ 904b36f0-10b4-4db6-9252-21668305de9c
grid_ϕ, sol_ϕ, bt_ϕ = trsolve(grid_2d(; nref = nref); ϕ = [1.0e-3, 1], tend = tend);

# ╔═╡ b260df8a-3721-4203-bc0c-a23bcab9a311
@test sum(bt_ϕ) ≈ 20.412256299447236

# ╔═╡ ce49bb25-b2d0-4d17-a8fe-d7b62e9b20be
begin
    p1 = GridVisualizer(;
        resolution = (500, 200),
        xlabel = "t",
        ylabel = "outflow",
        legend = :rb,
        title = "Breakthrough Curves"
    )
    scalarplot!(p1, sol_n.t, bt_n; label = "naive grid", color = :red)
    scalarplot!(
        p1,
        sol_1.t,
        bt_1;
        label = "1D grid",
        markershape = :x,
        markersize = 10,
        clear = false,
        color = :green
    )
    scalarplot!(
        p1,
        sol_f.t,
        bt_f;
        label = "grid with fix",
        markershape = :circle,
        color = :green,
        clear = false
    )
    scalarplot!(
        p1,
        sol_ϕ.t,
        bt_ϕ;
        label = "modified ϕ",
        markershape = :cross,
        color = :blue,
        clear = false
    )
    reveal(p1)
end

# ╔═╡ 78d92b4a-bdb1-4117-ab9c-b422eac403b1
md"""
### Reaction-Diffusion solver
"""

# ╔═╡ bb3a50ed-32e7-4305-87d8-4093c054a4d2
function rdsolve(grid; D = [1.0e-12, 1.0], R = [1, 0.1])
    function flux(y, u, edge, data)
        y[1] = D[edge.region] * (u[1, 1] - u[1, 2])
        return nothing
    end

    function rea(y, u, node, data)
        y[1] = R[node.region] * u[1]
        return nothing
    end
    function bc(y, u, bnode, data)
        boundary_dirichlet!(y, u, bnode, 1, Γ_in, 1)
        boundary_dirichlet!(y, u, bnode, 1, Γ_out, 0)
        return nothing
    end
    sys = VoronoiFVM.System(
        grid;
        flux = flux,
        reaction = rea,
        species = 1,
        bcondition = bc,
        check_allocs = true
    )
    dim = dim_space(grid)

    sol = VoronoiFVM.solve(sys)
    factory = TestFunctionFactory(sys)
    tf = testfunction(factory, [Γ_in, Γ_left, Γ_top, Γ_bot], [Γ_out])
    of = integrate(sys, tf, sol)
    fac = 1.0
    if dim == 1
        fac = W
    end
    return grid, sol[1, :], of[1] * fac
end

# ╔═╡ 2f560406-d169-4027-9cfe-7689494edf45
rdgrid_1, rdsol_1, of_1 = rdsolve(grid_1d(; nref = nref));

# ╔═╡ 40850999-12da-46cd-b86c-45808592fb9e
@test of_1 ≈ -0.013495959676585267

# ╔═╡ 34228382-4b1f-4897-afdd-19db7d5a7c59
scalarplot(rdgrid_1, rdsol_1; resolution = (300, 200))

# ╔═╡ a6714eac-9e7e-4bdb-beb7-aca354664ad6
rdgrid_n, rdsol_n, of_n = rdsolve(grid_2d(; nref = nref));

# ╔═╡ d1bfac0f-1f20-4c0e-9a9f-c7d36bc338ef
@test of_n ≈ -0.00023622450350365264

# ╔═╡ 5899df30-5198-4946-a148-108746cdde79
scalarplot(rdgrid_n, rdsol_n; resolution = (500, 200))

# ╔═╡ 20d7624b-f43c-4ac2-bad3-383a9e4e1b42
rdgrid_f, rdsol_f, of_f = rdsolve(grid_2d(; nref = nref, ε_fix = ε_fix));

# ╔═╡ 5d407d63-8a46-4480-94b4-80510eac5166
@test of_f ≈ -0.013466874615165499

# ╔═╡ 6a6d0e94-8f0d-4119-945c-dd48ec0798fd
scalarplot(rdgrid_f, rdsol_f; resolution = (500, 200))

# ╔═╡ c0fc1f71-52ba-41a9-92d1-74e82ac7826c
rdgrid_r, rdsol_r, of_r = rdsolve(grid_2d(; nref = nref); R = [0, 0.1]);

# ╔═╡ 43622531-b7d0-44d6-b840-782021eb2ef0
@test of_r ≈ -0.013495959676764535

# ╔═╡ c08e86f6-b5c2-4762-af23-382b1b153f45
md"""
We measure the outflow at the outlet. As a result, we obtain:
   - 1D case: $(of_1)
   - 2D case, naive grid: $(of_n)
   - 2D case, grid with "protective layer": $(of_f)
   - 2D case, naive grid, "modified" R: $(of_r)
"""

# ╔═╡ 0cc1c511-f351-421f-991a-a27f26a8db4f
html"<hr>"

# ╔═╡ 5aa4f27f-9d57-4c10-8f77-47c8b44963e3
html"<hr>"

# ╔═╡ Cell order:
# ╟─5e13b3db-570c-4159-939a-7e2268f0a102
# ╟─556480e0-94f1-4e47-be9a-3e1e0e99555c
# ╠═9f57cb53-5fd0-41e5-8d8c-a68a92038767
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╟─fae47c55-eef8-4428-bb5f-45824978753d
# ╟─8ba2300c-17ff-44e1-b33a-c5bdf1ce12fe
# ╟─51c9517c-8797-4406-b053-301694fb0484
# ╟─99341e32-9c78-4e31-bec0-d1ffbc85ec32
# ╟─46a0f078-4165-4e37-9e69-e69af8584f6e
# ╟─cc325b2c-6174-4b8d-8e39-202ac68b5705
# ╟─3f693666-4026-4c01-a7aa-8c7dcbc32372
# ╟─cd013964-f329-4d2c-ae4b-305093f0ac56
# ╠═afccbadb-2ca8-4c3e-8c6d-c78df59d8d7e
# ╠═dd9f8d38-8812-40ba-88c8-f873ec7d6121
# ╠═5f6ac608-b1a0-450e-910e-d7d8ea2ffae0
# ╠═cd88123a-b042-43e2-99b9-ec925a8794ed
# ╠═1cf0db37-42cc-4dd9-9da3-ebb94ff63b1b
# ╠═c52ed973-2250-423a-b427-e91972f7ce74
# ╠═b0ad0adf-6f6c-4fb3-b58e-e05cc8c0c796
# ╠═02330841-fdf9-4ebe-9da6-cf96529b223c
# ╠═76b77ec0-27b0-4a02-9ae4-43d756eb09dd
# ╠═d23d6634-266c-43e3-9493-b61fb390bbe7
# ╠═904b36f0-10b4-4db6-9252-21668305de9c
# ╠═b260df8a-3721-4203-bc0c-a23bcab9a311
# ╠═ce49bb25-b2d0-4d17-a8fe-d7b62e9b20be
# ╟─5b60c7d4-7bdb-4989-b055-6695b9fdeedc
# ╟─e36d2aef-1b5a-45a7-9289-8d1e544bcedd
# ╟─98ae56dd-d42d-4a93-bb0b-5956b6e981a3
# ╠═732e79fa-5b81-4401-974f-37ea3427e770
# ╠═f42d4eb6-3e07-40c9-a8b3-dc772e674222
# ╟─99c3b54b-d458-482e-8aa0-d2c2b51fdf25
# ╟─eef85cd7-eba4-4c10-9e1d-38411179314d
# ╠═2f560406-d169-4027-9cfe-7689494edf45
# ╠═40850999-12da-46cd-b86c-45808592fb9e
# ╠═a6714eac-9e7e-4bdb-beb7-aca354664ad6
# ╠═d1bfac0f-1f20-4c0e-9a9f-c7d36bc338ef
# ╠═20d7624b-f43c-4ac2-bad3-383a9e4e1b42
# ╠═5d407d63-8a46-4480-94b4-80510eac5166
# ╠═c0fc1f71-52ba-41a9-92d1-74e82ac7826c
# ╠═43622531-b7d0-44d6-b840-782021eb2ef0
# ╟─c08e86f6-b5c2-4762-af23-382b1b153f45
# ╠═34228382-4b1f-4897-afdd-19db7d5a7c59
# ╠═5899df30-5198-4946-a148-108746cdde79
# ╠═6a6d0e94-8f0d-4119-945c-dd48ec0798fd
# ╟─fcd066f1-bcd8-4479-a4e4-7b8c235336c4
# ╟─c9d92201-813c-499b-b863-b138c30eb634
# ╟─a372ac90-c871-4dc0-a44b-a5bddef71823
# ╟─124b2a0a-ef19-453e-9e3a-5b5ce7db5fac
# ╠═1ad18670-e7cb-4f7a-be0f-3db98cdeb6a4
# ╟─47bc8e6a-e296-42c9-bfc5-967edfb0feb7
# ╠═d1d5bad2-d282-4e7d-adb9-baf21f58155e
# ╠═9d736062-6821-46d9-9e49-34b43b78e814
# ╠═83b9931f-9020-4400-8aeb-31ad391184db
# ╠═c402f03c-746a-45b8-aaac-902a2f196094
# ╟─d772ac1b-3cda-4a2b-b0a9-b22b63b30653
# ╟─a63a655c-e48b-4969-9409-31cd3db3bdaa
# ╠═d7009231-4b43-44bf-96ba-9a203c0b5f5a
# ╟─26965e38-91cd-4022-bdff-4c503f724bfe
# ╠═c904c921-fa10-43eb-bd46-b2869fa7f431
# ╟─b143c846-2294-47f7-a2d1-8a6eabe942a3
# ╠═92e4e4ab-3485-4cb9-9b41-e702a211a477
# ╟─3df8bace-b4f1-4052-84f7-dff21d3a35f0
# ╠═e866db69-9388-4691-99f7-879cf0658418
# ╟─78d92b4a-bdb1-4117-ab9c-b422eac403b1
# ╠═bb3a50ed-32e7-4305-87d8-4093c054a4d2
# ╟─0cc1c511-f351-421f-991a-a27f26a8db4f
# ╟─5aa4f27f-9d57-4c10-8f77-47c8b44963e3
