### A Pluto.jl notebook ###
# v0.20.20

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

# ╔═╡ b285aca3-dee5-4b77-9276-537563e8643b
begin
    import Pkg as _Pkg
    haskey(ENV, "PLUTO_PROJECT") && _Pkg.activate(ENV["PLUTO_PROJECT"])
    using Revise
    using VoronoiFVM
    using ExtendableGrids
    using ExtendableSparse
    using Test
    using PlutoUI
    using GridVisualize
    using LinearSolve
    using ILUZero
    using LinearAlgebra
    if isdefined(Main, :PlutoRunner)
        using CairoMakie
        CairoMakie.activate!(; type = "png", visible = false)
        GridVisualize.default_plotter!(CairoMakie)
    end
end;

# ╔═╡ 4ed0c302-26e4-468a-a40d-0e6406f802d0
md"""
# API updates 

[Source](https://github.com/WIAS-PDELib/VoronoiFVM.jl/blob/master/pluto-examples/api-updates.jl)
"""

# ╔═╡ 3e6b4ffa-7b33-4b94-9fd6-75b030d5a115
md"""
Here we describe some updates for the API of `VoronoiFVM.jl`. These have been implemented mostly on top of the existing API, whose functionality is not affected.
"""

# ╔═╡ 7a104243-d3b9-421a-b494-5607c494b106
TableOfContents(; aside = false, depth = 5)

# ╔═╡ a2f1e6ba-80b2-4902-9c5d-2806a7fb16f6
md"""
## v0.19
"""

# ╔═╡ 56976316-ca4f-4d5c-b303-026edd8751c2
md"""
This is a breaking release. Implementations using default solver settings should continue to work (albeit possibly with deprecation and 
allocation warnings). Really breaking is control of iterative linear solvers and allocation checks.

"""

# ╔═╡ 740dd980-ca12-4f51-88b2-408930704952
md"""
### Solve now a method of CommonSolve.solve
"""

# ╔═╡ db216ebf-1869-4b0f-8dbd-2a8a244c3e4c
md"""
As a consequence, all VoronoiFVM.solve methods with signatures others than `solve(system; kwargs...)`  are now deprecated
"""

# ╔═╡ 9eae0139-379d-47ce-a8d2-932c68504317
n = 100

# ╔═╡ 8a31ed03-9cd8-4532-9fb9-8c060f777693
begin
    h = 1.0 / convert(Float64, n)
    const eps = 1.0e-2
    function reaction(f, u, node, data)
        f[1] = u[1]^2
        return nothing
    end

    function flux(f, u, edge, data)
        f[1] = eps * (u[1, 1]^2 - u[1, 2]^2)
        return nothing
    end

    function source(f, node, data)
        x1 = node[1] - 0.5
        x2 = node[2] - 0.5
        f[1] = exp(-20.0 * (x1^2 + x2^2))
        return nothing
    end

    function storage(f, u, node, data)
        f[1] = u[1]
        return nothing
    end

    function bcondition(f, u, node, data)
        boundary_dirichlet!(
            f,
            u,
            node;
            species = 1,
            region = 2,
            value = ramp(node.time; dt = (0, 0.1), du = (0, 1))
        )
        boundary_dirichlet!(
            f,
            u,
            node;
            species = 1,
            region = 4,
            value = ramp(node.time; dt = (0, 0.1), du = (0, 1))
        )
        return nothing
    end

    sys0 = VoronoiFVM.System(
        0.0:h:1.0,
        0.0:h:1.0;
        reaction,
        flux,
        source,
        storage,
        bcondition,
        species = [1],
    )
end

# ╔═╡ 4279ae2e-a948-4358-9037-0c6895ecb809
md"""
Deprecated call:
"""

# ╔═╡ 58cf3e44-ee53-42f7-9132-eacfd900dc3a
md"""
begin
    inival = unknowns(sys0; inival = 0.1)
    sol00 = unknowns(sys0)
    solve!(sol00, inival, sys0)
end
"""

# ╔═╡ b499a195-f06d-400f-9407-b07c3212c095
md"""
Replace this by:
"""

# ╔═╡ f660fc6e-0ab5-4fae-918f-39bf0c2153e5
sol0 = solve(sys0; inival = 0.1)

# ╔═╡ 13f284e8-4fe9-42ce-873a-c65febc7d7df
md"""
#### Docstring of solve
"""

# ╔═╡ c53cb54f-6099-4d09-9797-3da1f5428586
(@doc solve).content[end]

# ╔═╡ 2cebb072-bef1-4bce-9cbb-6b43b877b28b
md"""
#### Docstring of SolverControl
"""

# ╔═╡ 7b0f2021-a2fd-4bb2-a23a-432f61a38a07
@doc SolverControl

# ╔═╡ 182a9a9f-3f2c-4df3-abad-08488bb9fb33
md"""
### Rely on LinearSolve.jl for linear system solution
"""

# ╔═╡ 381c464e-a0cc-46ea-b0ba-089250c31555
md"""
This provides easy access to a large variety of linear solvers:
"""

# ╔═╡ ba4b8b92-617a-40b4-b3c9-1367067027fa
md"""
#### LU factorization from UMFPACK
"""

# ╔═╡ 00634dc4-e9c7-4165-a9ac-f3ffc8007e76
umf_sol = solve(sys0; inival = 0.1, method_linear = UMFPACKFactorization(), verbose = true)

# ╔═╡ d6c88c3e-6b8a-45e9-9344-71f0de3fff51
@test isapprox(umf_sol, sol0, atol = 1.0e-7)

# ╔═╡ cb2d42c4-8d5d-465c-9376-f568588ab453
md"""
#### LU factorization from Sparspak.jl
"""

# ╔═╡ 5908ae1d-b3b5-4681-a0b8-080f052af40f
sppk_sol = solve(sys0; inival = 0.1, method_linear = SparspakFactorization(), verbose = true)

# ╔═╡ 22c3cecc-12d2-4f7d-8a53-894a5ea513f0
@test isapprox(sppk_sol, sol0, atol = 1.0e-7)

# ╔═╡ 3107c644-5f43-4322-aad2-04a2ab0d576e
md"""
#### Iterative solvers
"""

# ╔═╡ 7e94ff44-4807-41b0-875d-526390903942
md"""
##### BICGstab from Krylov.jl with diagonal (Jacobi) preconditioner
The Jacobi preconditioner is defined in ExtendableSparse.jl.
"""

# ╔═╡ b70524fc-b8b1-4dee-b77d-e3f8d6d2837b
krydiag_sol = solve(
    sys0;
    inival = 0.1,
    method_linear = KrylovJL_BICGSTAB(),
    precon_linear = JacobiPreconditioner,
    verbose = true,
)

# ╔═╡ d21d3236-b3d7-4cf8-ab9d-b0be44c9970b
@test isapprox(krydiag_sol, sol0, atol = 1.0e-5)

# ╔═╡ f935061c-71db-4aaf-864b-bb566e68643e
md"""
##### BICGstab from Krylov.jl with delayed factorization preconditioner
"""

# ╔═╡ 4c7f8bbf-40c4-45b9-a62e-99ffaae30af1
krydel_sol = solve(
    sys0;
    inival = 0.1,
    method_linear = KrylovJL_BICGSTAB(),
    precon_linear = SparspakFactorization(),
    verbose = "nlad",
)

# ╔═╡ 4fa1c608-19b2-4eaa-8c0a-5881b373807c
@test isapprox(krydel_sol, sol0, atol = 1.0e-5)

# ╔═╡ 023ae612-6bef-4b0d-8d04-5c0461efbc18
md"""
##### BICGstab from Krylov.jl with ilu0 preconditioner 
`ILUZeroPreconditioner` is exported from ExtendableSparse and
wraps the predonditioner defined in  ILUZero.jl .
"""

# ╔═╡ 6895cdf9-8291-47ce-bd1d-4c5beec594ea
kryilu0_sol = solve(
    sys0;
    inival = 0.5,
    method_linear = KrylovJL_BICGSTAB(),
    precon_linear = ILUZeroPreconditioner,
    verbose = true,
)

# ╔═╡ d341f60e-191d-4cf9-9df9-fbe25c84a7da
@test isapprox(kryilu0_sol, sol0, atol = 1.0e-5)

# ╔═╡ 3bca80f4-ec88-41e8-829a-278071442d41
md"""
### New verbosity handling
"""

# ╔═╡ c2c3088a-6c6e-492e-957b-c5d20fec4244
md"""
- `verbose` can now be a Bool or a String of flag characters, allowing for control of different output categories. I would love to do this via  logging, but there is still a [long way to go](https://github.com/JuliaLang/julia/issues/33418) IMHO 
- Allocation check is active by default with warnings which can be muted by passing a `verbose` string without 'a'. This is now the only control in this respect. All `check_allocs` methods/kwargs, control via environment variables have been removed.
- Deprecation warnings can be switched off by passing a `verbose` string without 'd'.
- Improve iteration logging etc., allow for logging of linear iterations ('l' flag character)


"""

# ╔═╡ ba6b0e35-75ce-4540-861c-7654fd4dee63
md"""
The following example gives some information in this respect:
"""

# ╔═╡ d1440945-54ac-4c12-9f34-0db1c7dfac11
D = 0.1

# ╔═╡ 2ff77259-7a81-4c54-a291-cbc20ee56c5d
function xflux(f, u, edge, data)
    return f[1] = D * (u[1, 1]^2 - u[1, 2]^2)
end

# ╔═╡ b0a845d7-95d7-4212-9620-e1948698c596
xsys = VoronoiFVM.System(0:0.001:1; flux = xflux, species = [1])

# ╔═╡ 2a85a42c-1a79-425f-8887-71e2944cb0f3
solve(xsys; inival = 0.1, times = [0, 1]);

# ╔═╡ 1370a5dc-c434-423f-a0fa-b25f0f2878f9
md"""
If we find these warnings annoying, we can switch them off:
"""

# ╔═╡ 62f605a0-6d5c-4ce2-a131-9ab7b6188d23
solve(xsys; inival = 0.1, times = [0, 1], verbose = "");

# ╔═╡ 2b58318c-8e21-4f9a-97df-ff4af50c94b2
md"""
Or we get some more logging:
"""

# ╔═╡ 05381762-d8f8-46d2-8eb2-68275458787a
solve(xsys; inival = 0.1, times = [0, 1], verbose = "en");

# ╔═╡ f752f390-c98e-4832-b186-f484ebe5a4cb
md"""
But we can also look for the reasons of the allocations. Here, global values should be declared as constants.
"""

# ╔═╡ 0e8bbe4c-66d8-4545-8196-c7d8d9e30bfd
const D1 = 0.1

# ╔═╡ c40f1954-4fb7-48b2-ab4c-cdf6459b7383
function xflux1(f, u, edge, data)
    f[1] = D1 * (u[1, 1]^2 - u[1, 2]^2)
    return nothing
end

# ╔═╡ cec5152b-0597-4ea4-8387-b08e1e4ffcde
xsys1 = VoronoiFVM.System(0:0.001:1; flux = xflux1, species = [1])

# ╔═╡ 76c6ee48-8061-4ba0-b61e-0e6b68ad6435
solve(xsys1; inival = 0.1, times = [0, 1]);

# ╔═╡ d4ee4693-8ecd-4916-a722-79f54eb99d42
md"""
## v0.14
"""

# ╔═╡ 8a4f336c-2016-453e-9a9f-beac66533fc0
md"""
### `VoronoiFVM.System` constructor
"""

# ╔═╡ 78e7c000-3a83-446a-b577-3a1809c664d2
md"""
### Implicit creation of physics

The `VoronoiFVM.Physics` struct almost never was used outside of the constructor of `VoronoiFVM.System`. Now it is possible to specify the flux functions directly in the system constructor. By default, it is also possible to set a list of species which are attached to all interior and boundary regions of the grid.
"""

# ╔═╡ c59876bd-0cb1-4157-9ba4-bdbede151a44
grid1 = simplexgrid(0:0.1:1);

# ╔═╡ 90bbf212-c6c8-44f0-8132-4a98f094750e
function multispecies_flux(y, u, edge, data)
    for i in 1:(edge.nspec)
        y[i] = u[i, 1] - u[i, 2]
    end
    return nothing
end

# ╔═╡ adff41d1-9398-4a66-9a8e-e03809973fa6
function test_reaction(y, u, node, data)
    y[1] = u[1]
    y[2] = -u[1]
    return nothing
end

# ╔═╡ 5e6d83ab-65c7-4f33-b0a8-29cd5717b4d6
begin
    system1 = VoronoiFVM.System(
        grid1;
        flux = multispecies_flux,
        reaction = test_reaction,
        species = [1, 2]
    )
    boundary_dirichlet!(system1; species = 1, region = 1, value = 1)
    boundary_dirichlet!(system1; species = 2, region = 2, value = 0)
end;

# ╔═╡ f11c03a3-7899-42fd-a2da-a257715815dc
sol1 = solve(system1);

# ╔═╡ 2754d4c8-bbc1-4283-8156-c660c33cd62d
let
    vis = GridVisualizer(; resolution = (500, 300), legend = :rt)
    scalarplot!(vis, grid1, sol1[1, :]; color = :red, label = "species1")
    scalarplot!(vis, grid1, sol1[2, :]; color = :green, label = "species2", clear = false)
    reveal(vis)
end

# ╔═╡ fbd75cf1-64e4-4f07-b54f-f90626f3f6ba
@test isapprox(sum(sol1), 11.323894375033476, rtol = 1.0e-14)

# ╔═╡ b5f7e133-500d-4a27-8f78-11ff3582599c
md"""
### Boundary conditions as part of physics

This makes the API more consistent and opens an easy possibility to have
space and time dependent boundary conditions. One can specify them either in `breaction` or the synonymous `bcondition`.
"""

# ╔═╡ ec188c81-3374-4eed-9b7e-e22350886df2
function bcond2(y, u, bnode, data)
    boundary_neumann!(y, u, bnode; species = 1, region = 1, value = sin(bnode.time))
    boundary_dirichlet!(y, u, bnode; species = 2, region = 2, value = 0)
    return nothing
end;

# ╔═╡ c86e8a0f-299f-42ab-96f8-0cd62d50f196
system2 = VoronoiFVM.System(
    grid1;
    flux = multispecies_flux,
    reaction = test_reaction,
    species = [1, 2],
    bcondition = bcond2,
    check_allocs = false
);

# ╔═╡ b3d936fe-69ab-4013-b787-2f0b5410638a
sol2 = solve(system2; times = (0, 10), Δt_max = 0.01);

# ╔═╡ 17749697-d5d8-4629-a625-e96590a5f0ac
vis2 = GridVisualizer(; resolution = (500, 300), limits = (-2, 2), legend = :rt)

# ╔═╡ 0c916da5-2d6e-42df-ac4b-4a062f931ccd
md"""
time: $(@bind t2 PlutoUI.Slider(0:0.01:10; default=5,show_value=true))
"""

# ╔═╡ 783618f8-2470-4c7c-afc1-9800586625c1
let
    s = sol2(t2)
    scalarplot!(vis2, grid1, s[1, :]; color = :red, label = "species1")
    scalarplot!(
        vis2,
        grid1,
        s[2, :];
        color = :green,
        label = "species2",
        clear = false,
        title = "time=$(t2)"
    )
    reveal(vis2)
end

# ╔═╡ 4cbea340-9c02-4e69-8f5e-62bf45312bdd
@test isapprox(sum(sol2) / length(sol2), 2.4921650158811794, rtol = 1.0e-14)

# ╔═╡ 1c18b5a0-cca6-46a1-bb9f-b3d65b8043c5
md"""
### Implicit creation of grid
"""

# ╔═╡ 47280b56-e5ec-4345-b4a1-7c3c92536b2e
md"""
By passing data for grid creation (one to three abstract vectors) instead a grid, a tensor product grid is implicitly created.
This example also demonstrates position dependent boundary values.
"""

# ╔═╡ a71086fa-4ec6-4842-a4e1-6a6b60441fc2
function bcond3(y, u, bnode, data)
    boundary_dirichlet!(y, u, bnode; region = 4, value = bnode[2])
    boundary_dirichlet!(y, u, bnode; region = 2, value = -bnode[2])
    return nothing
end;

# ╔═╡ a514231a-e465-4f05-ba4c-b20aa968d96f
system3 = VoronoiFVM.System(
    -1:0.1:1,
    -1:0.1:1;
    flux = multispecies_flux,
    bcondition = bcond3,
    species = 1
);

# ╔═╡ d55f615c-d586-4ef7-adf9-5faf052b75ac
sol3 = solve(system3);

# ╔═╡ c17e5104-4d3a-4d54-81c1-d7253245a8bb
@test isapprox(sum(sol3), 0.0, atol = 1.0e-14)

# ╔═╡ 087ea16d-742e-4398-acf5-37248af1b5b4
md"""
### GridVisualize API extended to System
Instead of a grid, a system can be passed to `gridplot` and `scalarplot`.
"""

# ╔═╡ 34d465a5-7cc5-4348-b9ba-6d9381bb3a87
scalarplot(system3, sol3; resolution = (300, 300), levels = 10, colormap = :hot)

# ╔═╡ ac48f5bd-fd1e-4aa7-a2c9-90f0f427143c
md"""
### Parameters of `solve`
"""

# ╔═╡ 69974c02-57e6-4eb5-acf4-b2d480fbd67d
md"""
The `solve` API has been simplified and made more Julian. All entries of `VoronoiFVM.NewtonControl` can be now passed as keyword arguments to `solve`.

Another new keyword argument is `inival` which allows to pass an initial value which by default is initialized to zero. Therefore we now can write `solve(system)` as we already have seen above.
"""

# ╔═╡ 1e12afcf-cf46-4672-9434-44fa8af95ef7
reaction4(y, u, bnode, data) = y[1] = -bnode[1]^2 + u[1]^4;

# ╔═╡ 938ef63c-58c4-41a0-b3dd-4eb76987a4d7
bc4(f, u, node, data) = boundary_dirichlet!(f, u, node; value = 0);

# ╔═╡ fe424654-f070-46a9-850a-738b1d4aca8f
system4 = VoronoiFVM.System(
    -10:0.1:10;
    species = [1],
    reaction = reaction4,
    flux = multispecies_flux,
    bcondition = bc4
);

# ╔═╡ 37fc8816-5ccd-436e-8335-ebb1218d8a35
sol4 = solve(system4; log = true, damp_initial = 0.001, damp_growth = 3);

# ╔═╡ 6a256a29-f15f-4d82-8e84-7ceacb786715
scalarplot(
    system4,
    sol4;
    resolution = (500, 300),
    xlabel = "x",
    ylabel = "u",
    title = "solution"
)

# ╔═╡ 5c2a3836-dc81-4950-88e5-7f603514b1c0
@test isapprox(sum(sol4), 418.58515700568535, rtol = 1.0e-14)

# ╔═╡ fc0245fe-1bf2-45a3-aa7c-9cce8d7eef37
html"""<hr>"""

# ╔═╡ bef5d557-44f4-418b-935a-ebd0ffaf69d8
html"""<hr>"""


# ╔═╡ Cell order:
# ╟─4ed0c302-26e4-468a-a40d-0e6406f802d0
# ╟─3e6b4ffa-7b33-4b94-9fd6-75b030d5a115
# ╠═7a104243-d3b9-421a-b494-5607c494b106
# ╠═b285aca3-dee5-4b77-9276-537563e8643b
# ╟─a2f1e6ba-80b2-4902-9c5d-2806a7fb16f6
# ╟─56976316-ca4f-4d5c-b303-026edd8751c2
# ╟─740dd980-ca12-4f51-88b2-408930704952
# ╟─db216ebf-1869-4b0f-8dbd-2a8a244c3e4c
# ╠═9eae0139-379d-47ce-a8d2-932c68504317
# ╠═8a31ed03-9cd8-4532-9fb9-8c060f777693
# ╟─4279ae2e-a948-4358-9037-0c6895ecb809
# ╠═58cf3e44-ee53-42f7-9132-eacfd900dc3a
# ╟─b499a195-f06d-400f-9407-b07c3212c095
# ╠═f660fc6e-0ab5-4fae-918f-39bf0c2153e5
# ╟─13f284e8-4fe9-42ce-873a-c65febc7d7df
# ╟─c53cb54f-6099-4d09-9797-3da1f5428586
# ╟─2cebb072-bef1-4bce-9cbb-6b43b877b28b
# ╟─7b0f2021-a2fd-4bb2-a23a-432f61a38a07
# ╟─182a9a9f-3f2c-4df3-abad-08488bb9fb33
# ╟─381c464e-a0cc-46ea-b0ba-089250c31555
# ╟─ba4b8b92-617a-40b4-b3c9-1367067027fa
# ╠═00634dc4-e9c7-4165-a9ac-f3ffc8007e76
# ╠═d6c88c3e-6b8a-45e9-9344-71f0de3fff51
# ╟─cb2d42c4-8d5d-465c-9376-f568588ab453
# ╠═5908ae1d-b3b5-4681-a0b8-080f052af40f
# ╠═22c3cecc-12d2-4f7d-8a53-894a5ea513f0
# ╟─3107c644-5f43-4322-aad2-04a2ab0d576e
# ╟─7e94ff44-4807-41b0-875d-526390903942
# ╠═b70524fc-b8b1-4dee-b77d-e3f8d6d2837b
# ╠═d21d3236-b3d7-4cf8-ab9d-b0be44c9970b
# ╟─f935061c-71db-4aaf-864b-bb566e68643e
# ╠═4c7f8bbf-40c4-45b9-a62e-99ffaae30af1
# ╠═4fa1c608-19b2-4eaa-8c0a-5881b373807c
# ╟─023ae612-6bef-4b0d-8d04-5c0461efbc18
# ╠═6895cdf9-8291-47ce-bd1d-4c5beec594ea
# ╠═d341f60e-191d-4cf9-9df9-fbe25c84a7da
# ╟─3bca80f4-ec88-41e8-829a-278071442d41
# ╟─c2c3088a-6c6e-492e-957b-c5d20fec4244
# ╟─ba6b0e35-75ce-4540-861c-7654fd4dee63
# ╠═d1440945-54ac-4c12-9f34-0db1c7dfac11
# ╠═2ff77259-7a81-4c54-a291-cbc20ee56c5d
# ╠═b0a845d7-95d7-4212-9620-e1948698c596
# ╠═2a85a42c-1a79-425f-8887-71e2944cb0f3
# ╟─1370a5dc-c434-423f-a0fa-b25f0f2878f9
# ╠═62f605a0-6d5c-4ce2-a131-9ab7b6188d23
# ╟─2b58318c-8e21-4f9a-97df-ff4af50c94b2
# ╠═05381762-d8f8-46d2-8eb2-68275458787a
# ╟─f752f390-c98e-4832-b186-f484ebe5a4cb
# ╠═0e8bbe4c-66d8-4545-8196-c7d8d9e30bfd
# ╠═c40f1954-4fb7-48b2-ab4c-cdf6459b7383
# ╠═cec5152b-0597-4ea4-8387-b08e1e4ffcde
# ╠═76c6ee48-8061-4ba0-b61e-0e6b68ad6435
# ╟─d4ee4693-8ecd-4916-a722-79f54eb99d42
# ╟─8a4f336c-2016-453e-9a9f-beac66533fc0
# ╟─78e7c000-3a83-446a-b577-3a1809c664d2
# ╠═c59876bd-0cb1-4157-9ba4-bdbede151a44
# ╠═90bbf212-c6c8-44f0-8132-4a98f094750e
# ╠═adff41d1-9398-4a66-9a8e-e03809973fa6
# ╠═5e6d83ab-65c7-4f33-b0a8-29cd5717b4d6
# ╠═f11c03a3-7899-42fd-a2da-a257715815dc
# ╟─2754d4c8-bbc1-4283-8156-c660c33cd62d
# ╠═fbd75cf1-64e4-4f07-b54f-f90626f3f6ba
# ╟─b5f7e133-500d-4a27-8f78-11ff3582599c
# ╠═ec188c81-3374-4eed-9b7e-e22350886df2
# ╠═c86e8a0f-299f-42ab-96f8-0cd62d50f196
# ╠═b3d936fe-69ab-4013-b787-2f0b5410638a
# ╟─17749697-d5d8-4629-a625-e96590a5f0ac
# ╠═0c916da5-2d6e-42df-ac4b-4a062f931ccd
# ╟─783618f8-2470-4c7c-afc1-9800586625c1
# ╠═4cbea340-9c02-4e69-8f5e-62bf45312bdd
# ╟─1c18b5a0-cca6-46a1-bb9f-b3d65b8043c5
# ╟─47280b56-e5ec-4345-b4a1-7c3c92536b2e
# ╠═a71086fa-4ec6-4842-a4e1-6a6b60441fc2
# ╠═a514231a-e465-4f05-ba4c-b20aa968d96f
# ╠═d55f615c-d586-4ef7-adf9-5faf052b75ac
# ╠═c17e5104-4d3a-4d54-81c1-d7253245a8bb
# ╟─087ea16d-742e-4398-acf5-37248af1b5b4
# ╠═34d465a5-7cc5-4348-b9ba-6d9381bb3a87
# ╟─ac48f5bd-fd1e-4aa7-a2c9-90f0f427143c
# ╟─69974c02-57e6-4eb5-acf4-b2d480fbd67d
# ╠═1e12afcf-cf46-4672-9434-44fa8af95ef7
# ╠═938ef63c-58c4-41a0-b3dd-4eb76987a4d7
# ╠═fe424654-f070-46a9-850a-738b1d4aca8f
# ╠═37fc8816-5ccd-436e-8335-ebb1218d8a35
# ╟─6a256a29-f15f-4d82-8e84-7ceacb786715
# ╠═5c2a3836-dc81-4950-88e5-7f603514b1c0
# ╟─fc0245fe-1bf2-45a3-aa7c-9cce8d7eef37
# ╟─bef5d557-44f4-418b-935a-ebd0ffaf69d8
