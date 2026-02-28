### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ d5fcb04c-77f3-4a05-8497-a3eff38119a1
# Convert this cell to markdown in order to enable Pluto's inbuilt package manager
if isdefined(Main, :PlutoRunner)
    using Pkg
    docsdir = joinpath(@__DIR__, "..", "docs")
    if isdir(docsdir)
        Pkg.activate(docsdir)
    end
end

# ╔═╡ 18c423cc-18bf-41a0-a6e4-e30f91f39728
begin
    using Revise
    using VoronoiFVM
    using ExtendableGrids
    using GridVisualize
    using PlutoUI
    using HypertextLiteral
    using LinearAlgebra
    using LinearSolve
    using Test
    if isdefined(Main, :PlutoRunner)
        using CairoMakie
        CairoMakie.activate!(; type = "png")
        GridVisualize.default_plotter!(CairoMakie)
    end
end

# ╔═╡ 327af4a8-74cc-4834-ab19-d3a1d6873982
TableOfContents(; title = "", aside = false)

# ╔═╡ 3fa189e4-9e1c-470c-bf26-15b631945d2d
md"""
# Interface conditions in 1D
[Source](https://github.com/WIAS-PDELib/VoronoiFVM.jl/blob/master/pluto-examples/interfaces1d.jl)

This notebooks discusses handling of internal interfaces with VoronoiFVM.jl.

## Two subdomains
For a simple stationary diffusion equation with an interior interface, we discuss possible interface conditions between two subdomains.

Let ``\Omega=\Omega_1\cup\Omega_2`` where ``\Omega_1=(-1,0)`` and ``\Omega_2=(0,1)``.
Let ``\Gamma_1={-1}``,``\Gamma_2={1}``  and ``\Gamma_3={0}``.


Regard the following problem:


``\begin{aligned}
     -\Delta u_1 &= 0 & \text{in}\quad \Omega_1\\ 
     -\Delta u_2 &= 0 & \text{in}\quad \Omega_2\\ 
\end{aligned}
``

with exterior boundary conditions

``u_1|_{\Gamma_1} = g_1`` and ``u_2|_{\Gamma_2} = g_2`` 



For the interior boundary (interface) conditions we set 


``\nabla u_1|_{\Gamma_3}+f_1(u_1,u_2)=0``


``-\nabla u_2|_{\Gamma_3}+f_2(u_1,u_2)=0``

where ``f_1``, ``f_2`` are discussed later.
"""

# ╔═╡ d5a0ee0d-959d-476b-b3c5-79b741059992
md"""
### Set up
"""

# ╔═╡ f03ff283-c989-4b1a-b73e-2e616054e3db
md"""
Create a grid with two subdomins and an interface in the center.
"""

# ╔═╡ 670c78c1-d0be-4362-975b-2c944620681f
nref = 2

# ╔═╡ 79193d53-9dfa-47b7-aed2-c7eb43769b5f
begin
    hmax = 0.2 / 2.0^nref
    hmin = 0.05 / 2.0^nref
    X1 = geomspace(-1.0, 0.0, hmax, hmin)
    X2 = geomspace(0.0, 1.0, hmin, hmax)
    X = glue(X1, X2)
    grid = simplexgrid(X)

    bfacemask!(grid, [0.0], [0.0], 3)
    ## Material 1 left of 0
    cellmask!(grid, [-1.0], [0.0], 1)
    ## Material 2 right of 0
    cellmask!(grid, [0.0], [1.0], 2)
end;

# ╔═╡ 4cb07222-587b-4a74-a444-43f5433d5b03
gridplot(grid; legend = :rt, resolution = (600, 200))

# ╔═╡ 02ec3c0b-6e68-462b-84df-931370cbdcac
md"""
For later use (plotting) extract the two subgrids from the grid
"""

# ╔═╡ 592429a1-108c-4e9e-8961-497c2c31f319
subgrid1 = subgrid(grid, [1]);

# ╔═╡ 1fa5d9b3-0558-4558-a5f1-f34f70c8d9a0
subgrid2 = subgrid(grid, [2]);

# ╔═╡ 5b539ad2-46d6-43b0-8b3e-f8a7e1ae0a6d
md"""
Define the diffusion flux for the two species in their respective subdomains
"""

# ╔═╡ 6aabfbe1-de7d-49ba-8144-6d364b21b34f
function flux!(f, u, edge, data)
    if edge.region == 1
        f[1] = u[1, 1] - u[1, 2]
    end
    if edge.region == 2
        f[2] = u[2, 1] - u[2, 2]
    end
    return nothing
end

# ╔═╡ dbd27f10-9fd9-450d-9f78-89ea738d605b
md"""
Specify  the outer boundary values.
"""

# ╔═╡ e5b432b2-875e-49a5-8e78-68f7f47f06c3
const g_1 = 1.0

# ╔═╡ 8ae32292-df3e-4190-83a5-ec5ba529299e
const g_2 = 0.1

# ╔═╡ 677853a7-0c43-4d3c-bc74-799535f95aeb
md"""
Create the system. We pass the interface condition function as a parameter.
"""

# ╔═╡ 139058a9-44a0-43d5-a377-4fc72927fa28
function make_system(breaction)
    physics = VoronoiFVM.Physics(; flux = flux!, breaction = breaction)

    ## Create system
    sys = VoronoiFVM.System(grid, physics; unknown_storage = :sparse)

    ##  Enable species in their respective subregions
    enable_species!(sys, 1, [1])
    enable_species!(sys, 2, [2])

    ## Set boundary conditions
    for ispec in 1:2
        boundary_dirichlet!(sys, ispec, 1, g_1)
        boundary_dirichlet!(sys, ispec, 2, g_2)
    end
    return sys
end

# ╔═╡ de251b03-dd3a-4a44-9440-b7e654c32dac
md"""
Stationary solution with zero initial value
"""

# ╔═╡ 5b3a49d6-30ce-4d6e-94b4-a76645d7d8ce
function mysolve(sys)
    U = solve(sys)
    U1 = view(U[1, :], subgrid1)
    U2 = view(U[2, :], subgrid2)
    return U1, U2
end

# ╔═╡ 0735c061-68e1-429f-80f1-d8410989a91d
md"""
Plot the results
"""

# ╔═╡ 467dc381-3b3d-4de7-a7f9-bfc51300832b
function plot(U1, U2; title = "")
    vis = GridVisualizer(; resolution = (600, 300))
    scalarplot!(
        vis,
        subgrid1,
        U1;
        clear = false,
        show = false,
        color = :green,
        label = "u1"
    )
    return scalarplot!(
        vis,
        subgrid2,
        U2;
        clear = false,
        show = true,
        color = :blue,
        label = "u2",
        legend = :rt,
        title = title,
        flimits = (-0.5, 1.5)
    )
end

# ╔═╡ fa4fcc0d-1d3a-45a2-8857-50536bbe39cc
md"""
### No interface reaction

This means we set ``f_1(u_1,u_2)=0`` and ``f_2(u_1,u_2)=0``. 
"""

# ╔═╡ 8f210696-fcf4-47bc-a5a2-c561ad7efcbd
function noreaction(f, u, node, data)
    return nothing
end

# ╔═╡ 57e8515e-3be1-4478-af98-430501438ee7
system1 = make_system(noreaction);

# ╔═╡ 56136cd1-0c01-449d-9297-68924ac99ee7
plot(mysolve(system1)...)

# ╔═╡ fa1293ad-4df2-42e8-9855-5aa3ac664df2
md"""
The solution consists of two constants defined by the respective Dirichlet boundary conditions at the outer boundary.
"""

# ╔═╡ aad305a9-aac6-4aff-9f8e-08d6a2f756c8
md"""
### Mass action law reaction ``u_1 \leftrightharpoons u_2``

This is a rather general ansatz where we assume a backward-forward reaction between the two species meeting at the interface with reaction constants ``k_1`` and ``k_2``, respectively.

According to the mass action law, this translates to a reaction rate

``r(u_1,u_2)=k_1u_1 - k_2u_2``

and correspondingly

``f_1(u_1,u_2)=r``

``f_2(u_1,u_2)=-r`` 

Note, that ``f_i`` is monotonically increasing in ``u_i`` and monotonically decreasing in the respective other argument, leading to an M-Property of the overall discretization matrix.


Note that the "no reaction" case is just a special case where ``k_1,k_2=0``.
"""

# ╔═╡ d027ff24-3ad1-4528-b5cf-10814caf30db
begin
    const k1 = 0.1
    const k2 = 10
end

# ╔═╡ 1328b4bf-2d64-4b02-a910-1995da8be28b
function mal_reaction(f, u, node, data)
    if node.region == 3
        react = k1 * u[1] - k2 * u[2]
        f[1] = react
        f[2] = -react
    end
    return nothing
end

# ╔═╡ 610a0761-1c23-415d-a187-f7d93a1b7637
system2 = make_system(mal_reaction)

# ╔═╡ 87edce1f-df6d-4cd8-bce5-24fb666cd6b5
begin
    k1, k2
    U1, U2 = mysolve(system2)
    plot(U1, U2; title = "k1=$(k1), k2=$(k2)")
end

# ╔═╡ 2ecf2760-bb4a-4653-8c2d-f4d146e44cd4
md"""
The back reaction is 100 times stronger than the forward reaction. This means that species 2 is consumed, creating species 1.
"""

# ╔═╡ b82fc6b2-eee1-4a91-a115-61b86621f686
md"""
### Penalty enforcing continuity


Setting ``k_1,k_2`` to a large number leads to another special case of the above reaction - similar to the penalty method to implement the Dirichlet boundary conditions, this lets the reaction equation dominate, which in this case forces
``u_1-u_2=0`` at the interface, and thus continuity.
"""

# ╔═╡ 9eaea813-2628-47d0-9d36-54c367689142
function penalty_reaction(f, u, node, data)
    if node.region == 3
        react = 1.0e10 * (u[1] - u[2])
        f[1] = react
        f[2] = -react
    end
    return nothing
end

# ╔═╡ 817738c0-f1a3-4779-9075-7ea051a81e73
system3 = make_system(penalty_reaction);

# ╔═╡ 80019ef1-bf41-4a55-9262-613a2d20be1f
plot(mysolve(system3)...)

# ╔═╡ 2b474b25-d56e-4ba6-8c20-e07de9e803a3
md"""
### Penalty enforcing fixed jump

Instead of enforcing continuity, one can enforce a fixed jump.
"""

# ╔═╡ f99bf7c0-2246-4adf-9f6a-e5b7b3cbe0c0
const jump = 0.2

# ╔═╡ 7331db49-7ace-468e-87d8-56ab5d900905
function penalty_jump_reaction(f, u, node, data)
    if node.region == 3
        react = 1.0e10 * (u[1] - u[2] - jump)
        f[1] = react
        f[2] = -react
    end
    return nothing
end

# ╔═╡ 19b6dc1f-5e56-4487-be06-2ce90b030290
system3jump = make_system(penalty_jump_reaction);

# ╔═╡ f8bdf93d-7697-4d5c-92b5-976f8bcf605c
plot(mysolve(system3jump)...)

# ╔═╡ 31e00855-8906-4be5-8e69-e2d8d9539e04
md"""
### Interface recombination

Here, we implement an annihilation reaction ``u_1 + u_2 \to \emptyset``
According to the mass action law, this is implemented via

``r(u_1,u_2)=k_r u_1 u_2``

``f_1(u_1,u_2)=r``

``f_2(u_1,u_2)=r``



"""

# ╔═╡ f2490f99-04ca-4f42-af2a-53adae51ca68
const k_r = 1000

# ╔═╡ 39a0db1b-3a4e-4108-b43f-d4e578c92608
function recombination(f, u, node, data)
    if node.region == 3
        react = k_r * (u[1] * u[2])
        f[1] = react
        f[2] = react
    end
    return nothing
end;

# ╔═╡ 644149fb-2264-42bd-92c9-193ab07c08f6
system4 = make_system(recombination);

# ╔═╡ b479402f-ef00-4425-8b0f-45f2dae74d80
plot(mysolve(system4)...)

# ╔═╡ 661d3556-5520-4da3-bcdd-7882e4e36b1b
md"""
Bot species are consumed at the interface.
"""

# ╔═╡ ed068b51-92af-48d5-9230-debc178ec827
md"""
### Thin  conductive interface layer

Let us assume that the interface is of thickness $d$ which is however small with respect to ``\Omega`` that we want to derive an interface condition from the assumption of an exact continuous solution within the interface.

So let ``\Omega_I=(x_l,x_r)`` be  the interface region where
we have ``-\Delta u_I=0`` with values ``u_l``, ``u_r`` at the boundaries. 

Then we have for the flux  in the interface region, ``q_I=\nabla u = \frac1{d}(u_r - u_l)``

Continuity of fluxes then gives ``f_1=q_I`` and ``f_2=-q_I``.

Continuity of ``u`` gives ``u_{1,I}=u_l, u_{2,I}=u_r``
This gives

``r=q_I=\frac{1}{d}(u_1-u_{2})``

``f_1(u_1,v_1)=r``

``f_2(u_1,v_1)=-r``

and therefore another special case of the mass action law condition.
"""

# ╔═╡ a2d919a5-a395-40fb-8f93-db742f8a77c2
const d = 1

# ╔═╡ 58d8831b-ad66-4f77-a33a-933c15c46a52
function thinlayer(f, u, node, data)
    if node.region == 3
        react = (u[1] - u[2]) / d
        f[1] = react
        f[2] = -react
    end
    return nothing
end

# ╔═╡ 8c0b4ab5-09da-4d8f-b001-5e15f823423c
system5 = make_system(thinlayer);

# ╔═╡ d3d99b9c-ad18-4a04-bb3e-f17dd542f9f3
plot(mysolve(system5)...)

# ╔═╡ 5ca74233-6669-48c0-8842-a86449ac8e09
md"""
The solution looks very similar to the case of the jump condition, however here, the size of the jump is defined by the physics of the interface.
"""

# ╔═╡ eb9abf2e-372c-4f79-afbe-772b90eff9ad
md"""
## Multiple domains

From the above discussion it seems that discontinuous interface conditions can be formulated in a rather general way via linear or nonlinear robin boundary conditions for each of the adjacent discontinuous species. Technically, it is necessary to be able to access the adjacent bulk data.
"""

# ╔═╡ 9f3ae7b5-51b3-48bc-b4db-b7236ba30682
md"""
In order to streamline the handling of multiple interfaces,  we propose an API layer on top  of the species handling of VoronoiFVM. We call these "meta species" "quantities".
"""

# ╔═╡ d44407de-8c9c-42fa-b1a2-ae02b826eccc
N = 6

# ╔═╡ 2da8a5b1-b168-41d9-baa8-d65a4ef5c4c0
md"""
We define a grid with N=$(N) subregions
"""

# ╔═╡ ae268316-c058-4db8-9b71-57b0d9425274
begin
    XX = collect(0:0.1:1)
    local xcoord = XX
    for i in 1:(N - 1)
        xcoord = glue(xcoord, XX .+ i)
    end
    grid2 = simplexgrid(xcoord)
    for i in 1:N
        cellmask!(grid2, [i - 1], [i], i)
    end
    for i in 1:(N - 1)
        bfacemask!(grid2, [i], [i], i + 2)
    end
end

# ╔═╡ b53b9d28-4c25-4fb8-a3e4-599b0e385121
gridplot(grid2; legend = :lt, resolution = (600, 200))

# ╔═╡ e7ce7fd4-cfa4-4cc6-84a2-7e20ed2f4e5c
md"""
To work with quantities, we first introduce a new constructor call without the "physics" parameter:
"""

# ╔═╡ 29f36902-e355-4b02-b7b0-c4db12c47d33
system6 = VoronoiFVM.System(grid2)

# ╔═╡ 673e9320-ea30-4831-ad85-ba7936293ee2
md"""
First, we introduce a continuous quantity which we name "cspec". Note that the "species number" can be assigned automatically if not given explicitly.
"""

# ╔═╡ f35f419a-94dd-4051-a533-4b1ec9a4c9ec
const cspec = ContinuousQuantity(system6, 1:N; ispec = 1)

# ╔═╡ 9661e4fc-55e1-4c2c-a3ad-515cdac3b514
md"""
A discontinuous quantity can be introduced as well. by default, each reagion gets a new species number. This can be overwritten by the user. It is important that the speces numbers of neighboring regions differ.
"""

# ╔═╡ 90298676-fda7-4168-8a40-7ff53e7c761b
const dspec = DiscontinuousQuantity(system6, 1:N; regionspec = [2 + i % 2 for i in 1:N])

# ╔═╡ cebabf33-e769-47bd-b6f1-ddf525fea895
md"""
For both quantities, we define simple diffusion fluxes:
"""

# ╔═╡ 719f206a-5b9f-4d78-8778-1d89edb2bc4d
function flux2(f, u, edge, data)
    f[dspec] = u[dspec, 1] - u[dspec, 2]
    return f[cspec] = u[cspec, 1] - u[cspec, 2]
end

# ╔═╡ 1d7f442f-c057-4379-8a40-c6ce3646ad5c
md"""
Define a thin layer interface condition for `dspec` and an interface source for `cspec`.
"""

# ╔═╡ da41b22e-114d-4eee-81d0-73e6f3b45242
md"""
Add physics to the system, set dirichlet bc at both ends, and extract subgrids
for plotting (until there will be a plotting API for this...)
"""

# ╔═╡ 2867307e-1f46-4b62-8793-fa6668122bea
allsubgrids = subgrids(dspec, system6)

# ╔═╡ b8cd6ad1-d323-4888-bbd1-5deba5a5870d
const d1 = 0.1

# ╔═╡ 441a39a0-a7de-47db-8539-12dee30b8312
const q1 = 0.2

# ╔═╡ d6e1c6c7-060d-4c2f-8054-d8f33f54bd55
function breaction2(f, u, node, data)
    if node.region > 2
        react = (u[dspec, 1] - u[dspec, 2]) / d1
        f[dspec, 1] = react
        f[dspec, 2] = -react

        f[cspec] = -q1
    end
    return nothing
end

# ╔═╡ 59c83a22-a4cc-4b51-a1cc-5eb39588eacd
begin
    physics!(system6, VoronoiFVM.Physics(; flux = flux2, breaction = breaction2))

    ## Set boundary conditions
    boundary_dirichlet!(system6, dspec, 1, g_1)
    boundary_dirichlet!(system6, dspec, 2, g_2)
    boundary_dirichlet!(system6, cspec, 1, 0)
    boundary_dirichlet!(system6, cspec, 2, 0)

    # ensure that `solve` is called only after this cell
    # as mutating circumvents the reactivity of the notebook
    physics_ok = true
end;

# ╔═╡ de119a22-b695-4b4f-8e04-b7d68ec1e91b
if physics_ok
    sol6 = solve(system6; inival = 0.5)
end;

# ╔═╡ 83527778-76b2-4569-86c8-50dc6b48129f
function plot2(U, subgrids, system)
    dvws = VoronoiFVM.views(U, dspec, subgrids, system)
    cvws = VoronoiFVM.views(U, cspec, subgrids, system)
    vis = GridVisualizer(; resolution = (600, 300), legend = :rt)
    scalarplot!(
        vis,
        subgrids,
        grid2,
        dvws;
        flimits = (-0.5, 1.5),
        clear = false,
        color = :red,
        label = "discontinuous species"
    )
    scalarplot!(
        vis,
        subgrids,
        grid2,
        cvws;
        flimits = (-0.5, 1.5),
        clear = false,
        color = :green,
        label = "continuous species"
    )
    return reveal(vis)
end

# ╔═╡ d58407fe-dcd4-47bb-a65e-db5fedb58edc
if isdefined(Main, :PlutoRunner)
    plot2(sol6, allsubgrids, system6)
end

# ╔═╡ 8a435b96-4859-4452-b82a-e43a0c310a9a
md"""
## Testing
"""

# ╔═╡ 4d8e81c1-dbec-4379-9ab7-a585a369582d
if d1 == 0.1 && N == 6
    @test norm(system6, sol6, 2) ≈ 7.0215437706445245
end

# ╔═╡ 6ee90c4b-5ebf-48c4-a3b7-2efc32d16996
html"""<hr>"""

# ╔═╡ c344c0af-fb75-45f4-8977-45041a22b605
begin
    hrule() = html"""<hr>"""
    highlight(mdstring, color) = htl"""<blockquote style="padding: 10px; background-color: $(color);">$(mdstring)</blockquote>"""

    macro important_str(s)
        return :(highlight(Markdown.parse($s), "#ffcccc"))
    end
    macro definition_str(s)
        return :(highlight(Markdown.parse($s), "#ccccff"))
    end
    macro statement_str(s)
        return :(highlight(Markdown.parse($s), "#ccffcc"))
    end

    html"""
        <style>
         h1{background-color:#dddddd;  padding: 10px;}
         h2{background-color:#e7e7e7;  padding: 10px;}
         h3{background-color:#eeeeee;  padding: 10px;}
         h4{background-color:#f7f7f7;  padding: 10px;}

    	 pluto-log-dot-sizer  { max-width: 655px;}
         pluto-log-dot.Stdout { background: #002000;
    	                        color: #10f080;
                                border: 6px solid #b7b7b7;
                                min-width: 18em;
                                max-height: 300px;
                                width: 675px;
                                overflow: auto;
     	                       }
    	
        </style>
    """
end

# ╔═╡ b5a87200-eea5-4164-bfd5-dee1045a0464
html"""<hr>"""


# ╔═╡ Cell order:
# ╠═d5fcb04c-77f3-4a05-8497-a3eff38119a1
# ╠═18c423cc-18bf-41a0-a6e4-e30f91f39728
# ╟─327af4a8-74cc-4834-ab19-d3a1d6873982
# ╟─3fa189e4-9e1c-470c-bf26-15b631945d2d
# ╟─d5a0ee0d-959d-476b-b3c5-79b741059992
# ╟─f03ff283-c989-4b1a-b73e-2e616054e3db
# ╠═670c78c1-d0be-4362-975b-2c944620681f
# ╠═79193d53-9dfa-47b7-aed2-c7eb43769b5f
# ╠═4cb07222-587b-4a74-a444-43f5433d5b03
# ╟─02ec3c0b-6e68-462b-84df-931370cbdcac
# ╠═592429a1-108c-4e9e-8961-497c2c31f319
# ╠═1fa5d9b3-0558-4558-a5f1-f34f70c8d9a0
# ╟─5b539ad2-46d6-43b0-8b3e-f8a7e1ae0a6d
# ╠═6aabfbe1-de7d-49ba-8144-6d364b21b34f
# ╟─dbd27f10-9fd9-450d-9f78-89ea738d605b
# ╠═e5b432b2-875e-49a5-8e78-68f7f47f06c3
# ╠═8ae32292-df3e-4190-83a5-ec5ba529299e
# ╟─677853a7-0c43-4d3c-bc74-799535f95aeb
# ╠═139058a9-44a0-43d5-a377-4fc72927fa28
# ╟─de251b03-dd3a-4a44-9440-b7e654c32dac
# ╠═5b3a49d6-30ce-4d6e-94b4-a76645d7d8ce
# ╟─0735c061-68e1-429f-80f1-d8410989a91d
# ╠═467dc381-3b3d-4de7-a7f9-bfc51300832b
# ╟─fa4fcc0d-1d3a-45a2-8857-50536bbe39cc
# ╠═8f210696-fcf4-47bc-a5a2-c561ad7efcbd
# ╠═57e8515e-3be1-4478-af98-430501438ee7
# ╠═56136cd1-0c01-449d-9297-68924ac99ee7
# ╟─fa1293ad-4df2-42e8-9855-5aa3ac664df2
# ╟─aad305a9-aac6-4aff-9f8e-08d6a2f756c8
# ╠═1328b4bf-2d64-4b02-a910-1995da8be28b
# ╠═610a0761-1c23-415d-a187-f7d93a1b7637
# ╠═d027ff24-3ad1-4528-b5cf-10814caf30db
# ╟─87edce1f-df6d-4cd8-bce5-24fb666cd6b5
# ╟─2ecf2760-bb4a-4653-8c2d-f4d146e44cd4
# ╟─b82fc6b2-eee1-4a91-a115-61b86621f686
# ╠═9eaea813-2628-47d0-9d36-54c367689142
# ╠═817738c0-f1a3-4779-9075-7ea051a81e73
# ╠═80019ef1-bf41-4a55-9262-613a2d20be1f
# ╟─2b474b25-d56e-4ba6-8c20-e07de9e803a3
# ╠═f99bf7c0-2246-4adf-9f6a-e5b7b3cbe0c0
# ╠═7331db49-7ace-468e-87d8-56ab5d900905
# ╠═19b6dc1f-5e56-4487-be06-2ce90b030290
# ╠═f8bdf93d-7697-4d5c-92b5-976f8bcf605c
# ╟─31e00855-8906-4be5-8e69-e2d8d9539e04
# ╠═39a0db1b-3a4e-4108-b43f-d4e578c92608
# ╠═644149fb-2264-42bd-92c9-193ab07c08f6
# ╠═f2490f99-04ca-4f42-af2a-53adae51ca68
# ╠═b479402f-ef00-4425-8b0f-45f2dae74d80
# ╟─661d3556-5520-4da3-bcdd-7882e4e36b1b
# ╟─ed068b51-92af-48d5-9230-debc178ec827
# ╠═a2d919a5-a395-40fb-8f93-db742f8a77c2
# ╠═58d8831b-ad66-4f77-a33a-933c15c46a52
# ╠═8c0b4ab5-09da-4d8f-b001-5e15f823423c
# ╠═d3d99b9c-ad18-4a04-bb3e-f17dd542f9f3
# ╟─5ca74233-6669-48c0-8842-a86449ac8e09
# ╟─eb9abf2e-372c-4f79-afbe-772b90eff9ad
# ╟─9f3ae7b5-51b3-48bc-b4db-b7236ba30682
# ╟─2da8a5b1-b168-41d9-baa8-d65a4ef5c4c0
# ╠═d44407de-8c9c-42fa-b1a2-ae02b826eccc
# ╠═ae268316-c058-4db8-9b71-57b0d9425274
# ╠═b53b9d28-4c25-4fb8-a3e4-599b0e385121
# ╟─e7ce7fd4-cfa4-4cc6-84a2-7e20ed2f4e5c
# ╠═29f36902-e355-4b02-b7b0-c4db12c47d33
# ╟─673e9320-ea30-4831-ad85-ba7936293ee2
# ╠═f35f419a-94dd-4051-a533-4b1ec9a4c9ec
# ╟─9661e4fc-55e1-4c2c-a3ad-515cdac3b514
# ╠═90298676-fda7-4168-8a40-7ff53e7c761b
# ╟─cebabf33-e769-47bd-b6f1-ddf525fea895
# ╠═719f206a-5b9f-4d78-8778-1d89edb2bc4d
# ╟─1d7f442f-c057-4379-8a40-c6ce3646ad5c
# ╠═d6e1c6c7-060d-4c2f-8054-d8f33f54bd55
# ╟─da41b22e-114d-4eee-81d0-73e6f3b45242
# ╠═59c83a22-a4cc-4b51-a1cc-5eb39588eacd
# ╠═2867307e-1f46-4b62-8793-fa6668122bea
# ╠═de119a22-b695-4b4f-8e04-b7d68ec1e91b
# ╠═b8cd6ad1-d323-4888-bbd1-5deba5a5870d
# ╠═441a39a0-a7de-47db-8539-12dee30b8312
# ╠═83527778-76b2-4569-86c8-50dc6b48129f
# ╠═d58407fe-dcd4-47bb-a65e-db5fedb58edc
# ╟─8a435b96-4859-4452-b82a-e43a0c310a9a
# ╠═4d8e81c1-dbec-4379-9ab7-a585a369582d
# ╟─6ee90c4b-5ebf-48c4-a3b7-2efc32d16996
# ╟─c344c0af-fb75-45f4-8977-45041a22b605
# ╟─b5a87200-eea5-4164-bfd5-dee1045a0464
