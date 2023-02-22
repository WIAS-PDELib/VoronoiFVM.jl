### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 18c423cc-18bf-41a0-a6e4-e30f91f39728
begin
    if initialized
        using VoronoiFVM
        using ExtendableGrids
        using GridVisualize
        using PlutoVista
        using PlutoUI
        using HypertextLiteral
        using LinearAlgebra
        using LinearSolve
        using Test
        GridVisualize.default_plotter!(PlutoVista)
    end
end

# ╔═╡ 327af4a8-74cc-4834-ab19-d3a1d6873982
TableOfContents(; title = "", aside = false)

# ╔═╡ 3fa189e4-9e1c-470c-bf26-15b631945d2d
md"""
# Interface conditions in 1D
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



For the interior boundary (interface) conditons we set 


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
    grid = VoronoiFVM.Grid(X)

    bfacemask!(grid, [0.0], [0.0], 3)
    ## Material 1 left of 0
    cellmask!(grid, [-1.0], [0.0], 1)
    ## Material 2 right of 0
    cellmask!(grid, [0.0], [1.0], 2)
end;

# ╔═╡ 4cb07222-587b-4a74-a444-43f5433d5b03
gridplot(grid, legend = :rt, resolution = (600, 200))

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
function flux!(f, u, edge)
    if edge.region == 1
        f[1] = u[1, 1] - u[1, 2]
    end
    if edge.region == 2
        f[2] = u[2, 1] - u[2, 2]
    end
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
    physics = VoronoiFVM.Physics(flux = flux!, breaction = breaction)

    ## Create system
    sys = VoronoiFVM.System(grid, physics, unknown_storage = :sparse)

    ##  Enable species in their respective subregions
    enable_species!(sys, 1, [1])
    enable_species!(sys, 2, [2])

    ## Set boundary conditions
    for ispec = 1:2
        boundary_dirichlet!(sys, ispec, 1, g_1)
        boundary_dirichlet!(sys, ispec, 2, g_2)
    end
    sys
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
    U1, U2
end

# ╔═╡ 0735c061-68e1-429f-80f1-d8410989a91d
md"""
Plot the results
"""

# ╔═╡ 467dc381-3b3d-4de7-a7f9-bfc51300832b
function plot(U1, U2; title = "")
    vis = GridVisualizer(resolution = (600, 300))
    scalarplot!(
        vis,
        subgrid1,
        U1,
        clear = false,
        show = false,
        color = :green,
        label = "u1",
    )
    scalarplot!(
        vis,
        subgrid2,
        U2,
        clear = false,
        show = true,
        color = :blue,
        label = "u2",
        legend = :rt,
        title = title,
        flimits = (-0.5, 1.5),
    )
end

# ╔═╡ fa4fcc0d-1d3a-45a2-8857-50536bbe39cc
md"""
### No interface reaction

This means we set ``f_1(u_1,u_2)=0`` and ``f_2(u_1,u_2)=0``. 
"""

# ╔═╡ 8f210696-fcf4-47bc-a5a2-c561ad7efcbd
function noreaction(f, u, node) end

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

# ╔═╡ 1328b4bf-2d64-4b02-a910-1995da8be28b
function mal_reaction(f, u, node)
    if node.region == 3
        react = k1 * u[1] - k2 * u[2]
        f[1] = react
        f[2] = -react

    end
end

# ╔═╡ 610a0761-1c23-415d-a187-f7d93a1b7637
system2 = make_system(mal_reaction)

# ╔═╡ d027ff24-3ad1-4528-b5cf-10814caf30db
begin
    const k1 = 0.1
    const k2 = 10
end

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
function penalty_reaction(f, u, node)
    if node.region == 3
        react = 1.0e10 * (u[1] - u[2])
        f[1] = react
        f[2] = -react
    end
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
function penalty_jump_reaction(f, u, node)
    if node.region == 3
        react = 1.0e10 * (u[1] - u[2] - jump)
        f[1] = react
        f[2] = -react
    end
end

# ╔═╡ 19b6dc1f-5e56-4487-be06-2ce90b030290
system3jump = make_system(penalty_jump_reaction);

# ╔═╡ f8bdf93d-7697-4d5c-92b5-976f8bcf605c
plot(mysolve(system3jump)...)

# ╔═╡ 31e00855-8906-4be5-8e69-e2d8d9539e04
md"""
### Interface recombination

Here, we implement an annihilation reaction ``u_1 + u_2 \to \emptyset``
Acoording to the mass action law, this is implemented via

``r(u_1,u_2)=k_r u_1 u_2``

``f_1(u_1,u_2)=r``

``f_2(u_1,u_2)=r``



"""

# ╔═╡ 39a0db1b-3a4e-4108-b43f-d4e578c92608
function recombination(f, u, node)
    if node.region == 3
        react = k_r * (u[1] * u[2])
        f[1] = react
        f[2] = react
    end
end;

# ╔═╡ 644149fb-2264-42bd-92c9-193ab07c08f6
system4 = make_system(recombination);

# ╔═╡ f2490f99-04ca-4f42-af2a-53adae51ca68
const k_r = 1000

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
function thinlayer(f, u, node)
    if node.region == 3
        react = (u[1] - u[2]) / d
        f[1] = react
        f[2] = -react
    end
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

# ╔═╡ 2da8a5b1-b168-41d9-baa8-d65a4ef5c4c0
md"""
We define a grid with N=$(N) subregions
"""

# ╔═╡ d44407de-8c9c-42fa-b1a2-ae02b826eccc
N = 6

# ╔═╡ ae268316-c058-4db8-9b71-57b0d9425274
begin
    XX = collect(0:0.1:1)
    local xcoord = XX
    for i = 1:N-1
        xcoord = glue(xcoord, XX .+ i)
    end
    grid2 = simplexgrid(xcoord)
    for i = 1:N
        cellmask!(grid2, [i - 1], [i], i)
    end
    for i = 1:N-1
        bfacemask!(grid2, [i], [i], i + 2)
    end
end

# ╔═╡ b53b9d28-4c25-4fb8-a3e4-599b0e385121
gridplot(grid2, legend = :lt, resolution = (600, 200))

# ╔═╡ e7ce7fd4-cfa4-4cc6-84a2-7e20ed2f4e5c
md"""
To work with quantities, we first introduce a new constructor call without the "physics" parameter:
"""

# ╔═╡ 29f36902-e355-4b02-b7b0-c4db12c47d33
system6 = VoronoiFVM.System(grid2)

# ╔═╡ 673e9320-ea30-4831-ad85-ba7936293ee2
md"""
First, we introduce a continuous quantity which we name "cspec". Note that the "species number" can be assigned automatically if not given explicitely.
"""

# ╔═╡ f35f419a-94dd-4051-a533-4b1ec9a4c9ec
const cspec = ContinuousQuantity(system6, 1:N, ispec = 1)

# ╔═╡ 9661e4fc-55e1-4c2c-a3ad-515cdac3b514
md"""
A discontinuous quantity can be introduced as well. by default, each reagion gets a new species number. This can be overwritten by the user. It is important that the speces numbers of neighboring regions differ.
"""

# ╔═╡ 90298676-fda7-4168-8a40-7ff53e7c761b
const dspec = DiscontinuousQuantity(system6, 1:N; regionspec = [2 + i % 2 for i = 1:N])

# ╔═╡ cebabf33-e769-47bd-b6f1-ddf525fea895
md"""
For both quantities, we define simple diffusion fluxes:
"""

# ╔═╡ 719f206a-5b9f-4d78-8778-1d89edb2bc4d
function flux2(f, u, edge)
    f[dspec] = u[dspec, 1] - u[dspec, 2]
    f[cspec] = u[cspec, 1] - u[cspec, 2]
end

# ╔═╡ 1d7f442f-c057-4379-8a40-c6ce3646ad5c
md"""
Define a thin layer interface condition for `dspec` and an interface source for `cspec`.
"""

# ╔═╡ d6e1c6c7-060d-4c2f-8054-d8f33f54bd55
function breaction2(f, u, node)
    if node.region > 2
        react = (u[dspec, 1] - u[dspec, 2]) / d1
        f[dspec, 1] = react
        f[dspec, 2] = -react

        f[cspec] = -q1
    end
end

# ╔═╡ da41b22e-114d-4eee-81d0-73e6f3b45242
md"""
Add physics to the system, set dirichlet bc at both ends, and extract subgrids
for plotting (until there will be a plotting API for this...)
"""

# ╔═╡ 59c83a22-a4cc-4b51-a1cc-5eb39588eacd
begin
    physics!(system6, VoronoiFVM.Physics(flux = flux2, breaction = breaction2))


    ## Set boundary conditions
    boundary_dirichlet!(system6, dspec, 1, g_1)
    boundary_dirichlet!(system6, dspec, 2, g_2)
    boundary_dirichlet!(system6, cspec, 1, 0)
    boundary_dirichlet!(system6, cspec, 2, 0)
	
	# ensure that `solve` is called only after this cell
	# as mutating circumvents the reactivity of the notebook
	physics_ok=true 
end;

# ╔═╡ 2867307e-1f46-4b62-8793-fa6668122bea
allsubgrids = subgrids(dspec, system6)

# ╔═╡ de119a22-b695-4b4f-8e04-b7d68ec1e91b
if physics_ok
   sol6 = solve(system6, inival = 0.5)
end;

# ╔═╡ b8cd6ad1-d323-4888-bbd1-5deba5a5870d
const d1 = 0.1

# ╔═╡ 441a39a0-a7de-47db-8539-12dee30b8312
const q1 = 0.2

# ╔═╡ 83527778-76b2-4569-86c8-50dc6b48129f
function plot2(U, subgrids, system6)
    dvws = VoronoiFVM.views(U, dspec, allsubgrids, system6)
    cvws = VoronoiFVM.views(U, cspec, allsubgrids, system6)
    vis = GridVisualizer(resolution = (600, 300), legend = :rt)
    scalarplot!(
        vis,
        allsubgrids,
        grid2,
        dvws,
        flimits = (-0.5, 1.5),
        clear = false,
        color = :red,
        label = "discontinuous species",
    )
    scalarplot!(
        vis,
        allsubgrids,
        grid2,
        cvws,
        flimits = (-0.5, 1.5),
        clear = false,
        color = :green,
        label = "continuous species",
    )
    reveal(vis)
end

# ╔═╡ d58407fe-dcd4-47bb-a65e-db5fedb58edc
plot2(sol6, subgrids, system6)

# ╔═╡ 8a435b96-4859-4452-b82a-e43a0c310a9a
md"""
## Testing
"""

# ╔═╡ 4d8e81c1-dbec-4379-9ab7-a585a369582d
if d1 == 0.1 && N == 6
    @test norm(system6, sol6, 2) ≈ 7.0215437706445245
end

# ╔═╡ 335940d8-fe16-4256-abb0-0a25da14922f
md"""
## Appendix: Development + Styling
"""

# ╔═╡ 6ee90c4b-5ebf-48c4-a3b7-2efc32d16996
md"""
This notebook is also run during the automatic unit tests.
Next cell activates a development environment if the notebook is loaded from a checked out VoronoiFVM.jl. Otherwise, Pluto's built-in package manager is used.
"""


# ╔═╡ 12ab1a8c-943c-41e5-959f-4e0b956b2532
begin
    import Pkg as _Pkg
    developing = false
    if isfile(joinpath(@__DIR__, "..", "src", "VoronoiFVM.jl"))
        _Pkg.activate(@__DIR__)
        _Pkg.instantiate()
        using Revise
        developing = true
    end
    initialized = true
end;


# ╔═╡ bb45263a-9f7b-4322-8d62-de6dcf624c91
if developing
    @info "Developing VoronoiFVM at  $(pathof(VoronoiFVM))"
else
    @info "Loaded VoronoiFVM from  $(pathof(VoronoiFVM))"
end

# ╔═╡ c344c0af-fb75-45f4-8977-45041a22b605
begin
    hrule() = html"""<hr>"""
    highlight(mdstring, color) =
        htl"""<blockquote style="padding: 10px; background-color: $(color);">$(mdstring)</blockquote>"""

    macro important_str(s)
        :(highlight(Markdown.parse($s), "#ffcccc"))
    end
    macro definition_str(s)
        :(highlight(Markdown.parse($s), "#ccccff"))
    end
    macro statement_str(s)
        :(highlight(Markdown.parse($s), "#ccffcc"))
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


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ExtendableGrids = "cfc395e8-590f-11e8-1f13-43a2532b2fa8"
GridVisualize = "5eed8a63-0fb0-45eb-886d-8d5a387d12b8"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
LinearSolve = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PlutoVista = "646e1f28-b900-46d7-9d87-d554eb38a413"
Revise = "295af30f-e4ad-537b-8983-00126c2a3abe"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
VoronoiFVM = "82b139dc-5afc-11e9-35da-9b9bdfd336f3"

[compat]
ExtendableGrids = "~0.9.17"
GridVisualize = "~1.0.2"
HypertextLiteral = "~0.9.4"
LinearSolve = "~1.37.0"
PlutoUI = "~0.7.50"
PlutoVista = "~0.8.23"
Revise = "~3.5.1"
VoronoiFVM = "~0.19.5"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "95e834ebb7e9667419dbdb9ea1966e17e0e7a9e9"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Markdown", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "29e65c331f97db9189ef00a4c7aed8127c2fd2d4"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.27.10"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "0310e08cb19f5da31d08341c6120c047598f5b9c"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.5.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SnoopPrecompile", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "4d9946e51e24f5e509779e3e2c06281a733914c2"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.1.0"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SnoopPrecompile", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e5f08b5689b1aad068e01751889f2f615c7db36d"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.29"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "4aff5fa660eb95c2e0deb6bcdabe4d9a96bc4667"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "0.8.18"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "SnoopPrecompile", "SparseArrays"]
git-tree-sha1 = "ee75410471c18f40d57eb53840bc705a74566f23"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "0.17.16"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bijections]]
git-tree-sha1 = "fe4f8c5ee7f76f2198d5c2a06d3961c249cce7bd"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.4"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "0c5f81f47bbbcf4aea7b2959135713459170798b"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.5"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "Static"]
git-tree-sha1 = "2c144ddb46b552f72d7eafe7cc2f50746e41ea21"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.2"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "485193efd2176b88e6622a39a246f8c5b600e74e"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.6"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "70232f82ffaab9dc52585e0dd043b5e0c6b714f1"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.12"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "0e5c14c3bb8a61b3d53b2c0620570c332c8d0663"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.2.0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "9c209fb7536406834aa938fb149964b985de6c83"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.1"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random", "SnoopPrecompile"]
git-tree-sha1 = "aa3edc8f8dea6cbfa176ee12f7c2fc82f0608ed3"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.20.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "9441451ee712d1aec22edad62db1a9af3dc8d852"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.3"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "61fdd77467a5c3ad071ef8277ac6bd6af7dd4c04"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "02d2316b7ffceff992f3096ae48c7829a8aa0638"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.3"

[[deps.Configurations]]
deps = ["ExproniconLite", "OrderedCollections", "TOML"]
git-tree-sha1 = "62a7c76dbad02fdfdaa53608104edf760938c4ca"
uuid = "5218b696-f38b-4ac9-8b61-a12ec717816d"
version = "0.17.4"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "c5b6685d53f933c11404a3ae9822afe30d522494"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.12.2"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "9258430c176319dc882efa4088e2ff882a0cb1f1"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.81"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "988e2db482abeb69efc76ae8b6eba2e93805ee70"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.5.15"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPolynomials]]
deps = ["DataStructures", "Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "8b84876e31fa39479050e2d3395c4b3b210db8b0"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.4.6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.ElasticArrays]]
deps = ["Adapt"]
git-tree-sha1 = "e1c40d78de68e9a2be565f0202693a158ec9ad85"
uuid = "fdbdab4c-e67f-52f5-8c3f-e7b388dad3d4"
version = "1.2.11"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[deps.ExproniconLite]]
deps = ["Pkg", "TOML"]
git-tree-sha1 = "c2eb763acf6e13e75595e0737a07a0bec0ce2147"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.7.11"

[[deps.ExtendableGrids]]
deps = ["AbstractTrees", "Dates", "DocStringExtensions", "ElasticArrays", "InteractiveUtils", "LinearAlgebra", "Printf", "Random", "SparseArrays", "StaticArrays", "Test", "WriteVTK"]
git-tree-sha1 = "2921bf0ffab4c8b7eda6a36c7b06a0dde6df0137"
uuid = "cfc395e8-590f-11e8-1f13-43a2532b2fa8"
version = "0.9.17"

[[deps.ExtendableSparse]]
deps = ["DocStringExtensions", "ILUZero", "LinearAlgebra", "Printf", "Requires", "SparseArrays", "Sparspak", "SuiteSparse", "Test"]
git-tree-sha1 = "075dfc8c0049b676a6af6ce2d7e9e84f56371808"
uuid = "95c220a8-a1cf-11e9-0c77-dbfce5f500b3"
version = "0.9.6"

[[deps.Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[deps.FastLapackInterface]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c1293a93193f0ae94be7cf338d33e162c39d8788"
uuid = "29a986be-02c6-4525-aec4-84b980013641"
version = "1.2.9"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "7be5f99f7d15578798f338f5433b6c432ea8037b"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "d3ba08ab64bdfd27234d3f61956c966266757fe6"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.7"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "ed1b56934a2f7a65035976985da71b6a65b4f2cf"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.18.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "a69dd6db8a809f78846ff259298678f0d6212180"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.34"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.FuzzyCompletions]]
deps = ["REPL"]
git-tree-sha1 = "e16dd964b4dfaebcded16b2af32f05e235b354be"
uuid = "fb4132e2-a121-4a70-b8a1-d5b831dcdcc2"
version = "0.5.1"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "1cd7f0af1aa58abc02ea1d872953a97359cb87fa"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.4"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "e07a1b98ed72e3cdd02c6ceaab94b8a606faca40"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.2.1"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "fe9aea4ed3ec6afdfbeb5a4f39a2208909b162a6"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.5"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "1cf1d7dcb4bc32d7b4a5add4232db3750c27ecb4"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.8.0"

[[deps.GridVisualize]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "ElasticArrays", "ExtendableGrids", "GeometryBasics", "GridVisualizeTools", "HypertextLiteral", "LinearAlgebra", "Observables", "OrderedCollections", "PkgVersion", "Printf", "StaticArrays"]
git-tree-sha1 = "12165cfe9b04b67f0e349bf3bf370f08ea4cecbf"
uuid = "5eed8a63-0fb0-45eb-886d-8d5a387d12b8"
version = "1.0.2"

[[deps.GridVisualizeTools]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "StaticArraysCore"]
git-tree-sha1 = "7c892c426f8d03a180366411566d0f6ac1790f6c"
uuid = "5573ae12-3b76-41d9-b48c-81d0b6e61cc5"
version = "0.3.0"

[[deps.Groebner]]
deps = ["AbstractAlgebra", "Combinatorics", "Logging", "MultivariatePolynomials", "Primes", "Random"]
git-tree-sha1 = "47f0f03eddecd7ad59c42b1dd46d5f42916aff63"
uuid = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
version = "0.2.11"

[[deps.GroupsCore]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9e1a5e9f3b81ad6a5c613d181664a0efc6fe6dd7"
uuid = "d5909c97-4eac-4ecc-a3dc-fdd0858a4120"
version = "0.4.0"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "37e4657cd56b11abe3d10cd4a1ec5fbdb4180263"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.7.4"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "734fd90dd2f920a2f1921d5388dcebe805b262dc"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.14"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions", "Test"]
git-tree-sha1 = "709d864e3ed6e3545230601f94e11ebc65994641"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.11"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.ILUZero]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "b007cfc7f9bee9a958992d2301e9c5b63f332a90"
uuid = "88f59080-6952-5380-9ea5-54057fb9a43f"
version = "0.2.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "f366daebdfb079fd1fe4e3d560f99a0c892e15bc"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "3f91cd3f56ea48d4d2a75c2a65455c5fc74fa347"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.3"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "c3244ef42b7d4508c638339df1bdbf4353e144db"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.30"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "d9ae7a9081d9b1a3b2a5c1d3dac5e2fdaafbd538"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.22"

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "764164ed65c30738750965d55652db9c94c59bfe"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.4.0"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "dd90aacbfb622f898a97c2a4411ac49101ebab8a"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.9.0"

[[deps.KrylovKit]]
deps = ["ChainRulesCore", "GPUArraysCore", "LinearAlgebra", "Printf"]
git-tree-sha1 = "1a5e1d9941c783b0119897d29f2eb665d876ecf3"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.6.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "ForwardDiff", "LinearAlgebra", "MacroTools", "PreallocationTools", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "cd04158424635efd05ff38d5f55843397b7416a9"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.14.0"

[[deps.LambertW]]
git-tree-sha1 = "c5ffc834de5d61d00d2b0e18c96267cffc21f648"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.6"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "2422f47b34d4b127720a18f86fa7b1aa2e141f29"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.18"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "88b8f66b604da079a627b6fb2860d3704a6729a1"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.14"

[[deps.LazilyInitializedFields]]
git-tree-sha1 = "410fe4739a4b092f2ffe36fcb0dcc3ab12648ce1"
uuid = "0e77f7df-68c5-4e49-93ce-4cd80f5598bf"
version = "1.2.1"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.LightXML]]
deps = ["Libdl", "XML2_jll"]
git-tree-sha1 = "e129d9391168c677cd4800f5c0abb1ed8cb3794f"
uuid = "9c8b4983-aa76-5018-a973-4c85ecc9e179"
version = "0.9.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "DocStringExtensions", "FastLapackInterface", "GPUArraysCore", "IterativeSolvers", "KLU", "Krylov", "KrylovKit", "LinearAlgebra", "Preferences", "RecursiveFactorization", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "SnoopPrecompile", "SparseArrays", "Sparspak", "SuiteSparse", "UnPack"]
git-tree-sha1 = "d1fce810e9a4213607f0182cf25ffd6ce13e19b6"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "1.37.0"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "cedb76b37bc5a6c702ade66be44f831fa23c681e"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.0"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "ArrayInterfaceCore", "CPUSummary", "ChainRulesCore", "CloseOpenIntervals", "DocStringExtensions", "ForwardDiff", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "SIMDTypes", "SLEEFPirates", "SnoopPrecompile", "SpecialFunctions", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "d407ea0d7c354f5765914d0982c233328523c82f"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.151"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "60168780555f3e663c536500aa790b6368adc02a"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.3.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.MsgPack]]
deps = ["Serialization"]
git-tree-sha1 = "fc8c15ca848b902015bd4a745d350f02cf791c2a"
uuid = "99f44e22-a591-53d1-9472-aa23ef4bd671"
version = "1.2.0"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "eaa98afe2033ffc0629f9d0d83961d66a021dfcc"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.4.7"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "3295d296288ab1a0a2528feb424b854418acff57"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.2.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "6862738f9796b3edc1c09d0890afce4eca9e7e93"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.4"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "82d7c9e310fe55aa54996e6f7f94674e2a38fcb4"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.9"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "6503b77492fd7fcb9379bf73cd31035670e3c509"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.3.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9ff31d101d987eb9d66bd8b176ac7c277beccd09"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.20+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "cf494dca75a69712a72b80bc48f59dcf3dea63ec"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.16"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "6f4fbcd1ad45905a5dee3f4256fabb49aa2110c6"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.7"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f6cf8e7944e50901594838951729a1861e668cb8"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.2"

[[deps.Pluto]]
deps = ["Base64", "Configurations", "Dates", "Distributed", "FileWatching", "FuzzyCompletions", "HTTP", "HypertextLiteral", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "MsgPack", "Pkg", "PrecompileSignatures", "REPL", "RegistryInstances", "RelocatableFolders", "SnoopPrecompile", "Sockets", "TOML", "Tables", "URIs", "UUIDs"]
git-tree-sha1 = "0ee5bd226e5b95e2232229f7c4a97309ccd8158b"
uuid = "c3e4b0f8-55cb-11ea-2926-15256bba5781"
version = "0.19.22"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "5bb5129fdd62a2bbbe17c2756932259acf467386"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.50"

[[deps.PlutoVista]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "GridVisualizeTools", "HypertextLiteral", "Pluto", "UUIDs"]
git-tree-sha1 = "ce33a948b2dd496cf7db76ee466b726509e80d20"
uuid = "646e1f28-b900-46d7-9d87-d554eb38a413"
version = "0.8.23"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "0fe4e7c4d8ff4c70bfa507f0dd96fa161b115777"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.3"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "240d7170f5ffdb285f9427b92333c3463bf65bf6"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.1"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff", "Requires"]
git-tree-sha1 = "f739b1b3cc7b9949af3b35089931f2b58c289163"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.12"

[[deps.PrecompileSignatures]]
git-tree-sha1 = "18ef344185f25ee9d51d80e179f8dad33dc48eb1"
uuid = "91cefc8d-f054-46dc-8f8c-26e11d7c5411"
version = "3.0.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "311a2aa90a64076ea0fac2ad7492e914e6feeb81"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "786efa36b7eff813723c4849c90456609cf06661"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.8.1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RandomExtensions]]
deps = ["Random", "SparseArrays"]
git-tree-sha1 = "062986376ce6d394b23d5d90f01d81426113a3c9"
uuid = "fb686558-2515-59ef-acaa-46db3789a887"
version = "0.4.3"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "261dddd3b862bd2c940cf6ca4d1c8fe593e457c8"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.3"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "ChainRulesCore", "DocStringExtensions", "FillArrays", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "Requires", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "ZygoteRules"]
git-tree-sha1 = "3dcb2a98436389c0aac964428a5fa099118944de"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.38.0"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "SnoopPrecompile", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "9088515ad915c99026beb5436d0a09cd8c18163e"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.18"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegistryInstances]]
deps = ["LazilyInitializedFields", "Pkg", "TOML", "Tar"]
git-tree-sha1 = "ffd19052caf598b8653b99404058fce14828be51"
uuid = "2792f1a3-b283-48e8-9a74-f99dce5104f3"
version = "0.1.0"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "90cb983381a9dc7d3dff5fb2d1ee52cd59877412"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "50314d2ef65fce648975a8e80ae6d8409ebbf835"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.5"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "cda0aece8080e992f6370491b08ef3909d1c04e7"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.38"

[[deps.SciMLBase]]
deps = ["ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Preferences", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "33f031423eedc1f9e43f6112da6f13d5b49ea7da"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.86.2"

[[deps.SciMLOperators]]
deps = ["ArrayInterface", "DocStringExtensions", "Lazy", "LinearAlgebra", "Setfield", "SparseArrays", "StaticArraysCore", "Tricks"]
git-tree-sha1 = "8419114acbba861ac49e1ab2750bae5c5eda35c4"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.1.22"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SparseDiffTools]]
deps = ["Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "e19ac47477c9a8fcca06dab5e5471417d5d9d723"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.31.0"

[[deps.Sparspak]]
deps = ["Libdl", "LinearAlgebra", "Logging", "OffsetArrays", "Printf", "SparseArrays", "Test"]
git-tree-sha1 = "342cf4b449c299d8d1ceaf00b7a49f4fbc7940e7"
uuid = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
version = "0.3.9"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "d0435ba43ab5ad1cbb5f0d286ca4ba67029ed3ee"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.4"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "Requires", "SnoopPrecompile", "SparseArrays", "Static", "SuiteSparse"]
git-tree-sha1 = "5589ab073f8a244d2530b36478f53806f9106002"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.2.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "2d7d9e1ddadc8407ffd460e24218e37ef52dd9a3"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.16"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "ab6083f09b3e617e34a956b43e9d51b824206932"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.1.1"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "2842f1dbd12d59f2728ba79f4002cd6b61808f8b"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.4.8"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "b03a3b745aa49b566f128977a7dd1be8711c5e71"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.14"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+0"

[[deps.SymbolicIndexingInterface]]
deps = ["DocStringExtensions"]
git-tree-sha1 = "6b764c160547240d868be4e961a5037f47ad7379"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.2.1"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TimerOutputs", "Unityper"]
git-tree-sha1 = "ca0dbe8434ace322cea02fc8cce0dea8d5308e87"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "1.0.3"

[[deps.Symbolics]]
deps = ["ArrayInterface", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "Groebner", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Markdown", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TreeViews"]
git-tree-sha1 = "fce1fd0b13f860128c8b8aab0bab475eeeeb7994"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "5.1.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "c97f60dd4f2331e1a495527f80d242501d2f9865"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.1"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "f2fd3f288dfc6f507b0c3a2eb3bac009251e548b"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.22"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "94f38103c984f89cf77c402f2a68dbd870f8165f"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.11"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "31eedbc0b6d07c08a700e26d31298ac27ef330eb"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.19"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "d5f4ec8c22db63bd3ccb239f640e895cfde145aa"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.2"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "7bdcd1b36993026f91e61c3cc671c7127770be84"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.59"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.VoronoiFVM]]
deps = ["BandedMatrices", "CommonSolve", "DiffResults", "DocStringExtensions", "ExtendableGrids", "ExtendableSparse", "ForwardDiff", "GridVisualize", "InteractiveUtils", "JLD2", "LinearAlgebra", "LinearSolve", "Printf", "Random", "RecursiveArrayTools", "SnoopPrecompile", "SparseArrays", "SparseDiffTools", "StaticArrays", "Statistics", "SuiteSparse", "Symbolics", "Test"]
git-tree-sha1 = "8e4b60a100eb3bc1cd6e3931e7841e97d15728b5"
uuid = "82b139dc-5afc-11e9-35da-9b9bdfd336f3"
version = "0.19.5"

[[deps.WriteVTK]]
deps = ["Base64", "CodecZlib", "FillArrays", "LightXML", "TranscodingStreams"]
git-tree-sha1 = "49353f30da65f377cff0f934bb9f562a2c0441b9"
uuid = "64499a7a-5c06-52f2-abe2-ccb03c286192"
version = "1.17.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
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
# ╟─335940d8-fe16-4256-abb0-0a25da14922f
# ╟─6ee90c4b-5ebf-48c4-a3b7-2efc32d16996
# ╠═12ab1a8c-943c-41e5-959f-4e0b956b2532
# ╠═bb45263a-9f7b-4322-8d62-de6dcf624c91
# ╟─c344c0af-fb75-45f4-8977-45041a22b605
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002