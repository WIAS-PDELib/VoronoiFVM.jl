### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ e00d0175-866e-4f0f-8121-49e7bbda6fb6
begin
	using Pkg
	inpluto=isdefined(Main,:PlutoRunner)
	developing=false	
    if inpluto && isfile(joinpath(@__DIR__,"..","src","VoronoiFVM.jl"))
		# We have to outsmart Plutos cell parser here.
	    eval(:(Pkg.activate(joinpath(@__DIR__))))
		using Revise
		eval(:(Pkg.develop(path=joinpath(@__DIR__,".."))))
	    developing=true
    end
end;

# ╔═╡ b285aca3-dee5-4b77-9276-537563e8643b
begin 
	using VoronoiFVM
    using ExtendableGrids
	using Test
	if inpluto
 	  using PlutoUI
	  using GridVisualize
	  using PlutoVista
	  GridVisualize.default_plotter!(PlutoVista)
   end
end;

# ╔═╡ 4ed0c302-26e4-468a-a40d-0e6406f802d0
md"""
# v0.14.0 API update 
"""

# ╔═╡ 7a104243-d3b9-421a-b494-5607c494b106
inpluto && TableOfContents(title="")

# ╔═╡ 3e6b4ffa-7b33-4b94-9fd6-75b030d5a115
md"""
Here we describe some updates for the API of `VoronoiFVM.jl`. These have been implemented mostly on top of the existing API, whose functionality is not affected.
"""

# ╔═╡ 8a4f336c-2016-453e-9a9f-beac66533fc0
md"""
## `VoronoiFVM.System` constructor
"""

# ╔═╡ 78e7c000-3a83-446a-b577-3a1809c664d2
md"""
### Implicit creation of physics

The `VoronoiFVM.Physics` struct almost never was used outside of the constructor of `VoronoiFVM.System`. Now it is possible to specify the flux functions directly in the system constructor. By default, it is als possible to set a list of species which are attached to all interior and boundary regions of the grid.
"""

# ╔═╡ c59876bd-0cb1-4157-9ba4-bdbede151a44
grid1=simplexgrid(0:0.1:1);

# ╔═╡ 90bbf212-c6c8-44f0-8132-4a98f094750e
multispecies_flux(y,u,edge)= @views @. y= u[:,1]-u[:,2]

# ╔═╡ adff41d1-9398-4a66-9a8e-e03809973fa6
function test_reaction(y,u,node)
	y[1]=u[1]
	y[2]=-u[1]
end

# ╔═╡ 5e6d83ab-65c7-4f33-b0a8-29cd5717b4d6
begin
   system1=VoronoiFVM.System(grid1,flux=multispecies_flux, 
   reaction=test_reaction,species=[1,2])
   boundary_dirichlet!(system1,species=1,region=1,value=1)
   boundary_dirichlet!(system1,species=2,region=2,value=0)
end;

# ╔═╡ f11c03a3-7899-42fd-a2da-a257715815dc
sol1=solve(system1);

# ╔═╡ 2754d4c8-bbc1-4283-8156-c660c33cd62d
if inpluto 
	let
	vis=GridVisualizer(resolution=(500,300),legend=:rt)
	scalarplot!(vis,grid1,sol1[1,:],color=:red,label="species1")
	scalarplot!(vis,grid1,sol1[2,:],color=:green,label="species2",clear=false)
	reveal(vis)
end
end

# ╔═╡ fbd75cf1-64e4-4f07-b54f-f90626f3f6ba
@test isapprox(sum(sol1),11.323894375033476, rtol=1.0e-14) 

# ╔═╡ b5f7e133-500d-4a27-8f78-11ff3582599c
md"""
### Boundary conditions as part of physics

This makes the API more consistent and opens an easy possibility to have
space and time dependent boundary conditions. One can specify them either in `breaction` or the synomymous `bcondition`.
"""

# ╔═╡ ec188c81-3374-4eed-9b7e-e22350886df2
function bcond2(y,u,bnode)
    boundary_neumann!(y,u,bnode,species=1,region=1,value=sin(bnode.time))
    boundary_dirichlet!(y,u,bnode,species=2,region=2,value=0)
end; 

# ╔═╡ c86e8a0f-299f-42ab-96f8-0cd62d50f196
system2=VoronoiFVM.System(grid1,flux=multispecies_flux, 
   reaction=test_reaction,species=[1,2],bcondition=bcond2,check_allocs=true);

# ╔═╡ b3d936fe-69ab-4013-b787-2f0b5410638a
sol2=solve(system2, times=(0,10), Δt_max=0.01);

# ╔═╡ 17749697-d5d8-4629-a625-e96590a5f0ac
if inpluto
   vis2=GridVisualizer(resolution=(500,300),limits=(-2,2),legend=:rt)
end

# ╔═╡ 0c916da5-2d6e-42df-ac4b-4a062f931ccd
if inpluto
	md"""time: $(@bind t2 Slider(0:0.01:10; default=5,show_value=true))"""
end

# ╔═╡ 783618f8-2470-4c7c-afc1-9800586625c1
if inpluto
    let
	s=sol2(t2)
	scalarplot!(vis2,grid1,s[1,:],color=:red,label="species1")
	scalarplot!(vis2,grid1,s[2,:],color=:green,label="species2",clear=false,title="time=$(t2)")
	reveal(vis2)
end
end

# ╔═╡ 4cbea340-9c02-4e69-8f5e-62bf45312bdd
@test isapprox(sum(sol2)/length(sol2),2.4908109508494247, rtol=1.0e-14) 

# ╔═╡ 1c18b5a0-cca6-46a1-bb9f-b3d65b8043c5
md"""
### Implicit creation of grid
"""

# ╔═╡ 47280b56-e5ec-4345-b4a1-7c3c92536b2e
md"""
By passing data for grid creation (one  to thre abstract vectors) instead a grid, a tensor product grid is implicitely created.
This example also demonstrates position dependent boundary values.
"""

# ╔═╡ a71086fa-4ec6-4842-a4e1-6a6b60441fc2
function bcond3(y,u,bnode)
	boundary_dirichlet!(y,u,bnode,region=4,value=bnode[2])
	boundary_dirichlet!(y,u,bnode,region=2,value=-bnode[2])
end;

# ╔═╡ a514231a-e465-4f05-ba4c-b20aa968d96f
system3=VoronoiFVM.System(-1:0.1:1, -1:0.1:1, flux=multispecies_flux,bcondition=bcond3,species=1);

# ╔═╡ d55f615c-d586-4ef7-adf9-5faf052b75ac
sol3=solve(system3);

# ╔═╡ c17e5104-4d3a-4d54-81c1-d7253245a8bb
@test isapprox(sum(sol3), 0.0, atol=1.0e-14)

# ╔═╡ 087ea16d-742e-4398-acf5-37248af1b5b4
md"""
## GridVisualize API extended to System
Instead of a grid, a system can be passed to `gridplot` and `scalarplot`.
"""

# ╔═╡ 34d465a5-7cc5-4348-b9ba-6d9381bb3a87
inpluto && scalarplot(system3,sol3,resolution=(300,300),levels=10,colormap=:hot)

# ╔═╡ ac48f5bd-fd1e-4aa7-a2c9-90f0f427143c
md"""
## Parameters of `solve`
"""

# ╔═╡ 69974c02-57e6-4eb5-acf4-b2d480fbd67d
md"""
The `solve` API has been simplified and made more Julian. All entries of `VoronoiFVM.NewtonControl` can be now passed as keyword arguments to `solve`.

Another new keyword argument is `inival` which allows to pass an initial value which by default is initialized to zero. Therefore we now can write `solve(system)` as we already have seen above.
"""

# ╔═╡ 1e12afcf-cf46-4672-9434-44fa8af95ef7
reaction4(y,u,bnode)= y[1]=-bnode[1]^2+u[1]^4;

# ╔═╡ 938ef63c-58c4-41a0-b3dd-4eb76987a4d7
bc4(args...)=boundary_dirichlet!(args...,value=0);

# ╔═╡ fe424654-f070-46a9-850a-738b1d4aca8f
system4=VoronoiFVM.System(-10:0.1:10,species=[1],
	reaction=reaction4,flux=multispecies_flux,bcondition=bc4);

# ╔═╡ 37fc8816-5ccd-436e-8335-ebb1218d8a35
sol4,history=solve(system4,log=true,damp_initial=0.001,damp_growth=3);

# ╔═╡ 6a256a29-f15f-4d82-8e84-7ceacb786715
inpluto && scalarplot(system4,sol4,resolution=(500,300),xlabel="x",ylabel="u",title="solution")

# ╔═╡ ea974ff6-6ed4-4320-972a-067a4372f547
inpluto && scalarplot(1:length(history),history,yscale=:log,resolution=(500,200),xlabel="steps",ylabel="residual",title="convergence history")

# ╔═╡ 5c2a3836-dc81-4950-88e5-7f603514b1c0
@test isapprox(sum(sol4),418.58515700568535, rtol=1.0e-14)

# ╔═╡ fc0245fe-1bf2-45a3-aa7c-9cce8d7eef37
html"""<hr><hr><hr>"""

# ╔═╡ f50c6497-cba3-491a-bedd-5f94f88f76fb
md"""
# Appendix: Tests & Development
"""

# ╔═╡ ad899a81-baab-4433-8b7f-1e5c3b18dae6
md"""
This notebook is also run during the automatic unit tests. In this case, all interactive elements and visualizations should be deactivated.
For this purposes, the next cell detects if the notebook is running under Pluto
and sets the `inpluto` flag accordingly.

Furthermore, the cell activates a development environment if the notebook is loaded from a checked out VoronoiFVM.jl. Otherwise, Pluto's built-in package manager is used.
"""

# ╔═╡ bdbe6513-70b1-4d97-a79c-71534caad2b7
if developing 
	md""" Developing VoronoiFVM at  $(pathof(VoronoiFVM))"""
else
	md""" Loaded VoronoiFVM from  $(pathof(VoronoiFVM))"""
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ExtendableGrids = "cfc395e8-590f-11e8-1f13-43a2532b2fa8"
GridVisualize = "5eed8a63-0fb0-45eb-886d-8d5a387d12b8"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PlutoVista = "646e1f28-b900-46d7-9d87-d554eb38a413"
Revise = "295af30f-e4ad-537b-8983-00126c2a3abe"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
VoronoiFVM = "82b139dc-5afc-11e9-35da-9b9bdfd336f3"

[compat]
ExtendableGrids = "~0.8.11"
GridVisualize = "~0.4.4"
PlutoUI = "~0.7.23"
PlutoVista = "~0.8.12"
Revise = "~3.2.1"
VoronoiFVM = "~0.13.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "abb72771fd8895a7ebd83d5632dc4b989b022b5b"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.2"

[[AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "265b06e2b1f6a216e0e8f183d28e4d354eab3220"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.2.1"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Cassette]]
git-tree-sha1 = "6ce3cd755d4130d43bab24ea5181e77b89b51839"
uuid = "7057c7e9-c182-5462-911a-8362d720325c"
version = "0.3.9"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "4c26b4e9e91ca528ea212927326ece5918a04b47"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.2"

[[ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "9aa8a5ebb6b5bf469a7e0e2b5202cf6f8c291104"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.0.6"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "a851fec56cb73cfdf43762999ec72eff5b86882a"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.15.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[DiffRules]]
deps = ["LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "9bc5dac3c8b6706b58ad5ce24cffd9861f07c94f"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.9.0"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[ElasticArrays]]
deps = ["Adapt"]
git-tree-sha1 = "a0fcc1bb3c9ceaf07e1d0529c9806ce94be6adf9"
uuid = "fdbdab4c-e67f-52f5-8c3f-e7b388dad3d4"
version = "1.2.9"

[[ExtendableGrids]]
deps = ["AbstractTrees", "Dates", "DocStringExtensions", "ElasticArrays", "InteractiveUtils", "LinearAlgebra", "Printf", "Random", "SparseArrays", "Test"]
git-tree-sha1 = "fbb0efd29f2ba5e25eeaf73b76257acfc1a28630"
uuid = "cfc395e8-590f-11e8-1f13-43a2532b2fa8"
version = "0.8.11"

[[ExtendableSparse]]
deps = ["DocStringExtensions", "LinearAlgebra", "Printf", "Requires", "SparseArrays", "SuiteSparse", "Test"]
git-tree-sha1 = "793bd32bb280668e80c476ce4a3d0f171c8122d5"
uuid = "95c220a8-a1cf-11e9-0c77-dbfce5f500b3"
version = "0.6.6"

[[FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "2db648b6712831ecb333eae76dbfd1c156ca13bb"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.11.2"

[[FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "8756f9935b7ccc9064c6eef0bff0ad643df733a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.7"

[[FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "8b3c09b56acaf3c0e581c66638b85c8650ee9dca"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.8.1"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "2b72a5624e289ee18256111657663721d59c143e"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.24"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "92243c07e786ea3458532e199eb3feee0e7e08eb"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.4.1"

[[GridVisualize]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "ElasticArrays", "ExtendableGrids", "GeometryBasics", "HypertextLiteral", "LinearAlgebra", "Observables", "OrderedCollections", "PkgVersion", "Printf", "StaticArrays"]
git-tree-sha1 = "eef34bda67d8ea865d7467a141b2c7d8eccd5592"
uuid = "5eed8a63-0fb0-45eb-886d-8d5a387d12b8"
version = "0.4.4"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLD2]]
deps = ["DataStructures", "FileIO", "MacroTools", "Mmap", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "46b7834ec8165c541b0b5d1c8ba63ec940723ffb"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.15"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "e273807f38074f033d94207a201e6e827d8417db"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.8.21"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "491a883c4fef1103077a7f648961adbf9c8dd933"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.1.2"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NaNMath]]
git-tree-sha1 = "f755f36b19a5116bb580de457cda0c140153f283"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.6"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Observables]]
git-tree-sha1 = "fe29afdef3d0c4a8286128d4e45cc50621b1e43d"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.4.0"

[[OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "d7fa6237da8004be601e19bd6666083056649918"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.3"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "a7a7e1a88853564e551e4eba8650f8c38df79b37"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.1.1"

[[PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "5152abbdab6488d5eec6a01029ca6697dff4ec8f"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.23"

[[PlutoVista]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "GridVisualize", "HypertextLiteral", "UUIDs"]
git-tree-sha1 = "2435d1d3e02db324414f268f30999b5c06a0d10f"
uuid = "646e1f28-b900-46d7-9d87-d554eb38a413"
version = "0.8.12"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[RecursiveArrayTools]]
deps = ["ArrayInterface", "ChainRulesCore", "DocStringExtensions", "FillArrays", "LinearAlgebra", "RecipesBase", "Requires", "StaticArrays", "Statistics", "ZygoteRules"]
git-tree-sha1 = "c944fa4adbb47be43376359811c0a14757bdc8a8"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.20.0"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "8f82019e525f4d5c669692772a6f4b0a58b06a6a"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.2.0"

[[Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "e55f4c73ec827f96cd52db0bc6916a3891c726b5"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.2.1"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SparseDiffTools]]
deps = ["Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "f33a0f6721b270cdf417f0c986e93c973e5913c8"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.19.1"

[[SparsityDetection]]
deps = ["Cassette", "DocStringExtensions", "LinearAlgebra", "SparseArrays", "SpecialFunctions"]
git-tree-sha1 = "9e182a311d169cb9fe0c6501aa252983215fe692"
uuid = "684fba80-ace3-11e9-3d08-3bc7ed6f96df"
version = "0.3.4"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f0bccf98e16759818ffc5d97ac3ebf87eb950150"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.8.1"

[[Static]]
deps = ["IfElse"]
git-tree-sha1 = "7f5a513baec6f122401abfc8e9c074fdac54f6c1"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.4.1"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "2ce41e0d042c60ecd131e9fb7154a3bfadbf50d3"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.3"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "bb1064c9a84c52e277f1096cf41434b675cd368b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.1"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[VoronoiFVM]]
deps = ["DiffResults", "DocStringExtensions", "ExtendableGrids", "ExtendableSparse", "ForwardDiff", "GridVisualize", "IterativeSolvers", "JLD2", "LinearAlgebra", "Printf", "RecursiveArrayTools", "SparseArrays", "SparseDiffTools", "SparsityDetection", "StaticArrays", "SuiteSparse", "Test"]
git-tree-sha1 = "3e994ef993900aa9adbe666c39bfcae7f5675a6e"
uuid = "82b139dc-5afc-11e9-35da-9b9bdfd336f3"
version = "0.13.4"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─4ed0c302-26e4-468a-a40d-0e6406f802d0
# ╟─7a104243-d3b9-421a-b494-5607c494b106
# ╟─3e6b4ffa-7b33-4b94-9fd6-75b030d5a115
# ╠═b285aca3-dee5-4b77-9276-537563e8643b
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
# ╟─0c916da5-2d6e-42df-ac4b-4a062f931ccd
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
# ╟─ea974ff6-6ed4-4320-972a-067a4372f547
# ╠═5c2a3836-dc81-4950-88e5-7f603514b1c0
# ╟─fc0245fe-1bf2-45a3-aa7c-9cce8d7eef37
# ╟─f50c6497-cba3-491a-bedd-5f94f88f76fb
# ╟─ad899a81-baab-4433-8b7f-1e5c3b18dae6
# ╠═e00d0175-866e-4f0f-8121-49e7bbda6fb6
# ╟─bdbe6513-70b1-4d97-a79c-71534caad2b7
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
