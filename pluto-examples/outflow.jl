### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 9d1c6bef-b302-4c3c-b90c-1490b1d63568
# Convert this cell to markdown in order to enable Pluto's inbuilt package manager
if isdefined(Main, :PlutoRunner)
    using Pkg
    docsdir = joinpath(@__DIR__, "..", "docs")
    if isdir(docsdir)
        Pkg.activate(docsdir)
    end
end

# ╔═╡ 4fc7fda6-423b-48ea-8f86-6718a9050ee0
begin
    using Revise
    using Test
    using VoronoiFVM
    using GridVisualize
    using SimplexGridFactory, Triangulate
    using ExtendableGrids
    using PlutoUI, HypertextLiteral, UUIDs
    if isdefined(Main, :PlutoRunner)
        using CairoMakie
        CairoMakie.activate!(; type = "png", visible = false)
        GridVisualize.default_plotter!(CairoMakie)
    end
end

# ╔═╡ 7ea92fa4-272f-40b1-ac5e-a5f4808c8300
md"""
# Outflow boundary conditions

[Source](https://github.com/WIAS-PDELib/VoronoiFVM.jl/blob/master/pluto-examples/outflow.jl)
"""

# ╔═╡ 0d1b7d24-4dc7-4f1d-944f-e279e18f151c
md"""
We show how to implement outflow boundary conditions when the velocities at the boundary are calculated by another equation in the system. A typical case is solute transport in porous media where fluid flow is calculated by Darcy's law which defines the convective velocity in the solute transport equation.
"""

# ╔═╡ f951520e-3edc-4353-910c-7f6ad5678e6b
md"""
Regard the following system of equations in domain ``\Omega\subset \mathbb R^d``:
```math
\begin{aligned}
    \nabla \cdot \vec v &=0\\
	\vec v&=-k\nabla p\\
    \nabla \cdot \vec j &=0\\
   \vec j&= - D\nabla c - c\vec v 	
\end{aligned}
```
The variable `p` can be seen as a the pressure of a fluid in porous medium.
`c` is the concentration of a transported species.

We subdivide the boundary: ``\partial\Omega=Γ_{in}\cup Γ_{out}\cup Γ_{noflow}`` abs set 
```math
\begin{aligned}
   p=&1 	\quad &      c&=c_{in} & \text{on}\quad Γ_{in}\\
   p=&0 	\quad &      \vec j\cdot \vec n &= c\vec v \cdot \vec n & \text{on}\quad Γ_{out}\\
   \vec v\cdot \vec n &=0	\quad &      \vec j\cdot \vec n &= 0  & \text{on}\quad Γ_{noflow}\\
\end{aligned}
```
"""

# ╔═╡ 019239aa-7995-4bd7-9e9b-61a49e972f71
md"""
## Discretization data
"""

# ╔═╡ f2113479-4ace-4e47-8db3-e4c9cab95ecd
Base.@kwdef struct FlowTransportData
    k = 1.0
    v_in = 1.0
    c_in = 0.5
    D = 1.0
    Γ_in = 1
    Γ_out = 2
    ip = 1
    ic = 2
end

# ╔═╡ 97c582f7-4e95-41a6-85b5-4ec46dea1f48
X = 0:0.1:1

# ╔═╡ fdd26fa2-cef1-43bb-8b97-e8be39ba3a56
darcyvelo(u, data) = data.k * (u[data.ip, 1] - u[data.ip, 2])

# ╔═╡ dc2745d2-7616-4312-8321-af9086bc31f6
function flux(y, u, edge, data)
    vh = darcyvelo(u, data)
    y[data.ip] = vh

    bp, bm = fbernoulli_pm(vh / data.D)
    y[data.ic] = data.D * (bm * u[data.ic, 1] - bp * u[data.ic, 2])
    return nothing
end

# ╔═╡ addf62a3-9a71-49b8-b54a-441f65a01d01
function bcondition(y, u, bnode, data)
    boundary_neumann!(y, u, bnode; species = data.ip, region = data.Γ_in, value = data.v_in)
    boundary_dirichlet!(y, u, bnode; species = data.ip, region = data.Γ_out, value = 0)
    boundary_dirichlet!(y, u, bnode; species = data.ic, region = data.Γ_in, value = data.c_in)
    return nothing
end

# ╔═╡ e087d35c-47ca-47f4-ab20-32a7adf94f00
md"""
This function describes the outflow boundary condition.
It is called on edges (including interior ones) which have at least one ode
on one of the outflow boundaries. Within this function
`outflownode` can be used to identify
that node.  There is some ambiguity in the case that both nodes are outflow
nodes, in that case it is assumed that the contribution is zero. In the present case this is guaranteed by the constant Dirichlet boundary condition for the pressure.
"""

# ╔═╡ c5fb189b-e542-4313-bad2-d6fe64d70771
function boutflow(y, u, edge, data)
    y[data.ic] = -darcyvelo(u, data) * u[data.ic, outflownode(edge)]
    return nothing
end

# ╔═╡ 210aeda9-9c37-4278-8466-8d0a62347367
function flowtransportsystem(grid; kwargs...)
    data = FlowTransportData(; kwargs...)
    return VoronoiFVM.System(
        grid;
        flux,
        bcondition,
        boutflow,
        data,
        outflowboundaries = [data.Γ_out],
        species = [1, 2],
    )
end

# ╔═╡ 3e6b27b4-066f-475e-99f4-8eba666f9dc2
function checkinout(sys, sol)
    data = sys.physics.data
    tfact = TestFunctionFactory(sys)
    tf_in = testfunction(tfact, [data.Γ_out], [data.Γ_in])
    tf_out = testfunction(tfact, [data.Γ_in], [data.Γ_out])
    return (; in = integrate(sys, tf_in, sol), out = integrate(sys, tf_out, sol))
end

# ╔═╡ 3f273ed4-5522-4018-9551-b22a7db2e070
md"""
## 1D Case
"""

# ╔═╡ 00bb8b17-ce90-4b1c-9d8c-ef4a1294956c
grid = simplexgrid(X)

# ╔═╡ fa9fb64c-0600-426e-bb54-326b45d3e5de
sys1 = flowtransportsystem(grid);

# ╔═╡ 37c874cf-5fc7-47f8-99fe-1dafc01a5153
sol1 = solve(sys1; verbose = "n");

# ╔═╡ 1e3cac4e-1d53-4305-bbac-78b5ce6bc983
let
    vis = GridVisualizer(; size = (600, 300))
    scalarplot!(vis, grid, sol1[1, :])
    scalarplot!(vis, grid, sol1[2, :]; clear = false, color = :red)
    reveal(vis)
end

# ╔═╡ 30b8a9ef-6e12-4996-8c9a-35454555a835
t1 = checkinout(sys1, sol1)

# ╔═╡ 93ff29d0-249e-49fa-8375-6f41c494a6b4
@test t1.in ≈ -t1.out

# ╔═╡ 6b0fcdf7-5529-4eb3-b153-2e9e1c546cb9
@test maximum(sol1[2, :]) ≈ 0.5

# ╔═╡ b5e0928c-ab87-4d42-a6fc-3b04971c920b
@test minimum(sol1[2, :]) ≈ 0.5

# ╔═╡ fc29fa0f-c336-405b-9c9d-9e5052c2b453
md"""
## 2D Case
"""

# ╔═╡ 2067c00c-5056-475e-8ea1-261c8a44f1ad
begin
    g2 = simplexgrid(X, X)
    bfacemask!(g2, [1, 0.3], [1, 0.7], 5)
end

# ╔═╡ 0db7e6eb-2783-43f6-bcf4-fa1e16899b77
gridplot(g2; size = (300, 300))

# ╔═╡ 0aca8362-1fbe-4db6-9684-e0eacaa10fd0
sys2 = flowtransportsystem(g2; Γ_in = 4, Γ_out = 5);

# ╔═╡ 93ea6a39-4435-44d7-a2dd-091e8da15266
sol2 = solve(sys2; verbose = "n")

# ╔═╡ 2ebbdc56-806a-4590-a699-2a75eeb0b55c
let
    vis = GridVisualizer(; size = (700, 300), layout = (1, 2))
    scalarplot!(vis[1, 1], g2, sol2[1, :])
    scalarplot!(vis[1, 2], g2, sol2[2, :]; limits = (0, 1))
    reveal(vis)
end

# ╔═╡ c737b293-9595-4d1b-98a6-5b08ea6f8b90
t2 = checkinout(sys2, sol2)

# ╔═╡ 3a93a32b-45c8-4320-9679-2b094abac4e3
@test t2.in ≈ -t2.out

# ╔═╡ 6a872926-c62a-4425-a8f3-1ac8822b20c6
@test maximum(sol2[2, :]) ≈ 0.5

# ╔═╡ 839e214a-1e62-4c53-a879-fe5d199e1fce
@test minimum(sol2[2, :]) ≈ 0.5

# ╔═╡ 2f581185-b6ec-4e57-854d-a07653c484c2
md"""
## 3D Case
"""

# ╔═╡ 433f9542-2ea4-47c5-8cf5-aff0401dc19c
begin
    g3 = simplexgrid(X, X, X)
    bfacemask!(g3, [0.3, 0.3, 0], [0.7, 0.7, 0], 7)
end

# ╔═╡ f403351e-9174-471d-b4bd-6d5ce2f76ed5
gridplot(g3; size = (300, 300))

# ╔═╡ 9ad14ea9-359f-4cf4-a97c-cfa02855e5af
sys3 = flowtransportsystem(g3; Γ_in = 6, Γ_out = 7);

# ╔═╡ 46aab208-7ed7-4570-a0da-4e70665e9f25
sol3 = solve(sys3; verbose = "n")

# ╔═╡ d8bc8821-b9d7-4aae-8600-14db3af5e429
let
    vis = GridVisualizer(; size = (700, 300), layout = (1, 2))
    scalarplot!(vis[1, 1], g3, sol3[1, :])
    scalarplot!(vis[1, 2], g3, sol3[2, :]; limits = (0, 1))
    reveal(vis)
end

# ╔═╡ dd05ceb2-a42d-49fb-8e0a-f93ef3992ee3
t3 = checkinout(sys3, sol3)

# ╔═╡ 765bf153-5cc8-4ac6-885c-bb044b033969
@test t3.in ≈ -t3.out

# ╔═╡ e535aeca-f735-4dcf-88f3-9953221dcd01
@test maximum(sol3[2, :]) ≈ 0.5

# ╔═╡ 0a13160e-de6e-4625-931f-c61a8c8aa132
@test minimum(sol3[2, :]) ≈ 0.5

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781


# ╔═╡ f9b4d4dc-7def-409f-b40a-f4eba1163741
TableOfContents()

# ╔═╡ 7a93e9a8-8a2d-4b11-84ef-691706c0eb0f
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
        	/* Headers */
             h1{background-color:#dddddd;  padding: 10px;}
             h2{background-color:#e7e7e7;  padding: 10px;}
             h3{background-color:#eeeeee;  padding: 10px;}
             h4{background-color:#f7f7f7;  padding: 10px;}

    		/* "Terminal"  */
    	     pluto-log-dot-sizer  { max-width: 655px;}
             pluto-log-dot.Stdout { background: #002000;
    	                            color: #10f080;
                                    border: 6px solid #b7b7b7;
                                    min-width: 18em;
                                    max-height: 300px;
                                    width: 675px;
                                    overflow: auto;
     	                           }
            /* Standard cell width etc*/
    		main {
       			flex: 1;
    		    max-width: calc(700px + 25px + 6px); /* 700px + both paddings */
        		padding-top: 0px;
        		padding-bottom: 4rem;
        		padding-left: 25px;
        		padding-right: 6px;
        		align-content: center;
        		width: 100%;
    			}

           /* Cell width for slides*/
    		xmain {
    			margin: 0 auto;
    			max-width: 750px;
        		padding-left: max(20px, 3%);
        		padding-right: max(20px, 3%);
    	        }

    	
    	
        </style>
    """
end

# ╔═╡ 5beb3a0d-e57a-4aea-b7a0-59b8ce9ff5ce
hrule()

# ╔═╡ 86e1fffe-cf70-4d93-8131-153f0e98c94b
begin
    function floataside(text::Markdown.MD; top = 1)
        uuid = uuid1()
        return @htl(
            """
            		<style>


            		@media (min-width: calc(700px + 30px + 300px)) {
            			aside.plutoui-aside-wrapper-$(uuid) {

            	color: var(--pluto-output-color);
            	position:fixed;
            	right: 1rem;
            	top: $(top)px;
            	width: 400px;
            	padding: 10px;
            	border: 3px solid rgba(0, 0, 0, 0.15);
            	border-radius: 10px;
            	box-shadow: 0 0 11px 0px #00000010;
            	/* That is, viewport minus top minus Live Docs */
            	max-height: calc(100vh - 5rem - 56px);
            	overflow: auto;
            	z-index: 40;
            	background-color: var(--main-bg-color);
            	transition: transform 300ms cubic-bezier(0.18, 0.89, 0.45, 1.12);

            			}
            			aside.plutoui-aside-wrapper > div {
            #				width: 300px;
            			}
            		}
            		</style>

            		<aside class="plutoui-aside-wrapper-$(uuid)">
            		<div>
            		$(text)
            		</div>
            		</aside>

            		"""
        )
    end
    floataside(stuff; kwargs...) = floataside(md"""$(stuff)"""; kwargs...)
end;

# ╔═╡ 4652157f-3288-45b6-877e-17ed8dc7d18e
function slidemodeswitch()
    uuid = uuid1()
    return html"""
    	<script id=$(uuid)>

        const right = document.querySelector('button.changeslide.next')
        const left = document.querySelector('button.changeslide.prev')

        let fullScreen = false

        const func = (e) => {
            if (e.key == "F10") {
                e.preventDefault()
                window.present()
                if (fullScreen) {
                    document.exitFullscreen().then(() => fullScreen = false)
                } else {

                    document.documentElement.requestFullscreen().then(() => fullScreen = true)
                }
            }
            if (document.body.classList.contains('presentation')) {
             
                if (e.target.tagName == "TEXTAREA") return
                if (e.key == "PageUp") {
                    e.preventDefault()
                    left.click()
                    return
             }

                if (e.key == "PageDown") {
                    e.preventDefault()
                    right.click()
                    return
                }
                if (e.key == "Escape") {
                    window.present()
                    fullScreen = false
                document.exitFullscreen().catch(() => {return})
                }
            }
        }

        document.addEventListener('keydown',func)

        invalidation.then(() => {document.removeEventListener('keydown',func)})
    </script>
    """
end;

# ╔═╡ Cell order:
# ╠═9d1c6bef-b302-4c3c-b90c-1490b1d63568
# ╠═4fc7fda6-423b-48ea-8f86-6718a9050ee0
# ╟─7ea92fa4-272f-40b1-ac5e-a5f4808c8300
# ╟─0d1b7d24-4dc7-4f1d-944f-e279e18f151c
# ╟─f951520e-3edc-4353-910c-7f6ad5678e6b
# ╟─019239aa-7995-4bd7-9e9b-61a49e972f71
# ╠═f2113479-4ace-4e47-8db3-e4c9cab95ecd
# ╠═97c582f7-4e95-41a6-85b5-4ec46dea1f48
# ╠═fdd26fa2-cef1-43bb-8b97-e8be39ba3a56
# ╠═dc2745d2-7616-4312-8321-af9086bc31f6
# ╠═addf62a3-9a71-49b8-b54a-441f65a01d01
# ╟─e087d35c-47ca-47f4-ab20-32a7adf94f00
# ╠═c5fb189b-e542-4313-bad2-d6fe64d70771
# ╠═210aeda9-9c37-4278-8466-8d0a62347367
# ╠═3e6b27b4-066f-475e-99f4-8eba666f9dc2
# ╟─3f273ed4-5522-4018-9551-b22a7db2e070
# ╠═00bb8b17-ce90-4b1c-9d8c-ef4a1294956c
# ╠═fa9fb64c-0600-426e-bb54-326b45d3e5de
# ╠═37c874cf-5fc7-47f8-99fe-1dafc01a5153
# ╟─1e3cac4e-1d53-4305-bbac-78b5ce6bc983
# ╠═30b8a9ef-6e12-4996-8c9a-35454555a835
# ╠═93ff29d0-249e-49fa-8375-6f41c494a6b4
# ╠═6b0fcdf7-5529-4eb3-b153-2e9e1c546cb9
# ╠═b5e0928c-ab87-4d42-a6fc-3b04971c920b
# ╟─fc29fa0f-c336-405b-9c9d-9e5052c2b453
# ╠═2067c00c-5056-475e-8ea1-261c8a44f1ad
# ╠═0db7e6eb-2783-43f6-bcf4-fa1e16899b77
# ╠═0aca8362-1fbe-4db6-9684-e0eacaa10fd0
# ╠═93ea6a39-4435-44d7-a2dd-091e8da15266
# ╠═2ebbdc56-806a-4590-a699-2a75eeb0b55c
# ╠═c737b293-9595-4d1b-98a6-5b08ea6f8b90
# ╠═3a93a32b-45c8-4320-9679-2b094abac4e3
# ╠═6a872926-c62a-4425-a8f3-1ac8822b20c6
# ╠═839e214a-1e62-4c53-a879-fe5d199e1fce
# ╟─2f581185-b6ec-4e57-854d-a07653c484c2
# ╠═433f9542-2ea4-47c5-8cf5-aff0401dc19c
# ╠═f403351e-9174-471d-b4bd-6d5ce2f76ed5
# ╠═9ad14ea9-359f-4cf4-a97c-cfa02855e5af
# ╠═46aab208-7ed7-4570-a0da-4e70665e9f25
# ╠═d8bc8821-b9d7-4aae-8600-14db3af5e429
# ╠═dd05ceb2-a42d-49fb-8e0a-f93ef3992ee3
# ╠═765bf153-5cc8-4ac6-885c-bb044b033969
# ╠═e535aeca-f735-4dcf-88f3-9953221dcd01
# ╠═0a13160e-de6e-4625-931f-c61a8c8aa132
# ╟─5beb3a0d-e57a-4aea-b7a0-59b8ce9ff5ce
# ╟─60941eaa-1aea-11eb-1277-97b991548781
# ╟─f9b4d4dc-7def-409f-b40a-f4eba1163741
# ╟─7a93e9a8-8a2d-4b11-84ef-691706c0eb0f
# ╟─86e1fffe-cf70-4d93-8131-153f0e98c94b
# ╟─4652157f-3288-45b6-877e-17ed8dc7d18e
