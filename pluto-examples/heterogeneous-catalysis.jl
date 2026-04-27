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

# ╔═╡ 973f24a9-df94-46a9-8c54-34ea8e7e824c
# Convert this cell to markdown in order to enable Pluto's inbuilt package manager
if isdefined(Main, :PlutoRunner)
    using Pkg
    docsdir = joinpath(@__DIR__, "..", "docs")
    if isdir(docsdir)
        Pkg.activate(docsdir)
    end
    using Revise
end

# ╔═╡ e04d9162-e6ed-4e9f-86ce-51b5175f8103
begin
    using SciMLBase: ODEProblem, solve
    using OrdinaryDiffEqLowOrderRK: DP5
    using OrdinaryDiffEqRosenbrock: Rosenbrock23
    using OrdinaryDiffEqTsit5: Tsit5
    using Catalyst
    using VoronoiFVM: VoronoiFVM, enable_species!, enable_boundary_species!
    using VoronoiFVM: ramp, boundary_dirichlet!
    using ExtendableGrids: simplexgrid
    using GridVisualize: GridVisualize, GridVisualizer, reveal, scalarplot!, gridplot, available_kwargs
    doplots = isdefined(Main, :PlutoRunner)
    if doplots
        using Plots: Plots, plot, theme
        using PlotThemes
        Plots.gr()
        Plots.theme(:dark)
        GridVisualize.default_plotter!(Plots)
    end
    import PlutoUI
    using Test
    PlutoUI.TableOfContents(; depth = 4)
end

# ╔═╡ 784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
md"""
# Towards Heterogeneous Catalysis
"""

# ╔═╡ 32eca884-3ed2-4b03-b2d9-871e00e4c0b1
md"""
How to model and simulate heterogeneous catalysis with [Catalyst.jl](http://github.com/SciML/Catalyst.jl)
and [VoronoiFVM.jl](https://github.com/WIAS-PDELib/VoronoiFVM.jl).
"""

# ╔═╡ c8462f8e-9618-489e-9e02-7e2316c0e621
md"""
## Mass action kinetics
"""

# ╔═╡ 7d0c967a-c2b1-4cc5-aa83-95d070fd1e72
md"""
Sources: 
- General notation: [Horn/Jackson 1972](https://link.springer.com/content/pdf/10.1007/BF00251225.pdf)
- Textbook: [Érdi/Tóth 1989](https://www.google.de/books/edition/Mathematical_Models_of_Chemical_Reaction/iDu8AAAAIAAJ?hl=de&gbpv=1&dq=erdi+toth&pg=PA35&printsec=frontcover)
- [Wikipedia](https://en.wikipedia.org/wiki/Rate_equation#Stoichiometric_reaction_networks)
- [Catalyst.jl paper](https://doi.org/10.1371/journal.pcbi.1011530)
- [Catalyst.jl docs](https://docs.sciml.ai/Catalyst/stable/introduction_to_catalyst/math_models_intro/#math_models_in_catalyst)
"""

# ╔═╡ f972e726-573b-4910-95d1-66adb8bdd649
md"""
Assume ``j=1 … M`` reversible reactions of educt (substrate) species to product species
```math
	α_1^j S_1 + α_2^j S_2 + … + α_n^jS_n \underset{k_j^+}{\stackrel{k_j^-}{\longrightleftharpoons}} β_1^jS_1 + β_2^j S_2 + … + β_n^j S_n
```
or equivalently, 
```math 
∑_{i=1}^n α_{i}^j S_i  \underset{k_j^+}{\stackrel{k_j^-}{\longrightleftharpoons}} ∑_{i=1}^n β_i^j S_i 
```
The rate of these reactions depend on the concentrations
``[S_i]`` of the species. Within ``Catalyst.jl``, due to consistency with the derivation
from stochastic approaches, the  default "combinatoric rate law" is
```math
    r_j=k_j^+ ∏_{i=1}^n \frac{[S_i]^{α_i^j}}{α_i^j!} - k_j^- ∏_{i=1}^n \frac{[S_i]^{β_i^j}}{β_i^j!}  
```
while it appears that in most textbooks, the rate law is
```math
    r_j=k_j^+ ∏_{i=1}^n [S_i]^{α_i^j} - k_j^- ∏_{i=1}^n [S_i]^{β_i^j}.
```
See the corresponding remark in the [Catalyst.jl docs](https://docs.sciml.ai/Catalyst/stable/faqs/#faq_combinatoric_ratelaws)
and the [github issue](https://github.com/SciML/Catalyst.jl/issues/269). We will stick to the secobd version which can be
achieved by setting `combinatoric_ratelaws` to `false` at appropriate places.



Later in the project we will see that instead of the concentrations, we need to work with so called _activities_.
"""

# ╔═╡ a0e92e6c-f018-44ba-8f34-989818f5fc0c
md"""
The numbers ``σ_i^j=α_i^j-β_i^j`` are  the net stoichiometric coefficients of the system.
"""

# ╔═╡ 3bb9ce51-ed4a-485d-984e-6896e8083f6e
md"""
The rate differential equations then are (TODO: check this)
```math
	∂_t [S_i] + \sum_{j=1}^M \sigma_i^jr_j = 0
```
These assume that the reactions take place in a  continuously stirred tank reactor (CSTR)
which means  that we have complete mixing, 
and species concentrations are spatially constant, 
and we have just one concentration value for each species at a given point of time.
"""

# ╔═╡ e7d01b67-ef05-48f5-87b9-04944600e901
md"""
## Example 1: A ``\longrightleftharpoons`` B

```math
\begin{aligned}
   A&  \underset{k_1^+}{\stackrel{k_1^-}{\longrightleftharpoons}}  B\\
   r_1&= k_1^+ [A] - k_1^- [B]\\
   \partial_t[A] &= -r_1\\
   \partial_t[B] &= r_1
\end{aligned}   
```
"""

# ╔═╡ eb3379cb-1da1-4d43-8c76-4590dd90da04
md"""
### Solution via plain ODE problem using OrdinaryDiffEq.jl:
"""

# ╔═╡ b3d92a1b-a4cb-49bf-8bd7-912fd88652b9
md"""
Set the parameters such that the forward reaction is faster then the backward reaction:
"""

# ╔═╡ 7496dbdc-1892-4262-8a99-b0a19018685a
p1 = (k_p = 1, k_m = 0.1)

# ╔═╡ 8b1bd140-5fde-4cd5-ba54-06e4fa627b64
md"""
Define an ODE function describing the right hand side of the ODE System:
"""

# ╔═╡ 22b00d2e-f716-4a1c-a362-1a43481d5994
function example1(du, u, p, t)
    (; k_p, k_m) = p
    r1 = k_p * u[1] - k_m * u[2]
    du[1] = -r1
    return du[2] = r1
end

# ╔═╡ b68ec599-6e71-4522-ac1a-03dcd94b7668
md"""
Define some initial value:
"""

# ╔═╡ fba77ccc-f186-477c-bef6-e1d8ab76dfa4
u1_ini = [1.0, 0.0]

# ╔═╡ e2efd459-5b75-4722-9fc6-e553bee6484c
md"""
Create an solve the corresponding [ODEProblem](https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/)
"""

# ╔═╡ 87fb9170-ebd2-4168-92f9-b61e2dedd9e5
prob1 = ODEProblem(example1, u1_ini, (0, 10), p1)

# ╔═╡ 08605ca3-3b9f-4978-a64d-e66388238593
sol1 = solve(prob1, DP5())

# ╔═╡ 651e0888-a51f-4513-b31b-cde5099b45e9
doplots && plot(sol1; size = (600, 200))

# ╔═╡ 55d14ae8-add7-49d3-abf3-a02889e4c1d9
md"""
Mass conservation: adding the two reaction eqauations results in 
```math
	\partial_t ([A]+[B]) = 0,
```
therefore ``[A]+[B]`` must be constant:
"""

# ╔═╡ 20c7916f-f5e0-4916-a7ec-e4db80a7ef5d
all(s -> isapprox(s, sum(u1_ini)), sum(sol1; dims = 1))

# ╔═╡ c3acef74-6706-44a6-a9de-4af3f6002745
md"""
### `Catalyst.@reaction_network`
"""

# ╔═╡ 66b2995e-b510-4167-b8d3-5915ff7a8163
md"""
Catalyst.jl provides a convenient way to define a reaction network, and
the resulting reaction system. So we use this to build the same system:
"""

# ╔═╡ 6a246460-b389-4737-b264-cecfe19ecd4e
rn1 = @reaction_network rn1 begin
    @combinatoric_ratelaws false
    k_p, A --> B
    k_m, B --> A
end

# ╔═╡ f650db9b-7f4e-4777-9bc6-9445be108e7d
md"""
The corresponding ODE system is:
"""

# ╔═╡ e0203c0d-9866-407a-b680-5a0d3f0a1ecc
convert(ODESystem, rn1)

# ╔═╡ 4f1fa3bf-c3c9-49c3-8ddf-936fe0b38e70
md"""
Catalyst.jl adds a new method  to the ODEProblem constructor which allows to pass a reaction nerwork:
"""

# ╔═╡ 5baeeab7-570d-4094-b92a-027678165e15
prob1n = ODEProblem(rn1, u1_ini, (0, 10.0), Dict(pairs(p1)))

# ╔═╡ 3f0f2f0f-2931-403b-8c06-aabc0280a3ab
sol1n = solve(prob1n, Rosenbrock23())

# ╔═╡ fc602274-553c-4bc0-965b-dbf72ee2420f
doplots && plot(sol1n; idxs = [rn1.A, rn1.B, rn1.A + rn1.B], size = (600, 200))

# ╔═╡ 0f7fd013-c483-456f-b1e7-af80c2d05b0f
md"""
### Unraveling `@reaction_network`

Let us try to look behind the macro voodoo - this is necessary to build networks programmatically and is behind the translation from python to Julia in CatmapInterface.jl.
"""

# ╔═╡ 31ba9c9e-af50-4daa-a512-aa38c64eec7d
md"""
It is interesting if there is a "macro - less" way to define variables, parameters and species.
"""

# ╔═╡ 1da97997-69eb-442f-b215-fed8e0a0ea60
@variables t

# ╔═╡ 6c790fd1-285a-4e2c-a210-ccd8cbc04cb4
@parameters k_p k_m

# ╔═╡ 0b83cfb5-e043-4ff1-9a56-392cb08c3aa1
@species A(t) B(t)

# ╔═╡ 9a4c824d-3674-4247-bb51-3161556dad93
md"""
A reaction network can be combined from several reactions:
"""

# ╔═╡ 52bbcaa4-9fad-42d3-998b-74d227f1ac36
r1p = Reaction(k_p, [A], [B], [1], [1])

# ╔═╡ 7363a00e-9a29-4124-abde-effc8d8b3097
r1m = Reaction(k_m, [B], [A], [1], [1])

# ╔═╡ c6ea4d18-20e2-4ec1-b2ee-2342f840f31c
rn1x = complete(ReactionSystem([r1p, r1m], t; name = :example1))

# ╔═╡ eacf1e01-ed90-4370-a03e-1a40017d5b55
md"""
Once we are here, the rest remains the same.
"""

# ╔═╡ 720fa225-3572-4ced-a44e-61d7adbe083b
convert(ODESystem, rn1x)

# ╔═╡ b44ece71-09b2-4c47-a799-7c62bd3f0e28
prob1x = ODEProblem(rn1x, u1_ini, (0, 10.0), Dict(pairs(p1)))

# ╔═╡ cf56fe28-fb31-42c8-b3a2-003fe1aa11b9
sol1x = solve(prob1x, Rosenbrock23())

# ╔═╡ dd89b221-4360-4213-91d4-cad927ed32f7
doplots && plot(sol1x; size = (600, 200))

# ╔═╡ 1990a17b-38b0-46dc-a808-24ec51a50318
md"""
## Example 2: A + 2B ``\longrightleftharpoons`` AB_2
"""

# ╔═╡ 9ad624b2-4722-421e-9cae-8ffee1752ece
rn2 = @reaction_network rn2 begin
    @combinatoric_ratelaws false
    k_0A, ∅ --> A
    k_0B, ∅ --> B
    (k_1p, k_1m), A + 2B <--> AB_2
end

# ╔═╡ 4de8b3cd-9e9d-4952-a330-1aa467e44fa2
convert(ODESystem, rn2)

# ╔═╡ 7fb5e877-659c-45ac-92c8-ce0c0b35a515
p2 = (k_0A = 0.5, k_0B = 1, k_1p = 0.1, k_1m = 1.0e-5)

# ╔═╡ be9b595f-b201-42b5-8e1f-0ab25eb24df6
u2_ini = (A = 0, B = 0, AB_2 = 0)

# ╔═╡ 335055d9-5f33-4b27-97e4-98104d19f098
prob2 = ODEProblem(rn2, Dict(pairs(u2_ini)), (0, 20.0), Dict(pairs(p2)))

# ╔═╡ 3eb18099-fc92-4f60-8ebf-f9673c0bf0cd
sol2 = solve(prob2, Rosenbrock23())

# ╔═╡ e8efb099-439f-4f56-931a-1f8703ac3063
doplots && plot(sol2; legend = :topleft, size = (600, 300))

# ╔═╡ c9083f59-fd97-46f2-8020-fb7e86bcbca4
md"""
## Example 3: Catalysis for A + 2B ``\rightleftharpoons`` AB_2  
"""

# ╔═╡ 69b7757f-1856-4e46-bebd-b340acf5b757
md"""
The same reaction as in example 2, but now with a catalyst C.

The reaction between A and B takes places when A and B are bound (adsorbed) to the catalyst.  So we have adsorption reactions, reactions at the catalyst, and desorption. The overall of free and bound catalyst needs to be constant over time.
"""

# ╔═╡ 4e923331-67ef-4878-9bd9-0130030ce7e2
rn3 = @reaction_network rn3 begin
    @combinatoric_ratelaws false
    k_0A, ∅ --> A
    k_0B, ∅ --> B
    (k_1p, k_1m), A + C <--> CA
    (k_2p, k_2m), B + C <--> CB
    (k_3p, k_3m), CA + 2CB <--> CAB2 + 2C
    (k_4p, k_4m), CAB2 <--> AB2 + C
end

# ╔═╡ 7430aa42-1db3-445c-a547-961157de0d1a
convert(ODESystem, rn3)

# ╔═╡ 5cb92c3c-299f-4077-9ff8-012e24d3f9e8
p3 = (
    k_0A = 0.5, k_0B = 1,
    k_1p = 10, k_1m = 0.1,
    k_2p = 10, k_2m = 0.1,
    k_3p = 10, k_3m = 0.1,
    k_4p = 10, k_4m = 0.1,
)

# ╔═╡ 66500aa5-c2ec-4821-a926-658307614dd5
Cini = 40

# ╔═╡ 9f20138a-74fd-4468-8c4d-92109f39545a
u3_ini = (A = 0, B = 0, CA = 0, CB = 0, CAB2 = 0, AB2 = 0, C = Cini)

# ╔═╡ a7c5cb9d-99ca-4a8d-bf04-febdb7fa6a40
t3end = 200

# ╔═╡ c67f7488-bff8-42a8-bae2-c147b98f0ab7
prob3 = ODEProblem(rn3, Dict(pairs(u3_ini)), (0, t3end), Dict(pairs(p3)))

# ╔═╡ a56637ae-709d-4802-aac7-61c4b180081f
sol3 = solve(prob3, Rosenbrock23())

# ╔═╡ 6b6be455-f57b-477d-ba52-70c839d55cf7
doplots && plot(sol3; legend = :topleft, size = (600, 300))

# ╔═╡ 5da58dbd-b7b7-4151-b762-ccb11ecf0bd7
ctotal = rn3.C + rn3.CA + rn3.CB + rn3.CAB2

# ╔═╡ 227e8773-3835-42cf-8c79-08454e404149
sol3[ctotal]

# ╔═╡ 58a527cf-2873-4d2d-ae8f-ffbfe6b10e22
@test sol3[ctotal] ≈ fill(Cini, length(sol3))

# ╔═╡ 907c5c4d-e052-49ee-b96c-adcd37a9cf42
md"""
## Example 4: Heterogeneous catalysis
Heterogeneous catalysis assumes that the catalytic reaction takes place
at surface. This means that reacting species need to be transported 
towards or away from the surface, and one has to model coupled transport
and surface reaction.

Here we use VoronoiFVM.jl to model transport and Catalyst.jl to create
the surface reaction network.

### Problem specification
Assume ``\Omega=(0,1)`` where a catalytic reaction takes place at ``x=0``. 
We assume that the educts A, B, and the product AB2 are bulk species transported to the domain.
At ``x=1`` we set Dirichlet boundary conditions providing A,B and removing AB2.

A, B can adsorb  at the catalyst at ``x=0`` and  react to AB2 while adsorbed. 
The product desorbs and is released to the bulk. So we have 

- Mass transport in the interior of  ``\Omega``:
```math
\begin{aligned}
\partial_t c_A + \nabla \cdot D_A \nabla c_A &=0\\
\partial_t c_B + \nabla \cdot D_B \nabla c_B &=0\\
\partial_t c_{AB2} + \nabla \cdot D_{AB2} \nabla c_{AB2} &=0
\end{aligned}
```
- Coupled nonlinear robin boundary conditions at ``x=0``:
```math
\begin{aligned}
D_A\partial_n c_A + r_1 &= 0\\
D_B\partial_n c_A + r_2 &= 0\\
D_{AB2}\partial_n c_{AB2} - r_4 &= 0\\
\end{aligned}
```
- ``r_1, r_2`` and  ``r_4`` are asorption/desorption reactions:
```math
\begin{aligned}
  r_1&=k_{1p}c_A c_C - k_{1m}c_{CA}\\
  r_2&=k_{2p}c_B c_C - k_{2m}c_{CB}\\
  r_4&=k_{4p}c_{AB2}  - k_{4m}c_{C_C}c_{C_{AB2}}\\
\end{aligned}
```
"""

# ╔═╡ 1bb1ad0a-e89d-49ae-83ae-97b5e4e618b2
md"""
- The free catalyst sites C and  the catalyst coverages CA, CB, CAB2 behave according to:
"""

# ╔═╡ 1c1afb98-0e9e-40f5-b0d8-fa9d20bbe48a
md"""
- Dirichlet boundary conditions at  ``x=1`` :
```math
\begin{aligned}
	c_A&=1\\
        c_B&=1\\
	c_{AB2}&=0
\end{aligned}
```
"""

# ╔═╡ 90101f55-1d76-4ed2-ace9-77ad8b03bd66
md"""
Finally, we set all initial concentrations to zero besides of the catalyst concenration (number of catalyst sites) ``c_C|_{t=0}=C_0=1``.
"""

# ╔═╡ 78aeaad1-fa22-4276-9c73-e62be21133cb
md"""
### Implementation
"""

# ╔═╡ ad7673d8-6a8a-4827-b252-0f435da4d81f
md"""
#### Surface reaction network
"""

# ╔═╡ 96b7e52e-8d52-42af-b854-c0faf33040d8
md"""
Define a reaction network under the assumption that the supply of A and B
comes from the transport and does not need to be specified.
"""

# ╔═╡ d89fd2b1-78b6-43a9-974a-a97f081d2a2b
rnv = @reaction_network rnv begin
    @combinatoric_ratelaws false
    (k_1p, k_1m), A + C <--> CA
    (k_2p, k_2m), B + C <--> CB
    (k_3p, k_3m), CA + 2CB <--> CAB2 + 2C
    (k_4p, k_4m), CAB2 <--> AB2 + C
end

# ╔═╡ f0501240-75e9-425b-8690-81a8284aef28
odesys = convert(ODESystem, rnv)

# ╔═╡ 68a4d03f-a526-4d37-a45e-c951432e20e7
eqns = equations(odesys);

# ╔═╡ e82f039c-8de5-4f17-9187-3856c1a6aa0b
eqns[2]

# ╔═╡ 7be71fdc-304e-4801-ac53-b10bef0e7c51
eqns[3]

# ╔═╡ 0159929e-8086-4a06-8f62-ccfb3b86c308
eqns[5]

# ╔═╡ 96185f7f-99bc-437d-abdc-d8b21d54a9be
eqns[6]

# ╔═╡ 9d4a18cd-cc31-4606-9511-41e9cc6ca22a
md"""
For coupling with VoronoiFVM we need species numbers which need to correspond to the species in our network:
"""

# ╔═╡ 0c51422a-7dff-4d48-b575-5986ebbfb3a9
begin
    smap = speciesmap(rnv)
    const iA = smap[rnv.A]
    const iB = smap[rnv.B]
    const iC = smap[rnv.C]
    const iCA = smap[rnv.CA]
    const iCB = smap[rnv.CB]
    const iCAB2 = smap[rnv.CAB2]
    const iAB2 = smap[rnv.AB2]
end;

# ╔═╡ 45d917e4-829f-45c3-99a1-64e739ca171a
md"""
Grid:
"""

# ╔═╡ 7401a67f-5200-45aa-a075-377fab883ba5
grid = simplexgrid(0:0.01:1)

# ╔═╡ 8abc6bb4-afe6-42a8-84a6-a38c81d0bcc0
gridplot(grid, size = (600, 100))

# ╔═╡ 794d1fac-3c30-4c11-826b-50e198ddf874
md"""
The grid has two _boundary regions_: region 1 at x=0 and region 2 at x=1.
"""

# ╔═╡ 6c47e16f-9d7d-48db-baa8-73bcd5a0ca10
md"""
Reaction parameters:
"""

# ╔═╡ bcc319bc-6586-428d-ba56-7ea25b2deed9
pcat = (
    k_1p = 50, k_1m = 0.1,
    k_2p = 50, k_2m = 0.1,
    k_3p = 10, k_3m = 0.1,
    k_4p = 50, k_4m = 0.1,
)

# ╔═╡ 424a6774-9205-4da8-b8a5-943a50039934
md"""
Parameters for the VoronoiFVM system:
"""

# ╔═╡ 5b7ed6f1-0969-4bdd-9c58-676c2e3f3635
params = (
    D_A = 1.0,
    D_B = 1.0,
    D_AB2 = 1.0,
    pcat = pcat,
)

# ╔═╡ e08ce9e1-2dd7-4852-91e6-75745a67d69f
md"""
Initial values for the reaction network (needed only for the  definition of the ODE problem)
"""

# ╔═╡ bacdf186-607c-4771-b0ee-2e1d2b9987ba
C0 = 1.0

# ╔═╡ 170e7c33-9efb-4a4c-82f2-e18e5d0a7c13
uv_ini = (A = 0, B = 0, CA = 0, CB = 0, CAB2 = 0, AB2 = 0, C = C0)

# ╔═╡ 3717cddf-6bb4-4c71-84f0-1e6d17499481
tvend = 200.0

# ╔═╡ 6b570cde-a4c8-4dd0-86fb-4fe9daeadb9f
const probv = ODEProblem(rnv, Dict(pairs(uv_ini)), (0, tvend), Dict(pairs(pcat)))

# ╔═╡ 5719a8d0-0271-43a7-b217-96b4df9b249f
md"""
#### Callback functions for VoronoiFVM
"""

# ╔═╡ 6de1532c-c30c-482e-96eb-8a2dd98938e9
md"""
First, define flux and storage functions for the bulk process:
"""

# ╔═╡ dd667d53-abd5-4ab7-8158-94b9247cd575
function storage(y, u, node, p)
    y[iA] = u[iA]
    y[iB] = u[iB]
    return y[iAB2] = u[iAB2]
end

# ╔═╡ 249d1d2a-9dbd-4853-842b-a3efdb5d2012
function flux(y, u, edge, p)
    (; D_A, D_B, D_AB2) = p
    y[iA] = D_A * (u[iA, 1] - u[iA, 2])
    y[iB] = D_B * (u[iB, 1] - u[iB, 2])
    return y[iAB2] = D_A * (u[iAB2, 1] - u[iAB2, 2])
end

# ╔═╡ d5624271-dd74-480d-9b45-d276185243a9
md"""
Storage term for the surface reaction:
"""

# ╔═╡ a7872bef-409c-4c99-9ea7-d0792b3220e4
function bstorage(y, u, bnode, p)
    y[iC] = u[iC]
    y[iCA] = u[iCA]
    y[iCB] = u[iCB]
    return y[iCAB2] = u[iCAB2]
end

# ╔═╡ 3ca90bdc-ec53-4aea-b4ee-ecd2c348d262
md"""
Catalytic reaction. Here we use the right hand side function of the
ODE problem generated above. In VoronoiFVM, reaction term are a the 
left hand side, so we need to multiply by -1.

Note that we need to pass the parameter record as generated for the ODE problem
instead of `pcat`.
"""

# ╔═╡ 1cfb9878-453a-438b-9383-ada79b44ad1f
function catreaction(f, u, bnode, p)
    probv.f(f, u, probv.p, bnode.time)
    for i in 1:length(f)
        f[i] = -f[i]
    end
    return
end

# ╔═╡ 063f72b8-f021-49a3-bdb2-79aef69fdcff
md"""
Define the Dirichlet boundary condition  at x=1 (region 2):
"""

# ╔═╡ 37122f5b-cd43-4fc5-ac70-194a52cb5e34
function bulkbc(f, u, bnode, p)
    v = ramp(bnode.time; du = (0.0, 1.0), dt = (0.0, 0.01))
    boundary_dirichlet!(f, u, bnode; species = iA, value = v, region = 2)
    boundary_dirichlet!(f, u, bnode; species = iB, value = v, region = 2)
    return boundary_dirichlet!(f, u, bnode; species = iAB2, value = 0, region = 2)
end

# ╔═╡ 8374742e-d557-4fef-9ae1-2d8644bbcd6c
md"""
Dispatch the boundary conditions
"""

# ╔═╡ 461b603f-cee2-4731-825b-0a017fed162e
function breaction(f, u, bnode, p)
    return if bnode.region == 1
        catreaction(f, u, bnode, p)
    else
        bulkbc(f, u, bnode, p)
    end
end

# ╔═╡ 19de4726-63aa-47cd-8efa-8cc7df0ce8ad
md"""
#### Coupled transport-reaction system
"""

# ╔═╡ 76b7a30d-f957-456e-82b5-17f7bcab9d06
md"""
Define a VoronoiFVM system from grid, params and the callback functions
and enable the bulk and boundary species. 
`unknown_storage = :sparse` means that the solution is stored as a `nspecies x nnodes` 
sparse matrix in order to take into account that the surface species are non-existent
in the bulk. `unknown_storage = :dense` would store a full matrix and solve dummy equations
for the surface species values in the bulk.
"""

# ╔═╡ 41763843-a4fe-4adb-84af-c115a37bc265
begin
    sys = VoronoiFVM.System(
        grid; data = params,
        flux, breaction, bstorage, storage,
        unknown_storage = :sparse
    )
    enable_species!(sys, iA, [1])
    enable_species!(sys, iB, [1])
    enable_species!(sys, iAB2, [1])

    enable_boundary_species!(sys, iC, [1])
    enable_boundary_species!(sys, iCA, [1])
    enable_boundary_species!(sys, iCB, [1])
    enable_boundary_species!(sys, iCAB2, [1])
end;

# ╔═╡ 66a16772-2c9e-46ab-a608-78124fef0f48
md"""
Define an initial value for `sys`:
"""

# ╔═╡ da511ea3-15d9-4d26-86c7-f0453b6db2e9
begin
    u0 = VoronoiFVM.unknowns(sys; inival = 0)
    u0[iC, 1] = C0
end;

# ╔═╡ 4c5b1304-1ad1-4951-b0df-a86f8221f6d5
md"""
### Solution
"""

# ╔═╡ 02fd20eb-f6f0-4636-b805-52d6008552a5
md"""
Solve the time evolution
"""

# ╔═╡ bd958ed0-7500-4f04-b023-ba609cbc30e9
tsol = solve(sys; inival = u0, times = (1.0e-4, tvend));

# ╔═╡ 6b02ab1c-fde1-4f0b-b8f1-77e0f4a495f3
md"""
t: $(@bind log_t_plot PlutoUI.Slider(-4:0.1:log10(tvend), default = 0.4))
"""

# ╔═╡ 7f50d186-fad5-413b-bbc4-2d087474a72f
let
    t_plot = round(10^log_t_plot; sigdigits = 3)
    vis = GridVisualizer(; size = (600, 300), flimits = (0, 1), title = "Bulk concentrations: t=$t_plot", legend = :lt)
    sol = tsol(t_plot)
    scalarplot!(vis, grid, sol[iA, :]; color = :red, label = "A")
    scalarplot!(vis, grid, sol[iB, :]; color = :green, label = "B", clear = false)
    scalarplot!(vis, grid, sol[iAB2, :]; color = :blue, label = "AB2", clear = false)
    reveal(vis)
end

# ╔═╡ 8f8e2294-cb51-4649-a543-f6ef39384dad
Ctotalv = tsol[iC, 1, :] + tsol[iCA, 1, :] + tsol[iCB, 1, :] + tsol[iCAB2, 1, :]

# ╔═╡ a83e969e-ba25-4b6b-83e4-6febed1e8602
@test Ctotalv ≈ ones(length(tsol.t))

# ╔═╡ 2eeb7d1a-25bc-4c09-bc86-a998a7bf3ca7
let
    vis = GridVisualizer(;
        size = (600, 300),
        xlabel = "t",
        flimits = (0, 1), xlimits = (1.0e-3, tvend),
        legend = :lt, title = "Concentrations at x=0", xscale = :log10
    )
    t = tsol.t
    scalarplot!(vis, t, tsol[iA, 1, :]; color = :darkred, label = "A")
    scalarplot!(vis, t, tsol[iCA, 1, :]; color = :red, label = "CA")
    scalarplot!(vis, t, tsol[iB, 1, :]; color = :darkgreen, label = "B")
    scalarplot!(vis, t, tsol[iCB, 1, :]; color = :green, label = "CB")
    scalarplot!(vis, t, tsol[iAB2, 1, :]; color = :darkblue, label = "AB2")
    scalarplot!(vis, t, tsol[iCAB2, 1, :]; color = :blue, label = "CAB2")
    scalarplot!(vis, t, tsol[iC, 1, :] / C0; color = :orange, label = "C/C0")
    scalarplot!(vis, t, Ctotalv / C0; color = :darkorange, label = "Ctot/C0")

    reveal(vis)
end

# ╔═╡ Cell order:
# ╠═973f24a9-df94-46a9-8c54-34ea8e7e824c
# ╠═e04d9162-e6ed-4e9f-86ce-51b5175f8103
# ╟─784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
# ╟─32eca884-3ed2-4b03-b2d9-871e00e4c0b1
# ╟─c8462f8e-9618-489e-9e02-7e2316c0e621
# ╟─7d0c967a-c2b1-4cc5-aa83-95d070fd1e72
# ╟─f972e726-573b-4910-95d1-66adb8bdd649
# ╟─a0e92e6c-f018-44ba-8f34-989818f5fc0c
# ╟─3bb9ce51-ed4a-485d-984e-6896e8083f6e
# ╟─e7d01b67-ef05-48f5-87b9-04944600e901
# ╟─eb3379cb-1da1-4d43-8c76-4590dd90da04
# ╟─b3d92a1b-a4cb-49bf-8bd7-912fd88652b9
# ╠═7496dbdc-1892-4262-8a99-b0a19018685a
# ╟─8b1bd140-5fde-4cd5-ba54-06e4fa627b64
# ╠═22b00d2e-f716-4a1c-a362-1a43481d5994
# ╟─b68ec599-6e71-4522-ac1a-03dcd94b7668
# ╠═fba77ccc-f186-477c-bef6-e1d8ab76dfa4
# ╟─e2efd459-5b75-4722-9fc6-e553bee6484c
# ╠═87fb9170-ebd2-4168-92f9-b61e2dedd9e5
# ╠═08605ca3-3b9f-4978-a64d-e66388238593
# ╠═651e0888-a51f-4513-b31b-cde5099b45e9
# ╟─55d14ae8-add7-49d3-abf3-a02889e4c1d9
# ╠═20c7916f-f5e0-4916-a7ec-e4db80a7ef5d
# ╟─c3acef74-6706-44a6-a9de-4af3f6002745
# ╟─66b2995e-b510-4167-b8d3-5915ff7a8163
# ╠═6a246460-b389-4737-b264-cecfe19ecd4e
# ╟─f650db9b-7f4e-4777-9bc6-9445be108e7d
# ╠═e0203c0d-9866-407a-b680-5a0d3f0a1ecc
# ╟─4f1fa3bf-c3c9-49c3-8ddf-936fe0b38e70
# ╠═5baeeab7-570d-4094-b92a-027678165e15
# ╠═3f0f2f0f-2931-403b-8c06-aabc0280a3ab
# ╠═fc602274-553c-4bc0-965b-dbf72ee2420f
# ╟─0f7fd013-c483-456f-b1e7-af80c2d05b0f
# ╟─31ba9c9e-af50-4daa-a512-aa38c64eec7d
# ╠═1da97997-69eb-442f-b215-fed8e0a0ea60
# ╠═6c790fd1-285a-4e2c-a210-ccd8cbc04cb4
# ╠═0b83cfb5-e043-4ff1-9a56-392cb08c3aa1
# ╟─9a4c824d-3674-4247-bb51-3161556dad93
# ╠═52bbcaa4-9fad-42d3-998b-74d227f1ac36
# ╠═7363a00e-9a29-4124-abde-effc8d8b3097
# ╠═c6ea4d18-20e2-4ec1-b2ee-2342f840f31c
# ╟─eacf1e01-ed90-4370-a03e-1a40017d5b55
# ╠═720fa225-3572-4ced-a44e-61d7adbe083b
# ╠═b44ece71-09b2-4c47-a799-7c62bd3f0e28
# ╠═cf56fe28-fb31-42c8-b3a2-003fe1aa11b9
# ╠═dd89b221-4360-4213-91d4-cad927ed32f7
# ╟─1990a17b-38b0-46dc-a808-24ec51a50318
# ╠═9ad624b2-4722-421e-9cae-8ffee1752ece
# ╠═4de8b3cd-9e9d-4952-a330-1aa467e44fa2
# ╠═7fb5e877-659c-45ac-92c8-ce0c0b35a515
# ╠═be9b595f-b201-42b5-8e1f-0ab25eb24df6
# ╠═335055d9-5f33-4b27-97e4-98104d19f098
# ╠═3eb18099-fc92-4f60-8ebf-f9673c0bf0cd
# ╠═e8efb099-439f-4f56-931a-1f8703ac3063
# ╟─c9083f59-fd97-46f2-8020-fb7e86bcbca4
# ╟─69b7757f-1856-4e46-bebd-b340acf5b757
# ╠═4e923331-67ef-4878-9bd9-0130030ce7e2
# ╠═7430aa42-1db3-445c-a547-961157de0d1a
# ╠═5cb92c3c-299f-4077-9ff8-012e24d3f9e8
# ╠═66500aa5-c2ec-4821-a926-658307614dd5
# ╠═9f20138a-74fd-4468-8c4d-92109f39545a
# ╠═a7c5cb9d-99ca-4a8d-bf04-febdb7fa6a40
# ╠═c67f7488-bff8-42a8-bae2-c147b98f0ab7
# ╠═a56637ae-709d-4802-aac7-61c4b180081f
# ╠═6b6be455-f57b-477d-ba52-70c839d55cf7
# ╠═5da58dbd-b7b7-4151-b762-ccb11ecf0bd7
# ╠═227e8773-3835-42cf-8c79-08454e404149
# ╠═58a527cf-2873-4d2d-ae8f-ffbfe6b10e22
# ╟─907c5c4d-e052-49ee-b96c-adcd37a9cf42
# ╟─1bb1ad0a-e89d-49ae-83ae-97b5e4e618b2
# ╟─68a4d03f-a526-4d37-a45e-c951432e20e7
# ╟─e82f039c-8de5-4f17-9187-3856c1a6aa0b
# ╟─7be71fdc-304e-4801-ac53-b10bef0e7c51
# ╟─0159929e-8086-4a06-8f62-ccfb3b86c308
# ╟─96185f7f-99bc-437d-abdc-d8b21d54a9be
# ╟─1c1afb98-0e9e-40f5-b0d8-fa9d20bbe48a
# ╟─90101f55-1d76-4ed2-ace9-77ad8b03bd66
# ╟─78aeaad1-fa22-4276-9c73-e62be21133cb
# ╟─ad7673d8-6a8a-4827-b252-0f435da4d81f
# ╟─96b7e52e-8d52-42af-b854-c0faf33040d8
# ╠═d89fd2b1-78b6-43a9-974a-a97f081d2a2b
# ╠═f0501240-75e9-425b-8690-81a8284aef28
# ╟─9d4a18cd-cc31-4606-9511-41e9cc6ca22a
# ╠═0c51422a-7dff-4d48-b575-5986ebbfb3a9
# ╟─45d917e4-829f-45c3-99a1-64e739ca171a
# ╠═7401a67f-5200-45aa-a075-377fab883ba5
# ╠═8abc6bb4-afe6-42a8-84a6-a38c81d0bcc0
# ╟─794d1fac-3c30-4c11-826b-50e198ddf874
# ╟─6c47e16f-9d7d-48db-baa8-73bcd5a0ca10
# ╠═bcc319bc-6586-428d-ba56-7ea25b2deed9
# ╟─424a6774-9205-4da8-b8a5-943a50039934
# ╠═5b7ed6f1-0969-4bdd-9c58-676c2e3f3635
# ╟─e08ce9e1-2dd7-4852-91e6-75745a67d69f
# ╠═bacdf186-607c-4771-b0ee-2e1d2b9987ba
# ╠═170e7c33-9efb-4a4c-82f2-e18e5d0a7c13
# ╠═3717cddf-6bb4-4c71-84f0-1e6d17499481
# ╠═6b570cde-a4c8-4dd0-86fb-4fe9daeadb9f
# ╟─5719a8d0-0271-43a7-b217-96b4df9b249f
# ╟─6de1532c-c30c-482e-96eb-8a2dd98938e9
# ╠═dd667d53-abd5-4ab7-8158-94b9247cd575
# ╠═249d1d2a-9dbd-4853-842b-a3efdb5d2012
# ╟─d5624271-dd74-480d-9b45-d276185243a9
# ╠═a7872bef-409c-4c99-9ea7-d0792b3220e4
# ╟─3ca90bdc-ec53-4aea-b4ee-ecd2c348d262
# ╠═1cfb9878-453a-438b-9383-ada79b44ad1f
# ╟─063f72b8-f021-49a3-bdb2-79aef69fdcff
# ╠═37122f5b-cd43-4fc5-ac70-194a52cb5e34
# ╟─8374742e-d557-4fef-9ae1-2d8644bbcd6c
# ╠═461b603f-cee2-4731-825b-0a017fed162e
# ╟─19de4726-63aa-47cd-8efa-8cc7df0ce8ad
# ╟─76b7a30d-f957-456e-82b5-17f7bcab9d06
# ╠═41763843-a4fe-4adb-84af-c115a37bc265
# ╟─66a16772-2c9e-46ab-a608-78124fef0f48
# ╠═da511ea3-15d9-4d26-86c7-f0453b6db2e9
# ╟─4c5b1304-1ad1-4951-b0df-a86f8221f6d5
# ╟─02fd20eb-f6f0-4636-b805-52d6008552a5
# ╠═bd958ed0-7500-4f04-b023-ba609cbc30e9
# ╟─6b02ab1c-fde1-4f0b-b8f1-77e0f4a495f3
# ╠═7f50d186-fad5-413b-bbc4-2d087474a72f
# ╠═8f8e2294-cb51-4649-a543-f6ef39384dad
# ╠═a83e969e-ba25-4b6b-83e4-6febed1e8602
# ╠═2eeb7d1a-25bc-4c09-bc86-a998a7bf3ca7
