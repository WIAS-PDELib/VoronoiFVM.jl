### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin
    ENV["LC_NUMERIC"] = "C" # prevent pyplot from messing up string2float conversion
    using Triangulate
    using PlutoUI
    using PyPlot
    using Printf
    using SimplexGridFactory
    using GridVisualize
    using ExtendableGrids
    using VoronoiFVM
    GridVisualize.default_plotter!(PyPlot)
    PyPlot.svg(true)
end

# ╔═╡ cc49b312-7e21-49aa-b42d-32c184919aa2
TableOfContents(title = "", depth = 4, aside = false)

# ╔═╡ b61dce01-f84f-4e46-be76-8561bf0dddaa
md"""
# Handling of interfaces in VoronoiFVM

Let ``\Omega=(-0.5,1.5)\times (0,2)`` be subdivided into two
subregions, i.e. ``\Omega = \Omega_1\cap \Omega_2``. Let  ``\Gamma_{12}`` be the interface between them. Further, let ``Γ_1=\{-0.5\}\times (0,2)``,
        ``Γ_2=\{1.5\}\times (0,2)``, and ``Γ_N=\partial\Omega \setminus (Γ_1\cup Γ_2)``.
We assume that ``\Omega`` is triangulated with a Delaunay triangulation which conforms to interior
and exterior boundaries. The later property ensures, that subdomain boundaries are aligned
with triangle edges, and  that for $i=1,2$ the circumcenters
of triangles from ``\Omega_i`` are situated  in ``\bar\Omega_i``.
"""

# ╔═╡ f37e9de2-a681-45d8-9a0a-f53613f19319
md"""
Here we create such a domain using the Triangle mesh generator of J. Shewchuk via Julia's Triangulate.jl  package.
"""

# ╔═╡ b912c3ff-42aa-4484-8b84-f011c80382e9
begin
    const Γ_1 = 1
    const Γ_2 = 2
    const Γ_12 = 3
    const Γ_N = 4
    const Ω_1 = 1
    const Ω_2 = 2
end;

# ╔═╡ 8b806c22-3f34-11eb-0e03-4d3fa6461629
function tworegiongrid(; minangle = 20)
    triin = Triangulate.TriangulateIO()
    triin.pointlist = Cdouble[-0.5 0.0; 0.2 0.0; 1.5 0.0; 1.5 2.0; 0.8 2.0; -0.5 2.0]'
    triin.segmentlist = Cint[1 2; 2 3; 3 4; 4 5; 5 6; 6 1; 2 5]'
    triin.segmentmarkerlist = Cint[Γ_N, Γ_N, Γ_2, Γ_N, Γ_N, Γ_1, Γ_12]
    angle = @sprintf("%.15f", minangle)
    triin.regionlist = Cdouble[
        0.2 0.2 Ω_1 0.18
        0.8 0.2 Ω_2 0.18
    ]'
    return (triout, vorout) = triangulate("paAq$(angle)DQv", triin)
end

# ╔═╡ 40795d22-6091-4498-bd07-f716fb07682e
function plot_with_numbers(triout, vorout)
    PyPlot.clf()
    Triangulate.plot_triangulateio(PyPlot, triout, voronoi = vorout)
    PyPlot.scatter(triout.pointlist[1, :], triout.pointlist[2, :], color = :red, alpha = 1)
    dxy = [-0.03, 0.01]
    for ipoint in 1:size(triout.pointlist, 2)
        PyPlot.text(
            (dxy .+ triout.pointlist[:, ipoint])...,
            "$(ipoint)",
            ha = :right,
            color = :red,
        )
    end
    for ipoint in 1:size(vorout.pointlist, 2)
        PyPlot.text(
            (dxy .+ vorout.pointlist[:, ipoint])...,
            "$(ipoint)",
            ha = :right,
            color = :green,
        )
    end
    return PyPlot.gcf()
end

# ╔═╡ d55fdbfc-3f34-11eb-3f35-c3bcdb156288
triout, vorout = tworegiongrid()

# ╔═╡ 36ef21e0-078c-4287-a0fe-f942cd14dc08
plot_with_numbers(triout, vorout)

# ╔═╡ d6276724-6f55-426c-aba9-fb89d2c8166e
md"""
We see here collocation points in red, and triangle circumcenters aka Voronoi nodes in green.
"""

# ╔═╡ b4b335ef-de2b-4563-8f0f-cde4a04ac115
md"""
## Heterogeneous diffusion

### Problem formulations

```math
\begin{aligned}
-\nabla\cdot D\nabla u &=0, & \text{in}\quad Ω,\\
                     D\nabla u\cdot \vec n + αu &=g_1,& \text{on}\quad Γ_1, \\
                     D\nabla u\cdot \vec n + αu &=g_2,& \text{on}\quad Γ_2,\\
                     D\nabla u\cdot \vec n &=0,& \text{on}\quad Γ_N,
\end{aligned}
```

with ``D|_{\Omega_i}= D_i>0`` and ``α>0``. In fact, this is a sloppy version of a more precise formulation.

Let ``u_i`` be solely defined on ``\Omega_i, i = 1, 2`` and satisfy 
```math
\begin{aligned}
-\nabla\cdot D_i\nabla u_i &=0, & \text{in}\quad Ω_i,\\
                     D_i\nabla u_i\cdot \vec n_i+ α u_i &= g_i,& \text{on}\quad Γ_i, \\
                     D_i\nabla u_i\cdot \vec n_i &=0,& \text{on}\quad Γ_N\cap \partial Ω_i.
\end{aligned}
```
This means, two separate problems on ``\Omega_1, \Omega_2`` are considered. With the following interface conditions (in the sense of traces of functions) ``u_1, u_2`` are linked
```math
\begin{aligned}
    u_1 &= u_2,  & \text{on} \quad  Γ_{12},\\
 D_1\nabla u_1\cdot \vec n_1 &=-D_2\nabla u_2\cdot \vec n_2, & \text{on} \quad Γ_{12}.
\end{aligned}
```
From the later conditions we can directly follow the continuity of a solution ``u`` with ``u_{|\Omega_1} = u_1, u_{|\Omega_2} = u_2`` and the continuity of respective fluxes.

**Weak Formulation**

Let ``V=H^1(\Omega)`` and denote by ``v \in V`` a respective test function. Integrating over the domain, splitting the problem onto the subdomains, multiplying by the testfunction and using Greens first identity yields
```math
\begin{aligned}
  0 =&   \sum_{i=1}^2 \Bigl( \int_{Ω_i} D_i\nabla u \nabla v dx - \int_{Γ_i} (D_i\nabla u_i )\cdot \vec n_i v dγ  - \int_{Γ_N} (D_i\nabla u_i )\cdot \vec n_i v dγ\Bigr) \\
~
&  - \int_{Γ_{12}} \Bigl((D_1\nabla u_1)\cdot \vec n_1 + (D_2\nabla u_2 )\cdot \vec n_2 \Bigr)v dγ \quad ∀ v ∈ V
\end{aligned}
```
The integral over ``\Gamma_N`` vanishes due to the homogeneous Neumann boundary conditions, the integral over ``Γ_{12}`` likewise vanishes due to the interface conditions. Inserting the interface conditions on ``\Gamma_i, i = 1, 2`` brings us the following weak formulation:


Find ``u∈V`` such that
```math
    \sum_{i=1}^2 \int_{Ω_i} D_i\nabla u \nabla v dx + \sum_{i=1}^2 \int_{Γ_i} \alpha uv dγ =  \sum_{i=1}^2 \int_{Γ_i}  g_iv dγ \quad ∀ v ∈ V.
```


A P1 finite element discretization is straightforward if the interface is aligned with triangle edges.
"""

# ╔═╡ d95afcb5-68f9-4d2e-9111-06fc38fe949a
md"""
### Finite volume approach
We introduce the following notation: Let $\omega_K$ be a respective control volume. We denote by ``\sigma_{KL} = \partial \omega_K \cap \partial\omega_L`` the face between two adjacent control volumes. To each control volume we assign a collocation point ``x_K \in \omega_K``. The set of all neighbouring control volumes of ``\omega_K`` is given by ``\mathcal{N}(\omega_K)``. The Euclidean distance between two collocation points is given by ``h_{KL} = \Vert x_L - x_K \Vert_2``. Note that ``h_{KL}`` and ``\sigma_{KL}`` are orthogonal to each other.

In the following we consider the integral over a control volume ``\omega_K``
```math
0=-\int_{\omega_{K}} \nabla D\nabla u dx = -\int_{\partial\omega_{K}} (D\nabla u) \cdot \vec n_K d\gamma
```
and the respective discrete version of it.


#### Interior control volume, e.g. ``ω_{13}``
We can further rewrite the previous integral equation as a sum over shared faces, i.e.

```math
0= \sum_{\omega_L \in \mathcal{N}(\omega_K)} \int_{\sigma_{KL}} - D\nabla u \cdot \vec n_{KL} d\gamma \approx  \sum_{\omega_L \in \mathcal{N}(\omega_K)} \frac{|\sigma_{KL}|}{h_{KL}} D(u_K-u_L).
```

#### Boundary  control volume, e.g. ``ω_{8}``
We denote by ``\gamma_{K,j}=\partial\omega_K\cap Γ_j, j = 1, 2`` the faces of ``\omega_K`` which are located at the outer boundary. Thus, the reformulation of the integrated flux over a boundary control volume ``\omega_K`` for a fixed ``j`` reads

```math
\begin{aligned}
0=&\sum_{\substack{\omega_L \in \mathcal{N}(\omega_K)\\ \partial\omega_L \cap \Gamma_j = \emptyset}} \int_{\sigma_{KL}} -  D\nabla u \cdot \vec n_{KL} d\gamma + \int_{\gamma_{K, j}} - D\nabla u \cdot \vec n_{j}d\gamma  \\
  \approx&  \sum_{\substack{\omega_L \in \mathcal{N}(\omega_K)\\ \partial\omega_L \cap \Gamma_j = \emptyset}} \frac{|\sigma_{KL}|}{h_{KL}} D(u_K-u_L) +  |\gamma_{K, j}| (\alpha u_K-g_j)
\end{aligned}
```

#### Interface  control volume, e.g. ``ω_{7}``
This is the case we are interested to and where the magic happens. The control colume ``\omega_K`` is part of both domains ``\Omega_1, \Omega_2``. We can define ``\omega_{K,j} = \omega_K \cap \Omega_j, j = 1,2 `` as the part of ``\omega_K`` which is solely defined on one subdomain. Thus, we likewise introduce ``\sigma_{KL}^j = \partial \omega_{K,j} \cap \partial \omega_L`` as the respective face between a neighboring control volume and the sub-control volume. Lastly, we use the following notation for the face of ``\omega_{K,j}`` intersecting with the inner boundary ``Γ_{12}``, which is likewise the part of the interface intersecting ``\omega_K``: 
``\gamma_{K,12}=\omega_K\cap Γ_{12}``.


Integrating our PDE over this control volume ``\omega_K`` yields
```math
\begin{aligned}
    0  =& -\int_{\omega_K} -\nabla\cdot D\nabla u dx\\
    =& -\int_{\omega_{K,1}} -\nabla\cdot D_1\nabla u dx -\int_{\omega_{K,2}} -\nabla\cdot D_2\nabla udx\\
       =&  \sum_{\omega_L \in \mathcal{N}(\omega_{K,1})} \int_{\sigma_{KL}^1} - D_1\nabla u \cdot \vec n_{KL} d\gamma + \sum_{\omega_L \in \mathcal{N}(\omega_{K,2})} \int_{\sigma_{KL}^2} -D_2\nabla u \cdot \vec n_{KL} d\gamma \\
        & - \int_{\gamma_{K,12}}  D_1\nabla u \cdot \vec n_1 + D_2\nabla u \cdot \vec n_2 d\gamma \\
    \approx & \sum_{\omega_L \in \mathcal{N}(\omega_{K,1})} \frac{|\sigma_{KL}^1|}{h_{KL}} D_1(u_K-u_L) + \sum_{\omega_L \in \mathcal{N}(\omega_{K,2})} \frac{|\sigma_{KL}^2|}{h_{KL}} D_2(u_K-u_L)
\end{aligned}
```
The development used the flux continuity at the interface. Thus, interface control volumes are split into their
respective intersections with the adjacent domains, and integrals are performed separately. 

In practice, we attach the domain information to the triangles and perform the integration for each triangle separately.



"""

# ╔═╡ f9886dc9-5b1b-4c27-ac0d-5e65805e45bb
md"""
### Numerical solution with VoronoiFVM
"""

# ╔═╡ add94eb5-7859-474b-87fc-a4e4ab00e691
heterogrid = simplexgrid(
    triout.pointlist,
    triout.trianglelist,
    Int32.(triout.triangleattributelist[1, :]),
    triout.segmentlist,
    triout.segmentmarkerlist,
)

# ╔═╡ 427504f5-6c74-40d5-a4b9-9d2d7bd4f540
const ispec = 1

# ╔═╡ d908722f-a92d-49de-a3c2-aae33cc5a2e4
const D = [1, 100]

# ╔═╡ 6f108820-4d06-4fa3-b6ad-8f7454494670
flux(y, u, edge) = y[ispec] = D[edge.region] * (u[ispec, 1] - u[ispec, 2])

# ╔═╡ fca05141-71cb-4ecc-83b9-f50836dda100
md"""
Here we apply the Robin and Neumann boundary conditions in a `breaction`
"""

# ╔═╡ a8a69c7c-810a-4ea5-9a72-e0a2eea49782
const g = [0, 1, 0, 0]

# ╔═╡ bdc68399-19eb-48b5-ba69-5745a2541a17
const α = [100, 100, 0, 0]

# ╔═╡ ba273c41-f622-4f1a-9e1c-48b4143a6d43
breaction(y, u, bnode) = y[ispec] = α[bnode.region] * u[ispec] - g[bnode.region]

# ╔═╡ 92dadccc-44e9-49e5-ac94-a224dd8b61aa
begin
    heterophysics = VoronoiFVM.Physics(flux = flux, breaction = breaction)
    heterosys = VoronoiFVM.System(heterogrid, heterophysics)
    enable_species!(heterosys, ispec, [1, 2])
end

# ╔═╡ 40dc1b8d-9f11-4ae2-b566-793c3fea99fc
solution = solve(heterosys)

# ╔═╡ ead2425e-0109-4798-b79c-cbc8f76c1414
scalarplot(heterogrid, solution[1, :], resolution = (300, 300))

# ╔═╡ d5946f89-d8d7-4196-90a6-8e34723265c3
md"""
## Bulk reactions

In this case we have to evaluate integrals over reaction terms, which 
are split into the respective contributions from each subdomain which may
differ in the coefficients or in the physics.

TODO: elaborate
"""

# ╔═╡ b16061ae-8a08-4649-af8e-58ddf054996d
md"""
# Handling of a boundary flux

In this section, we will study the effect of the boundary flux. We consider the following test problem 
```math
- \nabla D \nabla u + k u = c
```
on two domains which shall in the limiting case coincide.
"""

# ╔═╡ 3d0eebe4-5eaa-4339-b6f6-bb5acd573033
md"""
The necessary parameters for the problem are defined here, likewise as the respective physics functions.
"""

# ╔═╡ d7910827-238a-4b71-ad44-7ce43ff1c0db
begin
    const dA_bulk = 2.0
    const dB_bulk = 2.0
    const dbB_interface = 5.0
    const k = 6.0
    const c = 1.0
end

# ╔═╡ 7b51cf5a-fa12-4756-a016-ed00cf90b1ba
md"""
## Domain A: Introducing a small interface layer

We introduce the height of this thin layer by following variable as a variable for the mesh refinement
"""


# ╔═╡ 3b429e0c-c07a-40df-ab0b-1f71d515ece3
n = 3

# ╔═╡ c770801b-d8f8-4095-98ff-e546f661e141
const thin_layer = 0.1 # just set 0.01, then you see agreement

# ╔═╡ 7d902e19-4ced-4685-ac5b-129ce101c9fe
begin
    function fluxA!(f, u, edge)
        if edge.region == 3
            f[1] = (dA_bulk + dbB_interface) * (u[1, 1] - u[1, 2])
        else
            f[1] = dA_bulk * (u[1, 1] - u[1, 2])
        end
        return nothing
    end

    function fluxB!(f, u, edge)
        f[1] = dB_bulk * (u[1, 1] - u[1, 2])
        return nothing
    end

    function reaction!(f, u, node)
        f[1] = k * u[1]
        return nothing
    end

    function source!(f, node)
        f[1] = c
        return nothing
    end

    function bfluxB!(f, u, bedge)
        if bedge.region == 3
            f[1] = thin_layer * dbB_interface * (u[1, 1] - u[1, 2])
        end
        return nothing
    end

    function bfluxA!(f, u, bedge)
        f[1] = 0.0
        return nothing
    end
end

# ╔═╡ fd59fa72-bc40-41ef-bdac-485b171da03e
begin
    h_pdopingA = 1.0
    h_intrinsic1A = 1.0 / 2
    h_thinlayerA = thin_layer
    h_intrinsic2A = 1.0 / 2 - thin_layer
    h_ndopingA = 1.0
    h_totalA = h_pdopingA + h_intrinsic1A + h_thinlayerA + h_intrinsic2A + h_ndopingA
    lengthA = 1.0

    regionAcceptorA = 1
    regionIntrinsic1A = 2
    regionThinLayerA = 3
    regionIntrinsic2A = 4
    regionDonorA = 5
    regionsA = [
        regionAcceptorA,
        regionIntrinsic1A,
        regionThinLayerA,
        regionIntrinsic2A,
        regionDonorA,
    ]

    # boundary region numbers
    bregionAcceptorA = 1
    bregionDonorA = 2

    coord_pA = collect(range(0.0, stop = h_pdopingA, length = n))
    coord_i_1A = collect(range(h_pdopingA, stop = h_pdopingA + h_intrinsic1A, length = n))
    coord_i_thinLA = collect(
        range(
            h_pdopingA + h_intrinsic1A,
            stop = h_pdopingA + h_intrinsic1A + h_thinlayerA,
            length = n,
        ),
    )
    coord_i_2A = collect(
        range(
            h_pdopingA + h_intrinsic1A + h_thinlayerA,
            stop = h_pdopingA + h_intrinsic1A + h_thinlayerA + h_intrinsic2A,
            length = n,
        ),
    )
    coord_nA = collect(range(h_totalA - h_ndopingA, stop = h_totalA, length = n))

    coordA = glue(coord_pA, coord_i_1A, tol = 1.0e-30)
    coordA = glue(coordA, coord_i_thinLA, tol = 1.0e-30)
    coordA = glue(coordA, coord_i_2A, tol = 1.0e-30)
    coordA = glue(coordA, coord_nA, tol = 1.0e-30)
    coord_lengthA = collect(range(0.0, stop = lengthA, length = 3 * n))

    gridA = simplexgrid(coordA, coord_lengthA)
    numberOfNodesA = size(gridA[Coordinates])[2]

    println("number of nodes in grid A are: ", numberOfNodesA)

    # specify inner regions
    cellmask!(gridA, [0.0, 0.0], [h_pdopingA, lengthA], regionAcceptorA, tol = 1.0e-30)
    cellmask!(
        gridA,
        [h_pdopingA, 0.0],
        [h_pdopingA + h_intrinsic1A, lengthA],
        regionIntrinsic1A,
        tol = 1.0e-30,
    )
    cellmask!(
        gridA,
        [h_pdopingA + h_intrinsic1A, 0.0],
        [h_pdopingA + h_intrinsic1A + h_thinlayerA, lengthA],
        regionThinLayerA,
        tol = 1.0e-30,
    )
    cellmask!(
        gridA,
        [h_pdopingA + h_intrinsic1A + h_thinlayerA, 0.0],
        [h_pdopingA + h_intrinsic1A + h_thinlayerA + h_intrinsic2A, lengthA],
        regionIntrinsic2A,
        tol = 1.0e-30,
    )
    cellmask!(
        gridA,
        [h_totalA - h_ndopingA, 0.0],
        [h_totalA, lengthA],
        regionDonorA,
        tol = 1.0e-30,
    )

    # specify outer regions
    # metal interfaces
    bfacemask!(
        gridA,
        [0.0, lengthA],
        [h_pdopingA, lengthA],
        bregionAcceptorA,
        tol = 1.0e-30,
    ) # BregionNumber = 1
    bfacemask!(gridA, [h_totalA, 0.0], [h_totalA, lengthA], bregionDonorA, tol = 1.0e-30) # BregionNumber = 2

    bfacemask!(
        gridA,
        [h_pdopingA + 1.0 / 2, 0.0],
        [h_pdopingA + 1.0 / 2, lengthA],
        3,
        tol = 1.0e-30,
    )

    # no flux interfaces [xmin, ymin], [xmax, ymax]
    # @Jürgen: If not defining this, then one of these boundaries is set to 1??
    bfacemask!(gridA, [0.0, 0.0], [0.0, lengthA], 4, tol = 1.0e-30)
    bfacemask!(gridA, [h_pdopingA, lengthA], [h_totalA, lengthA], 4, tol = 1.0e-30)
    bfacemask!(gridA, [0.0, 0.0], [h_totalA, 0.0], 4, tol = 1.0e-30)

    gridplot(gridA, Plotter = PyPlot, legend = :rt)
end

# ╔═╡ 0875ce0d-a22b-4e08-9b58-57673c986a67
md"""
## Domain B: Defining an interior interface

"""

# ╔═╡ a9d85817-df56-48ea-8aee-f700dc6310c8
begin
    h_pdopingB = 1.0
    h_intrinsicB = 1.0
    h_ndopingB = 1.0
    h_totalB = h_pdopingB + h_intrinsicB + h_ndopingB
    lengthB = 1.0

    # region numbers
    const regionAcceptorB = 1
    const regionIntrinsicB = 2
    const regionDonorB = 3
    const regionsB = [regionAcceptorB, regionIntrinsicB, regionDonorB]

    # boundary region numbers
    const bregionAcceptorB = 1
    const bregionDonorB = 2
    const binnerInterfaceB = 3

    coord_pB = collect(range(0.0, stop = h_pdopingB, length = n))
    coord_i1B = collect(range(h_pdopingB, stop = h_pdopingB + h_intrinsicB / 2, length = n))
    coord_i2B = collect(
        range(h_pdopingB + h_intrinsicB / 2, stop = h_pdopingB + h_intrinsicB, length = n),
    )
    coord_nB = collect(range(h_pdopingB + h_intrinsicB, stop = h_totalB, length = n))

    coordB = glue(coord_pB, coord_i1B)
    coordB = glue(coordB, coord_i2B)
    coordB = glue(coordB, coord_nB)
    coord_lengthB = collect(range(0.0, stop = lengthB, length = 4 * n))

    gridB = simplexgrid(coordB, coord_lengthB)
    numberOfNodesB = size(gridB[Coordinates])[2]

    println("number of nodes in grid B are: ", numberOfNodesB)

    # specify inner regions
    cellmask!(gridB, [0.0, 0.0], [h_pdopingB, lengthB], regionAcceptorB, tol = 1.0e-18)
    cellmask!(
        gridB,
        [h_pdopingB, 0.0],
        [h_pdopingB + h_intrinsicB, lengthB],
        regionIntrinsicB,
        tol = 1.0e-18,
    )
    cellmask!(
        gridB,
        [h_pdopingB + h_intrinsicB, 0.0],
        [h_totalB, lengthB],
        regionDonorB,
        tol = 1.0e-18,
    )

    # specify outer regions
    # metal interfaces
    bfacemask!(
        gridB,
        [0.0, lengthB],
        [h_pdopingB, lengthB],
        bregionAcceptorB,
        tol = 1.0e-30,
    )
    bfacemask!(gridB, [h_totalB, 0.0], [h_totalB, lengthB], bregionDonorB, tol = 1.0e-30)

    bfacemask!(
        gridB,
        [h_pdopingB + h_intrinsicB / 2, 0.0],
        [h_pdopingB + h_intrinsicB / 2, lengthB],
        binnerInterfaceB,
        tol = 1.0e-30,
    )

    # no flux interfaces [xmin, ymin], [xmax, ymax]
    # @Jürgen: If not defining this, then one of these boundaries is set to 1??
    bfacemask!(gridB, [0.0, 0.0], [0.0, lengthB], 4, tol = 1.0e-30)
    bfacemask!(gridB, [h_pdopingB, lengthB], [h_totalB, lengthB], 4, tol = 1.0e-30)
    bfacemask!(gridB, [0.0, 0.0], [h_totalB, 0.0], 4, tol = 1.0e-30)

    gridplot(gridB, Plotter = PyPlot, legend = :rt)
end

# ╔═╡ b7a40324-dea5-4fd8-9d7a-1031dceee305
md"""
## Solving the two test problems

"""

# ╔═╡ 905fd6fe-4858-4d85-aeba-7e57439531f5
begin
    const ispec_flux = 1

    sysA = VoronoiFVM.System(
        gridA,
        VoronoiFVM.Physics(
            flux = fluxA!,
            reaction = reaction!,
            source = source!,
            bflux = bfluxA!,
        ),
    )
    sysB = VoronoiFVM.System(
        gridB,
        VoronoiFVM.Physics(
            flux = fluxB!,
            reaction = reaction!,
            source = source!,
            bflux = bfluxB!,
        ),
    )

    # enable species in all regions
    enable_species!(sysA, ispec_flux, regionsA)
    enable_species!(sysB, ispec_flux, regionsB)

    # boundary conditions
    boundary_dirichlet!(sysA, ispec_flux, bregionAcceptorA, 1.0)
    boundary_dirichlet!(sysA, ispec_flux, bregionDonorA, 0.0)
    ####
    boundary_dirichlet!(sysB, ispec_flux, bregionAcceptorB, 1.0)
    boundary_dirichlet!(sysB, ispec_flux, bregionDonorB, 0.0)


    ## Stationary solution of both problems
    solA = solve(sysA)
    solB = solve(sysB)
end

# ╔═╡ 80cfa231-23b1-48a3-a5ad-ee0e122da746
begin
    p2 = GridVisualizer(;
        Plotter = PyPlot,
        dim = 2,
        layout = (1, 1),
        clear = false,
        resolution = (800, 500),
    )
    # this is for variable transformation, since we consider right outer boundary and want to transform to x-axis.
    function tran32!(a, b)
        return a[1] = b[2]
    end

    # note that if adjusting active_boundary to 3 or 4, then transform needs to be deleted.
    bgridA = subgrid(gridA, [3], boundary = true, transform = tran32!)
    bgridB = subgrid(gridB, [binnerInterfaceB], boundary = true, transform = tran32!)

    sol_boundA = view(solA[ispec, :], bgridA)
    sol_boundB = view(solB[ispec, :], bgridB)


    scalarplot!(
        p2[1, 1],
        bgridA,
        sol_boundA,
        show = true,
        clear = true,
        flimits = (0.1, 0.4),
        cellwise = true,
        color = :green,
        label = "Problem A ",
    )
    scalarplot!(
        p2[1, 1],
        bgridB,
        sol_boundB,
        show = true,
        clear = false,
        cellwise = true,
        color = :red,
        label = "Problem B ",
        legend = :lt,
    )
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ExtendableGrids = "cfc395e8-590f-11e8-1f13-43a2532b2fa8"
GridVisualize = "5eed8a63-0fb0-45eb-886d-8d5a387d12b8"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee"
SimplexGridFactory = "57bfcd06-606e-45d6-baf4-4ba06da0efd5"
Triangulate = "f7e6ffb2-c36d-4f8f-a77e-16e897189344"
VoronoiFVM = "82b139dc-5afc-11e9-35da-9b9bdfd336f3"

[compat]
ExtendableGrids = "~0.9.16"
GridVisualize = "~0.6.4"
PlutoUI = "~0.7.49"
PyPlot = "~2.11.0"
SimplexGridFactory = "~0.5.18"
Triangulate = "~2.2.0"
VoronoiFVM = "~0.19.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.3"
manifest_format = "2.0"
project_hash = "ee7394743850ff61a23f20a4f9bd275599b7492d"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Markdown", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "df23d15b1090a3332a09a7a51da45bd9f0a07f92"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.27.8"

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
deps = ["ArrayInterfaceCore", "Compat", "IfElse", "LinearAlgebra", "SnoopPrecompile", "Static"]
git-tree-sha1 = "dedc16cbdd1d32bead4617d27572f582216ccf23"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "6.0.25"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SnoopPrecompile", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e5f08b5689b1aad068e01751889f2f615c7db36d"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.29"

[[deps.ArrayInterfaceOffsetArrays]]
deps = ["ArrayInterface", "OffsetArrays", "Static"]
git-tree-sha1 = "3d1a9a01976971063b3930d1aed1d9c4af0817f8"
uuid = "015c0d05-e682-4f19-8f0a-679ce4c54826"
version = "0.1.7"

[[deps.ArrayInterfaceStaticArrays]]
deps = ["Adapt", "ArrayInterface", "ArrayInterfaceCore", "ArrayInterfaceStaticArraysCore", "LinearAlgebra", "Static", "StaticArrays"]
git-tree-sha1 = "f12dc65aef03d0a49650b20b2fdaf184928fd886"
uuid = "b0d46f97-bff5-4637-a19a-dd75974142cd"
version = "0.1.5"

[[deps.ArrayInterfaceStaticArraysCore]]
deps = ["Adapt", "ArrayInterfaceCore", "LinearAlgebra", "StaticArraysCore"]
git-tree-sha1 = "93c8ba53d8d26e124a5a8d4ec914c3a16e6a0970"
uuid = "dd5226c6-a4d4-4bc7-8575-46859f9c95b9"
version = "0.1.3"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "4aff5fa660eb95c2e0deb6bcdabe4d9a96bc4667"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "0.8.18"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "b277f5d9737fde7b5dce7f43faf49c80f542806c"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "0.17.10"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bijections]]
git-tree-sha1 = "fe4f8c5ee7f76f2198d5c2a06d3961c249cce7bd"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.4"

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
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "844b061c104c408b24537482469400af6075aae4"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.5"
weakdeps = ["ChainRulesCore"]

    [deps.ChangesOfVariables.extensions]
    ChainRulesCoreExt = "ChainRulesCore"

[[deps.CloseOpenIntervals]]
deps = ["ArrayInterface", "Static"]
git-tree-sha1 = "d61300b9895f129f4bd684b2aff97cf319b6c493"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.11"

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
version = "1.0.5+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "02d2316b7ffceff992f3096ae48c7829a8aa0638"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.3"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "6e47d11ea2776bc5627421d59cdcc1296c058071"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.7.0"

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
git-tree-sha1 = "74911ad88921455c6afcad1eefa12bd7b1724631"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.80"

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

[[deps.ExtendableGrids]]
deps = ["AbstractTrees", "Dates", "DocStringExtensions", "ElasticArrays", "InteractiveUtils", "LinearAlgebra", "Printf", "Random", "SparseArrays", "StaticArrays", "Test", "WriteVTK"]
git-tree-sha1 = "310b903a560b7b18f63486ff93da1ded9cae1f15"
uuid = "cfc395e8-590f-11e8-1f13-43a2532b2fa8"
version = "0.9.16"

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
git-tree-sha1 = "7fbaf9f73cd4c8561702ea9b16acf3f99d913fe4"
uuid = "29a986be-02c6-4525-aec4-84b980013641"
version = "1.2.8"

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
deps = ["ArrayInterfaceCore", "LinearAlgebra", "Requires", "Setfield", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "04ed1f0029b6b3af88343e439b995141cb0d0b8d"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.17.0"

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
git-tree-sha1 = "a5e6e7f12607e90d71b09e6ce2c965e41b337968"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.1"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "57f7cde02d7a53c9d1d28443b9f11ac5fbe7ebc9"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.3"

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
git-tree-sha1 = "ba2d094a88b6b287bd25cfa86f301e7693ffae2f"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.7.4"

[[deps.GridVisualize]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "ElasticArrays", "ExtendableGrids", "GeometryBasics", "GridVisualizeTools", "HypertextLiteral", "LinearAlgebra", "Observables", "OrderedCollections", "PkgVersion", "Printf", "StaticArrays"]
git-tree-sha1 = "8713676cbb58adb61cd3b29ec4f0c618f8b32ceb"
uuid = "5eed8a63-0fb0-45eb-886d-8d5a387d12b8"
version = "0.6.4"

[[deps.GridVisualizeTools]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "StaticArraysCore"]
git-tree-sha1 = "5964fd3e4080af45bfdbdaff75567759fd0367bd"
uuid = "5573ae12-3b76-41d9-b48c-81d0b6e61cc5"
version = "0.2.1"

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
deps = ["ArrayInterfaceCore", "ArrayInterfaceStaticArrays", "ArrayInterfaceStaticArraysCore", "ChainRulesCore", "ForwardDiff", "LinearAlgebra", "MacroTools", "PreallocationTools", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "0a92979c14dfa71adbf892f0cd073e34b7189197"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.13.0"

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
deps = ["ArrayInterface", "ArrayInterfaceOffsetArrays", "ArrayInterfaceStaticArrays", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static"]
git-tree-sha1 = "0ad6f0c51ce004dadc24a28a0dfecfb568e52242"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.13"

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
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearSolve]]
deps = ["ArrayInterfaceCore", "DocStringExtensions", "FastLapackInterface", "GPUArraysCore", "IterativeSolvers", "KLU", "Krylov", "KrylovKit", "LinearAlgebra", "Preferences", "RecursiveFactorization", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "SnoopPrecompile", "SparseArrays", "Sparspak", "SuiteSparse", "UnPack"]
git-tree-sha1 = "ed97c2b4e46d02d4c866d3ccfae039a6c09568b1"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "1.35.0"

    [deps.LinearSolve.extensions]
    LinearSolveHYPRE = "HYPRE"

    [deps.LinearSolve.weakdeps]
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "45b288af6956e67e621c5cbb2d75a261ab58300b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.20"
weakdeps = ["ChainRulesCore", "ChangesOfVariables", "InverseFunctions"]

    [deps.LogExpFunctions.extensions]
    ChainRulesCoreExt = "ChainRulesCore"
    ChangesOfVariablesExt = "ChangesOfVariables"
    InverseFunctionsExt = "InverseFunctions"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "ArrayInterfaceCore", "ArrayInterfaceOffsetArrays", "ArrayInterfaceStaticArrays", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "SIMDTypes", "SLEEFPirates", "SnoopPrecompile", "Static", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "9696a80c21a56b937e3fd89e972f8db5db3186e2"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.150"
weakdeps = ["ChainRulesCore", "ForwardDiff", "SpecialFunctions"]

    [deps.LoopVectorization.extensions]
    ForwardDiffExt = ["ChainRulesCore", "ForwardDiff"]
    SpecialFunctionsExt = "SpecialFunctions"

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

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.MeshIO]]
deps = ["ColorTypes", "FileIO", "GeometryBasics", "Printf"]
git-tree-sha1 = "8be09d84a2d597c7c0c34d7d604c039c9763e48c"
uuid = "7269a6da-0436-5bbc-96c2-40638cbb6118"
version = "0.4.10"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "393fc4d82a73c6fe0e2963dd7c882b09257be537"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.4.6"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "aa532179d4a643d4bd9f328589ca01fa20a0d197"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.1.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

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
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

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
git-tree-sha1 = "151d91d63d8d6c1a5789ecb7de51547e00480f1b"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.4"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f6cf8e7944e50901594838951729a1861e668cb8"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eadad7b14cf046de6eb41f13c9275e5aa2711ab6"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.49"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "e8e0fabcff4df8686c4267503887202a783d498e"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.2"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "240d7170f5ffdb285f9427b92333c3463bf65bf6"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.1"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ForwardDiff"]
git-tree-sha1 = "758f3283aba57c53960c8e1900b4c724bf24ba74"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.8"

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

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "62f417f6ad727987c755549e9cd88c46578da562"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.95.1"

[[deps.PyPlot]]
deps = ["Colors", "LaTeXStrings", "PyCall", "Sockets", "Test", "VersionParsing"]
git-tree-sha1 = "f9d953684d4d21e947cb6d642db18853d43cb027"
uuid = "d330b81b-6aea-500a-939a-2ce795aea3ee"
version = "2.11.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "de191bc385072cc6c7ed3ffdc1caeed3f22c74d4"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.7.0"

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
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceStaticArraysCore", "ChainRulesCore", "DocStringExtensions", "FillArrays", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "ZygoteRules"]
git-tree-sha1 = "f311e004143b4dc7c5492a2e9b9a1d945058fd8c"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.36.0"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "SnoopPrecompile", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "315b2c85818eea6ad1b6b84fd4ecb40cd4146665"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.17"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

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
deps = ["ArrayInterfaceCore", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Preferences", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "9a81b4a706217684f5dbffc22662d93659db96fa"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.82.0"

[[deps.SciMLOperators]]
deps = ["ArrayInterfaceCore", "DocStringExtensions", "Lazy", "LinearAlgebra", "Setfield", "SparseArrays", "StaticArraysCore", "Tricks"]
git-tree-sha1 = "60dc07c77bc831f97945ab1545a5e83252a85342"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.1.19"

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

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimplexGridFactory]]
deps = ["DocStringExtensions", "ElasticArrays", "ExtendableGrids", "FileIO", "GridVisualize", "LinearAlgebra", "MeshIO", "Printf", "Test"]
git-tree-sha1 = "4566d826852b7815d34c7a8829679c6f15f4b2e7"
uuid = "57bfcd06-606e-45d6-baf4-4ba06da0efd5"
version = "0.5.18"

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
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SparseDiffTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceStaticArrays", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "4245283bee733122a9cb4545748d64e0c63337c0"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.30.0"

[[deps.Sparspak]]
deps = ["Libdl", "LinearAlgebra", "Logging", "OffsetArrays", "Printf", "SparseArrays", "Test"]
git-tree-sha1 = "d844b30825ecfd478594d3d500ed8581e1bf03b8"
uuid = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
version = "0.3.8"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "c35b107b61e7f34fa3f124026f2a9be97dea9e1c"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.3"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "6954a456979f23d05085727adb17c4551c19ecd1"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.12"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

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
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "SIMDTypes", "Static", "ThreadingUtilities"]
git-tree-sha1 = "45190b743cdc6f761da1e079bb15ff103a89069c"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.4.6"

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
version = "5.10.1+6"

[[deps.SymbolicIndexingInterface]]
deps = ["DocStringExtensions"]
git-tree-sha1 = "6b764c160547240d868be4e961a5037f47ad7379"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.2.1"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TimerOutputs", "Unityper"]
git-tree-sha1 = "348ad5af9c916b6e1641c74378fac8bb49236688"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "1.0.1"

[[deps.Symbolics]]
deps = ["ArrayInterfaceCore", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "Groebner", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Markdown", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TreeViews"]
git-tree-sha1 = "eb34c1e20b225c1de5adeeff3a085c9e985df532"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "5.0.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

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
version = "1.10.0"

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

[[deps.Triangle_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "fe28e9a4684f6f54e868b9136afb8fd11f1734a7"
uuid = "5639c1d2-226c-5e70-8d55-b3095415a16a"
version = "1.6.2+0"

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "31eedbc0b6d07c08a700e26d31298ac27ef330eb"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.19"

[[deps.Triangulate]]
deps = ["DocStringExtensions", "Libdl", "Printf", "Test", "Triangle_jll"]
git-tree-sha1 = "bbca6ec35426334d615f58859ad40c96d3a4a1f9"
uuid = "f7e6ffb2-c36d-4f8f-a77e-16e897189344"
version = "2.2.0"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "ac00576f90d8a259f2c9d823e91d1de3fd44d348"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.1"

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
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static"]
git-tree-sha1 = "4c59c2df8d2676c4691a39fa70495a6db0c5d290"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.58"

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.VoronoiFVM]]
deps = ["BandedMatrices", "CommonSolve", "DiffResults", "DocStringExtensions", "ExtendableGrids", "ExtendableSparse", "ForwardDiff", "GridVisualize", "InteractiveUtils", "JLD2", "LinearAlgebra", "LinearSolve", "Printf", "Random", "RecursiveArrayTools", "SnoopPrecompile", "SparseArrays", "SparseDiffTools", "StaticArrays", "Statistics", "SuiteSparse", "Symbolics", "Test"]
git-tree-sha1 = "eed07777f5b4ab286e75997ddbfc87443aebd432"
uuid = "82b139dc-5afc-11e9-35da-9b9bdfd336f3"
version = "0.19.3"

[[deps.WriteVTK]]
deps = ["Base64", "CodecZlib", "FillArrays", "LightXML", "TranscodingStreams"]
git-tree-sha1 = "f1b5fce1c4849bc49212c5f471284abf11c57eec"
uuid = "64499a7a-5c06-52f2-abe2-ccb03c286192"
version = "1.17.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

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
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╟─cc49b312-7e21-49aa-b42d-32c184919aa2
# ╟─b61dce01-f84f-4e46-be76-8561bf0dddaa
# ╟─f37e9de2-a681-45d8-9a0a-f53613f19319
# ╠═b912c3ff-42aa-4484-8b84-f011c80382e9
# ╠═8b806c22-3f34-11eb-0e03-4d3fa6461629
# ╠═40795d22-6091-4498-bd07-f716fb07682e
# ╠═d55fdbfc-3f34-11eb-3f35-c3bcdb156288
# ╠═36ef21e0-078c-4287-a0fe-f942cd14dc08
# ╟─d6276724-6f55-426c-aba9-fb89d2c8166e
# ╟─b4b335ef-de2b-4563-8f0f-cde4a04ac115
# ╟─d95afcb5-68f9-4d2e-9111-06fc38fe949a
# ╟─f9886dc9-5b1b-4c27-ac0d-5e65805e45bb
# ╠═add94eb5-7859-474b-87fc-a4e4ab00e691
# ╠═427504f5-6c74-40d5-a4b9-9d2d7bd4f540
# ╠═d908722f-a92d-49de-a3c2-aae33cc5a2e4
# ╠═6f108820-4d06-4fa3-b6ad-8f7454494670
# ╟─fca05141-71cb-4ecc-83b9-f50836dda100
# ╠═a8a69c7c-810a-4ea5-9a72-e0a2eea49782
# ╠═bdc68399-19eb-48b5-ba69-5745a2541a17
# ╠═ba273c41-f622-4f1a-9e1c-48b4143a6d43
# ╠═92dadccc-44e9-49e5-ac94-a224dd8b61aa
# ╠═40dc1b8d-9f11-4ae2-b566-793c3fea99fc
# ╠═ead2425e-0109-4798-b79c-cbc8f76c1414
# ╟─d5946f89-d8d7-4196-90a6-8e34723265c3
# ╟─b16061ae-8a08-4649-af8e-58ddf054996d
# ╟─3d0eebe4-5eaa-4339-b6f6-bb5acd573033
# ╠═d7910827-238a-4b71-ad44-7ce43ff1c0db
# ╠═7d902e19-4ced-4685-ac5b-129ce101c9fe
# ╟─7b51cf5a-fa12-4756-a016-ed00cf90b1ba
# ╠═3b429e0c-c07a-40df-ab0b-1f71d515ece3
# ╠═c770801b-d8f8-4095-98ff-e546f661e141
# ╠═fd59fa72-bc40-41ef-bdac-485b171da03e
# ╟─0875ce0d-a22b-4e08-9b58-57673c986a67
# ╠═a9d85817-df56-48ea-8aee-f700dc6310c8
# ╟─b7a40324-dea5-4fd8-9d7a-1031dceee305
# ╠═905fd6fe-4858-4d85-aeba-7e57439531f5
# ╠═80cfa231-23b1-48a3-a5ad-ee0e122da746
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
