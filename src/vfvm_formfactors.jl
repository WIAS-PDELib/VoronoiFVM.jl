#=
Form factor + edge velocity calculation
=#

################################################
"""
$(SIGNATURES)

Calculate node volume  and voronoi surface contributions for cell.
"""

function cellfactors!(T::Type{Edge1D}, ::Type{Cartesian1D}, coord, cellnodes, icell, nodefac, edgefac)
    en = local_celledgenodes(T)
    K = cellnodes[en[1, 1], icell]
    L = cellnodes[en[2, 1], icell]
    xK = coord[1, K]
    xL = coord[1, L]
    d = abs(xL[1] - xK[1])
    nodefac[1] = d / 2
    nodefac[2] = d / 2
    edgefac[1] = 1 / d
    return nothing
end

function cellfactors!(T::Type{Edge1D}, ::Type{<:Polar1D}, coord, cellnodes, icell, nodefac, edgefac)
    en = local_celledgenodes(T)
    K = cellnodes[en[1, 1], icell]
    L = cellnodes[en[2, 1], icell]
    xK = coord[1, K]
    xL = coord[1, L]
    r0 = xK[1]
    r1 = xL[1]
    if r1 < r0
        r0 = xL[1]
        r1 = xK[1]
    end
    rhalf = 0.5 * (r1 + r0)
    # cpar[1]= π*(r1*r1-r0*r0);         # circular volume
    nodefac[1] = π * (rhalf * rhalf - r0 * r0)   # circular volume between midline and boundary
    nodefac[2] = π * (r1 * r1 - rhalf * rhalf)   # circular volume between midline and boundary
    edgefac[1] = 2.0 * π * rhalf / (r1 - r0)     # circular surface / width
    return nothing
end

function cellfactors!(T::Type{Edge1D}, ::Type{<:Spherical1D}, coord, cellnodes, icell, nodefac, edgefac)
    en = local_celledgenodes(T)
    K = cellnodes[en[1, 1], icell]
    L = cellnodes[en[2, 1], icell]
    xK = coord[1, K]
    xL = coord[1, L]
    r0 = xK[1]
    r1 = xL[1]
    if r1 < r0
        r0 = xL[1]
        r1 = xK[1]
    end
    rhalf = 0.5 * (r1 + r0)
    nodefac[1] = π * (rhalf^3 - r0^3) * 4.0 / 3.0   # sphere volume between midline and boundary
    nodefac[2] = π * (r1^3 - rhalf^3) * 4.0 / 3.0   # sphere volume between midline and boundary
    edgefac[1] = 4.0 * π * rhalf^2 / (r1 - r0)      # circular surface / width
    return nothing
end

function cellfactors!(T::Type{Triangle2D}, ::Type{Cartesian2D}, coord, cellnodes, icell, npar, epar)
    en = local_celledgenodes(T)
    @views n = cellnodes[:, icell]

    # Fill matrix of edge vectors
    V = @SMatrix[
        coord[1, n[en[1, 1]]] - coord[1, n[en[2, 1]]] coord[1, n[en[1, 2]]] - coord[1, n[en[2, 2]]] coord[1, n[en[1, 3]]] - coord[1, n[en[2, 3]]];
        coord[2, n[en[1, 1]]] - coord[2, n[en[2, 1]]] coord[2, n[en[1, 2]]] - coord[2, n[en[2, 2]]] coord[2, n[en[1, 3]]] - coord[2, n[en[2, 3]]]
    ]

    # Compute determinant
    det = V[1, 3] * V[2, 2] - V[1, 2] * V[2, 3]
    vol = abs(0.5 * det)
    ivol = 1.0 / vol

    # squares of edge lengths
    dd = (
        V[1, 1] * V[1, 1] + V[2, 1] * V[2, 1],
        V[1, 2] * V[1, 2] + V[2, 2] * V[2, 2],
        V[1, 3] * V[1, 3] + V[2, 3] * V[2, 3],
    )

    # contributions to \sigma_kl/h_kl
    epar[1] = (dd[2] + dd[3] - dd[1]) * 0.125 * ivol
    epar[2] = (dd[3] + dd[1] - dd[2]) * 0.125 * ivol
    epar[3] = (dd[1] + dd[2] - dd[3]) * 0.125 * ivol

    # contributions to \omega_k

    npar .= 0.0
    for i in 1:3
        npar[en[1, i]] += epar[i] * dd[i] * 0.25
        npar[en[2, i]] += epar[i] * dd[i] * 0.25
    end

    return nothing
end

function cellfactors!(T::Type{Triangle2D}, ::Type{Cylindrical2D}, coord, cellnodes, icell, npar, epar)
    en = local_celledgenodes(T)
    @views n = cellnodes[:, icell]

    # Fill matrix of edge vectors
    V = @SMatrix[
        coord[1, n[en[1, 1]]] - coord[1, n[en[2, 1]]] coord[1, n[en[1, 2]]] - coord[1, n[en[2, 2]]] coord[1, n[en[1, 3]]] - coord[1, n[en[2, 3]]];
        coord[2, n[en[1, 1]]] - coord[2, n[en[2, 1]]] coord[2, n[en[1, 2]]] - coord[2, n[en[2, 2]]] coord[2, n[en[1, 3]]] - coord[2, n[en[2, 3]]]
    ]

    # Compute determinant
    det = V[1, 3] * V[2, 2] - V[1, 2] * V[2, 3]
    vol = abs(0.5 * det)
    ivol = 1.0 / vol

    # squares of edge lengths
    dd = (
        V[1, 1] * V[1, 1] + V[2, 1] * V[2, 1],
        V[1, 2] * V[1, 2] + V[2, 2] * V[2, 2],
        V[1, 3] * V[1, 3] + V[2, 3] * V[2, 3],
    )

    # Edge midpoints
    emid = @SArray[
        0.5 * (coord[1, n[en[1, 1]]] + coord[1, n[en[2, 1]]]) 0.5 * (coord[1, n[en[1, 2]]] + coord[1, n[en[2, 2]]]) 0.5 * (coord[1, n[en[1, 3]]] + coord[1, n[en[2, 3]]]);
        0.5 * (coord[2, n[en[1, 1]]] + coord[2, n[en[2, 1]]]) 0.5 * (coord[2, n[en[1, 2]]] + coord[2, n[en[2, 2]]]) 0.5 * (coord[2, n[en[1, 3]]] + coord[2, n[en[2, 3]]])
    ]

    # Circumcenter; use epar as temp storage
    @views tricircumcenter!(epar, coord[:, n[1]], coord[:, n[2]], coord[:, n[3]])
    rcc = epar[1]

    # contributions to \sigma_kl/h_kl (cartesian)
    epar[1] = (dd[2] + dd[3] - dd[1]) * 0.125 * ivol
    epar[2] = (dd[3] + dd[1] - dd[2]) * 0.125 * ivol
    epar[3] = (dd[1] + dd[2] - dd[3]) * 0.125 * ivol

    # contributions to \omega_k
    # for integration we multiply cartesian areas with 2\pi the average radius
    npar .= 0.0
    for i in 1:3
        @views r1 = coord[1, n[en[1, i]]]
        @views r2 = coord[1, n[en[2, i]]]

        cylfac1 = 2 * π * (r1 + rcc + emid[1, i]) / 3
        cylfac2 = 2 * π * (r2 + rcc + emid[1, i]) / 3

        npar[en[1, i]] += epar[i] * dd[i] * 0.25 * cylfac1
        npar[en[2, i]] += epar[i] * dd[i] * 0.25 * cylfac2
    end

    # for angular integration we multiply interface lengths with 2\pi the average radius of the interface
    # need to do this after angular integrating volume contributions
    for i in 1:3
        rmid = (rcc + emid[1, i]) / 2
        epar[i] *= 2π * rmid
    end

    return nothing
end

function cellfactors!(T::Type{Tetrahedron3D}, ::Type{Cartesian3D}, coord, cellnodes, icell, npar, epar)
    # Transferred from WIAS/pdelib, (c) J. Fuhrmann, H. Langmach, I. Schmelzer
    @views n = cellnodes[:, icell]
    en = local_celledgenodes(T)
    @views pp1 = en[1, :]
    @views pp2 = en[2, :]

    pi1 = (5, 6, 5, 1)
    pi2 = (6, 3, 1, 4)
    pi3 = (4, 2, 3, 2)

    po1 = (2, 1, 2, 6)
    po2 = (1, 4, 6, 3)
    po3 = (3, 5, 4, 5)

    # use epar as intermediate memory
    for i in 1:6
        p1 = cellnodes[pp1[i], icell]
        p2 = cellnodes[pp2[i], icell]
        dx = coord[1, p1] - coord[1, p2]
        dy = coord[2, p1] - coord[2, p2]
        dz = coord[3, p1] - coord[3, p2]
        epar[i] = dx * dx + dy * dy + dz * dz
    end

    # non-allocating tuple
    dd = (epar[1], epar[2], epar[3], epar[4], epar[5], epar[6])
    epar .= 0

    x1 = coord[1, n[2]] - coord[1, n[1]]
    y1 = coord[2, n[2]] - coord[2, n[1]]
    z1 = coord[3, n[2]] - coord[3, n[1]]

    x2 = coord[1, n[3]] - coord[1, n[1]]
    y2 = coord[2, n[3]] - coord[2, n[1]]
    z2 = coord[3, n[3]] - coord[3, n[1]]

    x3 = coord[1, n[4]] - coord[1, n[1]]
    y3 = coord[2, n[4]] - coord[2, n[1]]
    z3 = coord[3, n[4]] - coord[3, n[1]]

    det = (x1 * (y2 * z3 - y3 * z2) + x2 * (y3 * z1 - y1 * z3) + x3 * (y1 * z2 - y2 * z1))

    if det < 0
        det = -det
    end
    vol = det / 6
    vv = 96 * 6 * vol

    for i in 1:4
        npar[i] = 0.0
        i1 = pi1[i]
        i2 = pi2[i]
        i3 = pi3[i]

        h1 = dd[i1] * (dd[i2] + dd[i3] - dd[i1])
        h2 = dd[i2] * (dd[i3] + dd[i1] - dd[i2])
        h3 = dd[i3] * (dd[i1] + dd[i2] - dd[i3])

        df = h1 + h2 + h3
        vf = (h1 * dd[po1[i]] + h2 * dd[po2[i]] + h3 * dd[po3[i]] - 2 * dd[i1] * dd[i2] * dd[i3]) / (vv * df)
        epar[i1] += h1 * vf
        epar[i2] += h2 * vf
        epar[i3] += h3 * vf
    end

    for i in 1:6
        npar[pp1[i]] += epar[i]
        npar[pp2[i]] += epar[i]
        epar[i] = 6 * epar[i] / dd[i]
    end
    return nothing
end

################################################
"""
$(SIGNATURES)

Calculate node volume  contributions for boundary face.
"""

function bfacefactors!(T::Type{Vertex0D}, ::Type{Cartesian1D}, coord, bfacenodes, ibface, nodefac, edgefac)
    nodefac[1] = 1.0
    return nothing
end

function bfacefactors!(T::Type{Vertex0D}, ::Type{<:Polar1D}, coord, bfacenodes, ibface, nodefac, edgefac)
    inode = bfacenodes[1, ibface]
    r = coord[1, inode]
    nodefac[1] = 2 * π * r
    return nothing
end

function bfacefactors!(T::Type{Vertex0D}, ::Type{<:Spherical1D}, coord, bfacenodes, ibface, nodefac, edgefac)
    inode = bfacenodes[1, ibface]
    r = coord[1, inode]
    nodefac[1] = 4 * π * r^2
    return nothing
end

function bfacefactors!(T::Type{Edge1D}, ::Type{Cartesian2D}, coord, bfacenodes, ibface, nodefac, edgefac)
    en = local_celledgenodes(T)
    i1 = bfacenodes[en[1, 1], ibface]
    i2 = bfacenodes[en[2, 1], ibface]
    dx = coord[1, i1] - coord[1, i2]
    dy = coord[2, i1] - coord[2, i2]
    d = sqrt(dx * dx + dy * dy)
    nodefac[1] = d / 2
    nodefac[2] = d / 2
    edgefac[1] = 1 / d
    return nothing
end

function bfacefactors!(T::Type{Edge1D}, ::Type{<:Cylindrical2D}, coord, bfacenodes, ibface, nodefac, edgefac)
    en = local_celledgenodes(T)
    i1 = bfacenodes[en[1, 1], ibface]
    i2 = bfacenodes[en[2, 1], ibface]
    r1 = coord[1, i1]
    r2 = coord[1, i2]
    z1 = coord[2, i1]
    z2 = coord[2, i2]
    dr = r1 - r2
    rmid = (r1 + r2) / 2
    dz = z1 - z2
    l = sqrt(dr * dr + dz * dz)
    nodefac[1] = π * (r1 + rmid) * l / 2
    nodefac[2] = π * (r2 + rmid) * l / 2
    return nothing
end

function bfacefactors!(T::Type{Triangle2D}, ::Type{<:Cartesian3D}, coord, bfacenodes, ibface, npar, epar)
    # Transferred from WIAS/pdelib, (c) J. Fuhrmann, H. Langmach, I. Schmelzer

    en = local_celledgenodes(T)
    @views n = bfacenodes[:, ibface]

    epar .= 0
    for j in 1:3
        d = coord[j, n[en[1, 1]]] - coord[j, n[en[2, 1]]]
        epar[1] += d * d
        d = coord[j, n[en[1, 2]]] - coord[j, n[en[2, 2]]]
        epar[2] += d * d
        d = coord[j, n[en[1, 3]]] - coord[j, n[en[2, 3]]]
        epar[3] += d * d
    end

    dd = (epar[1], epar[2], epar[3])

    # Kanten-Flaechenanteile (ohne Abschneiden); epar -> Hilfsfeld
    epar[1] = (dd[2] + dd[3] - dd[1]) * dd[1]
    epar[2] = (dd[3] + dd[1] - dd[2]) * dd[2]
    epar[3] = (dd[1] + dd[2] - dd[3]) * dd[3]
    vol = sqrt(epar[1] + epar[2] + epar[3]) * 0.25

    d = 1.0 / (8 * vol)

    # Knoten-Flaechenanteile (ohne Abschneiden)
    npar .= 0.0
    for i in 1:3
        npar[en[1, i]] += epar[i] * d * 0.25
        npar[en[2, i]] += epar[i] * d * 0.25
    end

    # Kantengewichte
    epar[1] = epar[1] * d / dd[1]
    epar[2] = epar[2] * d / dd[2]
    epar[3] = epar[3] * d / dd[3]

    return nothing
end

##################################################################

#
# TODO: this should be generalized for more quadrules
#
"""
$(SIGNATURES)

This is an internal function to integrate `velofunc` along the edge ``\\sigma=\\overline{\\mathtt{coordl}\\,\\mathtt{coordr}}``
between the ``x_K`` and ``x_L`` where ``\\mathtt{hnormal}=x_K-x_L`` using [Simpson's Rule](https://en.wikipedia.org/wiki/Simpson%27s_rule).
To be precise, compute for a cartesian coordinate system:
``\\int_{\\sigma} \\mathbf{v} \\cdot \\mathbf{n} \\,\\mathrm{ds} \\lvert x_K - x_L \\rvert / \\lvert\\sigma\\rvert``.
"""
function integrate(::Type{<:Cartesian2D}, coordl, coordr, hnormal, velofunc; kwargs...)
    wl = 1.0 / 6.0
    wm = 2.0 / 3.0
    wr = 1.0 / 6.0
    coordm = (0.5 * (coordl[1] + coordr[1]), 0.5 * (coordl[2] + coordr[2]))
    (vxl, vyl) = velofunc(coordl[1], coordl[2])
    (vxm, vym) = velofunc(coordm[1], coordm[2])
    (vxr, vyr) = velofunc(coordr[1], coordr[2])
    return (wl * vxl + wm * vxm + wr * vxr) * hnormal[1] + (wl * vyl + wm * vym + wr * vyr) * hnormal[2]
end

"""
$(SIGNATURES)

This is an internal function similar to `integrate(::Type{<:Cartesian2D},...)`, but computes instead
``\\int_{\\sigma} r \\, \\mathbf{v} \\cdot \\mathbf{n} \\,\\mathrm{ds} \\lvert x_K - x_L \\rvert / \\left ( \\lvert\\sigma\\rvert r(\\mathrm{mid}(\\sigma)) \\right )``
where ``r(\\mathrm{mid}(\\sigma))`` is the ``r``-coordinate of the mid-point of ``\\sigma``.
"""
function integrate(::Type{<:Cylindrical2D}, coordl, coordr, hnormal, velofunc; kwargs...)
    wl = 1.0 / 6.0
    wm = 2.0 / 3.0
    wr = 1.0 / 6.0
    coordm = (0.5 * (coordl[1] + coordr[1]), 0.5 * (coordl[2] + coordr[2]))

    rl = coordl[1]
    rm = coordm[1]
    rr = coordr[1]

    if abs(coordm[1]) < eps()
        return 0
    else
        (vxl, vyl) = velofunc(coordl[1], coordl[2])
        (vxm, vym) = velofunc(coordm[1], coordm[2])
        (vxr, vyr) = velofunc(coordr[1], coordr[2])

        vl = vxl * hnormal[1] + vyl * hnormal[2]
        vm = vxm * hnormal[1] + vym * hnormal[2]
        vr = vxr * hnormal[1] + vyr * hnormal[2]

        return (wl * vl * rl + wm * vm * rm + wr * vr * rr) / rm
    end
end

"""
$(SIGNATURES)

Project velocity onto grid edges.
That is, we compute the path integrals of the given `velofunc` along the 
Voronoi cell edges as provided by [`integrate`](@ref).
"""
function edgevelocities(grid, velofunc::F; kwargs...) where {F}
    @assert dim_space(grid) < 3

    cn = grid[CellNodes]
    ec = grid[EdgeCells]
    en = grid[EdgeNodes]
    coord = grid[Coordinates]
    coord_system = grid[CoordinateSystem]

    velovec = zeros(Float64, num_edges(grid))
    if dim_space(grid) == 1
        for iedge in 1:num_edges(grid)
            K = en[1, iedge]
            L = en[2, iedge]
            elen = coord[1, L] - coord[1, K]
            vx, vy = velofunc((coord[1, K] + coord[1, L]) / 2, 0)
            velovec[iedge] = -elen * vx
        end
    else
        hnormal = zeros(2)
        p1 = zeros(2)
        p2 = zeros(2)
        for iedge in 1:num_edges(grid)
            K = en[1, iedge]
            L = en[2, iedge]
            @views tricircumcenter!(
                p1,
                coord[:, cn[1, ec[1, iedge]]],
                coord[:, cn[2, ec[1, iedge]]],
                coord[:, cn[3, ec[1, iedge]]]
            )
            if ec[2, iedge] > 0
                @views tricircumcenter!(
                    p2,
                    coord[:, cn[1, ec[2, iedge]]],
                    coord[:, cn[2, ec[2, iedge]]],
                    coord[:, cn[3, ec[2, iedge]]]
                )
            else
                @views @. p2 = 0.5 * (coord[:, K] + coord[:, L])
            end
            @views @. hnormal .= coord[:, K] - coord[:, L]
            velovec[iedge] = integrate(coord_system, p1, p2, hnormal, velofunc; kwargs...)
        end
    end
    return velovec
end

"""
$(SIGNATURES)

Similar to [`edgevelocities`](@ref), but for boundary faces.
"""
function bfacevelocities(grid::ExtendableGrid{Tc, Ti}, velofunc::F; kwargs...) where {Tc, Ti, F}
    @assert dim_space(grid) < 3
    bfacenodes = grid[BFaceNodes]
    coord = grid[Coordinates]
    coord_system = grid[CoordinateSystem]
    bfacecells = grid[BFaceCells]
    bfacenormals::Matrix{Tc} = grid[BFaceNormals]
    bfr = grid[BFaceRegions]
    velovec = zeros(Float64, 2, num_bfaces(grid))
    if dim_space(grid) == 1
        for ibface in 1:num_bfaces(grid)
            vx, vy = velofunc(coord[1, bfacenodes[1, ibface]])
            velovec[ibface] = vx * bfacenormals[1, ibface]
        end
    else
        p1 = zeros(2)
        p2 = zeros(2)
        pm = zeros(2)
        for ibface in 1:num_bfaces(grid)
            @views @. p1 = coord[:, bfacenodes[1, ibface]]
            @views @. p2 = coord[:, bfacenodes[2, ibface]]
            @views @. pm = 0.5 * (p1 + p2)
            @views velovec[1, ibface] = integrate(coord_system, p1, pm, bfacenormals[:, ibface], velofunc; kwargs...)
            @views velovec[2, ibface] = integrate(coord_system, pm, p2, bfacenormals[:, ibface], velofunc; kwargs...)
        end
    end
    return velovec
end

"""
$(SIGNATURES)

Calculate for each Voronoi cell associated with a node of sys.grid
the divergence of the velocity field used to obtain `evelo` and `bfvelo` via 
[`edgevelocities`](@ref) and [`bfacevelocities`](@ref) by means of summing
all `evelos` and `bfvelos` per Voronoi cell.
"""
function calc_divergences(sys, evelo, bfvelo)
    boundary_assem_data = sys.boundary_assembly_data
    grid = sys.grid
    div4nodes = zeros(Float64, num_nodes(grid))

    edge = Edge(sys)
    for item in edgebatch(sys.assembly_data)
        for iedge in edgerange(sys.assembly_data, item)
            _fill!(edge, sys.assembly_data, iedge, item)
            node1 = edge.node[1]
            node2 = edge.node[2]
            div4nodes[node1] -= evelo[edge.index] * edge.fac
            div4nodes[node2] += evelo[edge.index] * edge.fac
        end
    end

    # assem_data=sys.assembly_data
    # edgefactors = assem_data.edgefactors
    # cellnodes = grid[CellNodes]
    # celledges = grid[CellEdges]
    # celledgesigns = grid[CellFaceSigns]
    # en=local_celledgenodes(Triangle2D)
    # for icell in 1:num_cells(grid)
    #     for localedge in 1 : 3
    #         node1 = cellnodes[en[1,localedge], icell]
    #         node2 = cellnodes[en[2,localedge], icell]
    #         sign = celledgesigns[localedge, icell]
    #         div4nodes[node1] -= sign*evelo[celledges[localedge,icell]] * edgefactors[localedge,icell]
    #         div4nodes[node2] += sign*evelo[celledges[localedge,icell]] * edgefactors[localedge,icell]
    #     end
    # end

    bfacenodes = grid[BFaceNodes]
    boundarynodefactors = boundary_assem_data.nodefactors
    for ibface in 1:num_bfaces(grid)
        node1 = bfacenodes[1, ibface]
        node2 = bfacenodes[2, ibface]
        div4nodes[node1] += boundarynodefactors[1, ibface] * bfvelo[1, ibface]
        div4nodes[node2] += boundarynodefactors[2, ibface] * bfvelo[2, ibface]
    end

    return div4nodes
end

function calc_divergences(sys, velofunc)
    grid = sys.grid
    evelo = edgevelocities(grid, velofunc)
    bfvelo = bfacevelocities(grid, velofunc)

    return calc_divergences(sys, evelo, bfvelo)
end

"""
$(SIGNATURES)

Node form factors per boundary face
"""
bfacenodefactors(sys) = sys.boundary_assembly_data.nodefactors
