"""
    time(edge_or_node)

Return actual simulation time stored in node or edge
"""
time(item::AbstractGeometryItem) = item.time

"""
    embedparam(edge_or_node)

Return embedding parameter stored in node or edge
"""
embedparam(item::AbstractGeometryItem) = item.embedparam

"""
   region(edge_or_node)

Return region number node or edge is belonging to
"""
region(item::AbstractGeometryItem) = item.region

"""
   $(TYPEDEF)

Abstract type for nodes. 

`node[idim]` gives the the corresponding coordinate.
"""
abstract type AbstractNode{Tc <: Number, Tp <: Number, Ti <: Integer} <: AbstractGeometryItem{Tc, Tp, Ti} end
Base.size(node::AbstractNode) = (size(node.coord)[1],)
Base.getindex(node::AbstractNode, idim) = @inbounds node.coord[idim, node.index]

"""
    $(TYPEDEF)

Abstract type for data on nodes.
`u[ispec]` accesses value of species at this node.
"""
abstract type AbstractNodeData{Tv <: Number} <: AbstractVector{Tv} end
Base.size(u::AbstractNodeData) = (u.nspec, 1)
Base.getindex(u::AbstractNodeData, i) = @inbounds u.val[i]
Base.setindex!(f::AbstractNodeData, v, i) = @inbounds f.val[i] = v

struct DParameters{Tv <: Number} <: AbstractVector{Tv}
    val::Vector{Tv}
    offset::Int32
end

Base.size(p::DParameters) = (length(p.val) - p.offset, 1)
Base.getindex(p::DParameters, i) = @inbounds p.val[p.offset + i]

"""
    parameters(edge_or_node)

Return abstract vector of parameters passed via vector of unknowns. 
This allows differentiation with respect to these parameters.
"""
function parameters(u::AbstractNodeData)
    return DParameters(u.val, u.nspec)
end

"""
   $(TYPEDEF)

Abstract type for edges 

`edge[idim,inode]` gives coordinate of node.
"""
abstract type AbstractEdge{Tv <: Number, Tp <: Number, Ti <: Integer} <: AbstractGeometryItem{Tv, Tp, Ti} end
Base.size(edge::AbstractEdge) = (size(edge.coord)[1], 2)
Base.getindex(edge::AbstractEdge, idim, inode) = @inbounds edge.coord[idim, edge.node[inode]]

"""
    $(TYPEDEF)

Abstract type for data on edges.
`u[ispec,inode]` accesses value of species at corresponding node.
"""
abstract type AbstractEdgeData{Tv <: Number} <: AbstractMatrix{Tv} end
Base.size(u::AbstractEdgeData) = (u.n1, 2)
Base.getindex(u::AbstractEdgeData, i, j) = @inbounds u.val[(j - 1) * u.n1 + i]

function parameters(u::AbstractEdgeData)
    return DParameters(u.val, u.n1 + u.n1)
end

##################################################################
"""
$(TYPEDEF)

Structure holding local node information.

$(TYPEDFIELDS)
"""
mutable struct Node{Tc, Tp, Ti} <: AbstractNode{Tc, Tp, Ti}
    """
    Index in grid

    """
    index::Ti
    """
    Inner region number
    """
    region::Ti

    """
    Number of species defined in node
    """
    nspec::Ti

    """
    Number of discretization cell the node is invoked from
    """
    icell::Ti

    """
    Grid coordinates
    """
    coord::Matrix{Tc}

    """
    Grid cell nodes
    """
    cellnodes::Array{Ti, 2}

    """
    Grid cell regions
    """
    cellregions::Vector{Ti}

    """
    System time
    """
    time::Float64

    """
    Current value of embedding parameter
    """
    embedparam::Float64

    """
    parameters (deprecated)
    """
    params::Vector{Tp}

    """
    Form factor
    """
    fac::Float64

    """
    Local loop index
    """
    _idx::Ti

    function Node{Tc, Tp, Ti}(sys::AbstractSystem{Tv, Tc, Ti, Tm}, time, embedparam, params::Vector{Tp}) where {Tv, Tc, Tp, Ti, Tm}
        return new(
            zero(Ti), 0,
            num_species(sys), 0,
            coordinates(sys.grid),
            sys.grid[CellNodes],
            sys.grid[CellRegions],
            time, embedparam, params, 0.0, 0
        )
    end
end

function Node(sys::AbstractSystem{Tv, Tc, Ti, Tm}, time, embedparam, params::Vector{Tp}) where {Tv, Tc, Tp, Ti, Tm}
    return Node{Tc, Tp, Ti}(sys, time, embedparam, params)
end

Node(sys) = Node(sys, 0, 0, zeros(0))

"""
    $(TYPEDEF)

Unknown data on node. 
"""
struct NodeUnknowns{Tv, Tc, Tp, Ti} <: AbstractNodeData{Tv}
    val::Vector{Tv}
    nspec::Ti
    geom::Node{Tc, Tp, Ti}
end

@inline function unknowns(node::Node{Tc, Tp, Ti}, u::AbstractVector{Tv}) where {Tv, Tc, Tp, Ti}
    return NodeUnknowns{Tv, Tc, Tp, Ti}(u, node.nspec, node)
end

"""
    $(TYPEDEF)

RHS data on node. 
"""
struct NodeRHS{Tv, Tc, Tp, Ti} <: AbstractNodeData{Tv}
    val::Vector{Tv}
    nspec::Ti
    geom::Node{Tc, Tp, Ti}
end

@inline rhs(node::Node{Tc, Tp, Ti}, f::AbstractVector{Tv}) where {Tv, Tc, Tp, Ti} = NodeRHS{Tv, Tc, Tp, Ti}(f, node.nspec, node)

##################################################################
"""
$(TYPEDEF)

Structure holding local boundary  node information.

$(TYPEDFIELDS)
"""
mutable struct BNode{Tv, Tc, Tp, Ti} <: AbstractNode{Tc, Tp, Ti}
    """
    Index in grid
    """
    index::Ti

    """
    BFace number it is called from
    """
    ibface::Ti

    """
    local node number
    """
    ibnode::Ti

    """
    Boundary region number
    """
    region::Ti

    cellregions::Vector{Ti}

    """
    Number of species defined in node
    """
    nspec::Ti

    """
    Grid coordinates
    """
    coord::Matrix{Tc}

    bfacenodes::Array{Ti, 2}

    bfaceregions::Vector{Ti}

    allcellregions::Vector{Ti}

    bfacecells::ExtendableGrids.Adjacency{Ti}

    Dirichlet::Tv

    """
    System time
    """
    time::Float64

    """
    Current value of embedding parameter
    """
    embedparam::Float64

    """
    Parameters (deprecated)
    """
    params::Vector{Tp}

    dirichlet_value::Vector{Tv}

    fac::Float64

    function BNode{Tv, Tc, Tp, Ti}(
            sys::AbstractSystem{Tv, Tc, Ti, Tm}, time, embedparam,
            params::Vector{Tp}
        ) where {Tv, Tc, Tp, Ti, Tm}
        return new(
            0, 0, 0, 0, zeros(Ti, 2),
            num_species(sys),
            coordinates(sys.grid),
            sys.grid[BFaceNodes],
            sys.grid[BFaceRegions],
            sys.grid[CellRegions],
            sys.grid[BFaceCells],
            Dirichlet(Tv), time, embedparam, params,
            zeros(Tv, num_species(sys)), 0.0
        )
    end
end
function BNode(sys::AbstractSystem{Tv, Tc, Ti, Tm}, time, embedparam, params::Vector{Tp}) where {Tv, Tc, Tp, Ti, Tm}
    return BNode{Tv, Tc, Tp, Ti}(sys, time, embedparam, params)
end
BNode(sys) = BNode(sys, 0, 0, zeros(0))

struct BNodeUnknowns{Tval, Tv, Tc, Tp, Ti} <: AbstractNodeData{Tv}
    val::Vector{Tval}
    nspec::Ti
    geom::BNode{Tv, Tc, Tp, Ti}
end

@inline function unknowns(bnode::BNode{Tv, Tc, Tp, Ti}, u::AbstractVector{Tval}) where {Tval, Tv, Tc, Tp, Ti}
    return BNodeUnknowns{Tval, Tv, Tc, Tp, Ti}(u, bnode.nspec, bnode)
end

struct BNodeRHS{Tval, Tv, Tc, Tp, Ti} <: AbstractNodeData{Tv}
    val::Vector{Tval}
    nspec::Ti
    geom::BNode{Tv, Tc, Tp, Ti}
end

@inline function rhs(bnode::BNode{Tv, Tc, Tp, Ti}, f::AbstractVector{Tval}) where {Tval, Tv, Tc, Tp, Ti}
    return BNodeRHS{Tval, Tv, Tc, Tp, Ti}(f, bnode.nspec, bnode)
end

##################################################################
"""
$(TYPEDEF)

Structure holding local edge information.

$(TYPEDFIELDS)
"""
mutable struct Edge{Tc, Tp, Ti} <: AbstractEdge{Tc, Tp, Ti}
    """
    Index in grid
    """
    index::Ti

    """
    Index 
    """
    node::Vector{Ti}

    """
    Inner region number corresponding to edge
    """
    region::Ti

    """
    Number of species defined in edge
    """
    nspec::Ti

    """
    Number of discretization cell the edge is invoked from
    """
    icell::Ti

    """
    Grid coordinates
    """
    coord::Matrix{Tc}

    cellx::Array{Ti, 2}
    edgenodes::Array{Ti, 2}
    cellregions::Vector{Ti}
    has_celledges::Bool

    """
    System time
    """
    time::Float64

    """
    Current value of embedding parameter
    """
    embedparam::Float64

    """
    Parameters (deprecated)
    """
    params::Vector{Tp}

    """
    Form factor
    """
    fac::Float64

    """
    Local loop index
    """
    _idx::Ti

    outflownoderegions::Union{Nothing, SparseMatrixCSC{Bool, Int}}

    """
    Outflow node
    """
    outflownode::Int

    Edge{Tc, Tp, Ti}(::Nothing) where {Tc, Tp, Ti} = new()
end

Edge(sys) = Edge(sys, 0, 0, zeros(0))

function Edge(sys::AbstractSystem{Tv, Tc, Ti, Tm}, time, embedparam, params::Vector{Tp}) where {Tv, Tc, Tp, Ti, Tm}
    edge = Edge{Tc, Tp, Ti}(nothing)

    edge.index = 0
    edge.node = [0, 0]
    edge.region = 0
    edge.nspec = num_species(sys)
    edge.icell = 0
    edge.coord = coordinates(sys.grid)
    geom = sys.grid[CellGeometries][1]
    if haskey(sys.grid, CellEdges)
        edge.cellx = sys.grid[CellEdges]
        edge.edgenodes = sys.grid[EdgeNodes]
        edge.has_celledges = true
    else
        edge.cellx = sys.grid[CellNodes]
        edge.edgenodes = local_celledgenodes(geom)
        edge.has_celledges = false
    end
    edge.cellregions = sys.grid[CellRegions]
    edge.time = time
    edge.embedparam = embedparam
    edge.params = params
    edge.fac = 0
    edge.outflownode = 0
    edge._idx = 0
    edge.outflownoderegions = sys.outflownoderegions
    return edge
end

struct EdgeUnknowns{Tv, Tc, Tp, Ti} <: AbstractEdgeData{Tv}
    val::Vector{Tv}
    n1::Ti
    geom::Edge{Tc, Tp, Ti}
end

@inline function unknowns(edge::Edge{Tc, Tp, Ti}, u::AbstractVector{Tv}) where {Tv, Tc, Tp, Ti}
    return EdgeUnknowns{Tv, Tc, Tp, Ti}(u, edge.nspec, edge)
end

struct EdgeRHS{Tv, Tc, Tp, Ti} <: AbstractNodeData{Tv}
    val::Vector{Tv}
    nspec::Ti
    geom::Edge{Tc, Tp, Ti}
end

@inline rhs(edge::Edge{Tc, Tp, Ti}, f::AbstractVector{Tv}) where {Tv, Tc, Tp, Ti} = EdgeRHS{Tv, Tc, Tp, Ti}(f, edge.nspec, edge)

"""
    hasoutflownode(edge)

Check if one node of the edge is situated on a boundary region listed in `outflowboundaries`, see
[`struct Physics`].
"""
hasoutflownode(edge) = isoutflownode(edge, 1) || isoutflownode(edge, 2)

"""
    isoutflownode(edge,inode)

Check if inode (1 or 2) is an outflow node.
"""
isoutflownode(edge, inode) = length(nzrange(edge.outflownoderegions, edge.node[inode])) > 0

"""
    isoutflownode(edge,inode,irefgion)

Check if inode (1 or 2) is an outflow node on boundary region `iregion`.
"""
isoutflownode(edge, inode, iregion) = edge.outflownoderegions[iregion, edge.node[inode]]

"""
    outflownode(edge)

Return outflow node of edge (1 or 2).
"""
outflownode(edge) = edge.outflownode

"""
    outflownode!(edge)
Set `edge.outflownode` entry.
"""
function outflownode!(edge)
    isoutflownode(edge, 1) ? edge.outflownode = 1 : true
    return isoutflownode(edge, 2) ? edge.outflownode = 2 : true
end

##################################################################
"""
$(TYPEDEF)

Structure holding local edge information.

$(TYPEDFIELDS)
"""
mutable struct BEdge{Tc, Tp, Ti} <: AbstractEdge{Tc, Tp, Ti}
    """
    Index in grid
    """
    index::Ti

    """
    Index 
    """
    node::Vector{Ti}

    """
    Inner region number corresponding to edge
    """
    region::Ti

    """
    Number of species defined in edge
    """
    nspec::Ti

    """
    Number of discretization cell the edge is invoked from
    """
    icell::Ti

    """
    Grid coordinates
    """
    coord::Matrix{Tc}

    bedgenodes::Array{Ti, 2}
    bfaceedges::Array{Ti, 2}
    bfaceregions::Vector{Ti}

    """
    System time
    """
    time::Float64

    """
    Current value of embedding parameter
    """
    embedparam::Float64

    """
    Parameters (deprecated)
    """
    params::Vector{Tp}

    fac::Float64
    BEdge{Tc, Tp, Ti}(::Nothing) where {Tc, Tp, Ti} = new()
end

BEdge(sys) = BEdge(sys, 0, 0, zeros(0))

function BEdge(sys::AbstractSystem{Tv, Tc, Ti, Tm}, time, embedparam, params::Vector{Tp}) where {Tv, Tc, Tp, Ti, Tm}
    bedge = BEdge{Tc, Tp, Ti}(nothing)

    bedge.index = 0
    bedge.node = [0, 0]
    bedge.region = 0
    bedge.nspec = num_species(sys)
    bedge.icell = 0
    bedge.coord = coordinates(sys.grid)

    bedge.bfaceedges = sys.grid[BFaceEdges] # !!! another bug in ExtendableGrids
    bedge.bedgenodes = sys.grid[BEdgeNodes]
    bedge.bfaceregions = sys.grid[BFaceRegions]
    bedge.time = time
    bedge.embedparam = embedparam
    bedge.params = params
    bedge.fac = 0.0
    return bedge
end

struct BEdgeUnknowns{Tv, Tc, Tp, Ti} <: AbstractEdgeData{Tv}
    val::Vector{Tv}
    n1::Ti
    geom::BEdge{Tc, Tp, Ti}
end

@inline function unknowns(edge::BEdge{Tc, Tp, Ti}, u::AbstractVector{Tv}) where {Tv, Tc, Tp, Ti}
    return BEdgeUnknowns{Tv, Tc, Tp, Ti}(u, edge.nspec, edge)
end

struct BEdgeRHS{Tv, Tc, Tp, Ti} <: AbstractNodeData{Tv}
    val::Vector{Tv}
    nspec::Ti
    geom::BEdge{Tc, Tp, Ti}
end

@inline rhs(edge::BEdge{Tc, Tp, Ti}, f::AbstractVector{Tv}) where {Tv, Tc, Tp, Ti} = BEdgeRHS{Tv, Tc, Tp, Ti}(f, edge.nspec, edge)

##################################################################
"""
$(TYPEDSIGNATURES)

Return number of species for edge
"""
@inline num_species(edge::AbstractEdge) = edge.nspec

##################################################################
"""
$(TYPEDSIGNATURES)
   
Calculate the length of an edge. 
"""
function meas(edge::AbstractEdge)
    l = 0.0
    for i in 1:size(edge.coord)[1]
        d = edge.coord[i, edge.node[1]] - edge.coord[i, edge.node[2]]
        l = l + d * d
    end
    return sqrt(l)
end

"""
    edgelength(edge)

Return length of edge
"""
edgelength(edge::AbstractEdge) = meas(edge)

"""
    project(edge, vector)

Project d-vector onto d-dimensional vector, i.e. calculate the dot product
of `vector` with the difference of the edge end coordinates.
"""
function project(edge::Edge, vec)
    vh = zero(eltype(vec))
    for i in 1:size(edge.coord)[1]
        vh += (edge.coord[i, edge.node[2]] - edge.coord[i, edge.node[1]]) * vec[i]
    end
    return vh
end

###############################################################
# Deprecation warnings here ?
# Wrapper struct for viewing unknowns passed to callback functions
struct VectorUnknowns{Tv} <: AbstractVector{Tv}
    val::Vector{Tv}
    n::Int64
    offset::Int64
end

unknowns(edge::AbstractEdge, u::AbstractVector{Tv}, i) where {Tv} = VectorUnknowns{Tv}(u, edge.nspec, (i - 1) * (edge.nspec))
Base.getindex(u::VectorUnknowns, i) = @inbounds u.val[u.offset + i]
Base.size(u::VectorUnknowns) = (u.n,)

# For backward compatibility
unknowns(edge, u::AbstractEdgeData) = u
# For backward compatibility
unknowns(edge::Edge, u::AbstractEdgeData{Tv}, i) where {Tv} = VectorUnknowns{Tv}(u.val, edge.nspec, (i - 1) * (edge.nspec))

"""
$(TYPEDSIGNATURES)

Solution view on first edge node
"""
viewK(edge::AbstractEdge, u) = unknowns(edge, u, 1)

"""
$(TYPEDSIGNATURES)

Solution view on second edge node
"""
viewL(edge::AbstractEdge, u) = unknowns(edge, u, 2)
