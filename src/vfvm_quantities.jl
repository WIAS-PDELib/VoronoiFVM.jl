##########################################################
"""
$(TYPEDEF)

Abstract supertype of quantities
"""
abstract type AbstractQuantity{Ti <: Integer} end

##########################################################
"""
$(TYPEDEF)

A continuous quantity is represented by exactly one species number

$(TYPEDFIELDS)
"""
struct ContinuousQuantity{Ti} <: AbstractQuantity{Ti}
    """
    Species number representing the quantity
    """
    ispec::Ti

    """
    Quantity identifier allowing to use the quantity as index
    in parameter fields
    """
    id::Ti
end

"""
     ContinuousQuantity(system,regions; ispec=0, id=0)


Add continuous quantity to the regions listed in `regions`.

Unless specified in `ispec`, the species number is generated
automatically.

Unless specified by `id`, the quantity ID is generated automatically.
"""
function ContinuousQuantity(sys::AbstractSystem{Tv, Tc, Ti, Tm}, regions; ispec = 0, id = 0) where {Tv, Tc, Ti, Tm}
    if ispec == 0
        nspec = num_species(sys)
        nspec = nspec + 1
    else
        nspec = ispec
    end
    sys.num_quantities += 1
    if id == 0
        id = sys.num_quantities
    end
    enable_species!(sys, nspec, regions)
    return ContinuousQuantity{Ti}(nspec, id)
end

##########################################################
"""
$(TYPEDEF)

An interface quantity is represented by exactly one species number

$(TYPEDFIELDS)
"""
struct InterfaceQuantity{Ti} <: AbstractQuantity{Ti}
    """
    Species number representing the quantity
    """
    ispec::Ti

    """
    boundary regions, where interface quantity is defined
    """
    bregspec::Vector{Ti}

    """
    Quantity identifier allowing to use the quantity as index
    in parameter fields
    """
    id::Ti
end

"""
     InterfaceQuantity(system,regions; ispec=0, id=0)


Add interface quantity to the boundary regions given in `breg`.

Unless specified in `ispec`, the species number is generated
automatically.

Unless specified by `id`, the quantity ID is generated automatically.
"""
function InterfaceQuantity(sys::AbstractSystem{Tv, Tc, Ti, Tm}, bregions::AbstractVector; ispec = 0, id = 0) where {Tv, Tc, Ti, Tm}
    if ispec == 0
        nspec = num_species(sys)
        nspec = nspec + 1
    else
        nspec = ispec
    end
    enable_boundary_species!(sys, nspec, bregions)
    sys.num_quantities += 1
    if id == 0
        id = sys.num_quantities
    end
    return InterfaceQuantity{Ti}(nspec, bregions, id)
end

###########################################################
"""
$(TYPEDEF)

A discontinuous quantity is represented by different species
in neighboring regions.

$(TYPEDFIELDS)
"""
struct DiscontinuousQuantity{Ti} <: AbstractQuantity{Ti}
    """
    Species numbers representing the quantity in each region
    """
    regionspec::Vector{Ti}

    """
    Quantity identifier allowing to use the quantity as index
    in parameter fields
    """
    id::Ti
end

"""
     DiscontinuousQuantity(system,regions; regionspec=nothing, id=0)


Add discontinuous quantity to the regions listed in `regions`.

Unless specified in `regionspec`, the species numbers for each region
are generated automatically.

Unless specified by `id`, the quantity ID is generated automatically.
"""
function DiscontinuousQuantity(
        sys::AbstractSystem{Tv, Tc, Ti, Tm}, regions::AbstractVector; regionspec = nothing,
        id = 0
    ) where {Tv, Tc, Ti, Tm}
    rspec = zeros(Ti, num_cellregions(sys.grid))
    if regionspec == nothing
        nspec = num_species(sys)
        for ireg in regions
            nspec = nspec + 1
            enable_species!(sys, nspec, [ireg])
            rspec[ireg] = nspec
        end
    else
        for ireg in regions
            enable_species!(sys, regionspec[ireg], [ireg])
            rspec[ireg] = regionspec[ireg]
        end
    end
    sys.num_quantities += 1
    if id == 0
        id = sys.num_quantities
    end
    quantity = DiscontinuousQuantity{Ti}(rspec, id)
    return quantity
end

"""
    $(SIGNATURES)

Number of quantities defined for system
"""
num_quantities(system::AbstractSystem) = system.num_quantities

"""
    subgrids(quantity, system)

Return a vector of subgrids containing a subgrid for each
region where discontinuous quantity is defined.
"""
function subgrids(quantity::DiscontinuousQuantity, sys)
    grid = sys.grid
    subgrids = Vector{typeof(sys.grid)}(undef, 0)
    for ireg in 1:num_cellregions(grid)
        ispec = quantity.regionspec[ireg]
        if ispec > 0
            push!(subgrids, subgrid(grid, [ireg]))
        end
    end
    return subgrids
end

"""
    subgrids(quantity, system)

Return the subgrid where interface quantity is defined.
"""
function subgrids(quantity::InterfaceQuantity, sys)
    grid = sys.grid
    return bgrid = subgrid(grid, quantity.bregspec; boundary = true)
end

"""
    views(quantity, subgrids,system)

Return a vector of solutions containing the solutions with respect tp
each region where discontinuous quantity is defined.
"""
function views(U, quantity::DiscontinuousQuantity, subgrids, sys)
    grid = sys.grid
    projections = Vector[]
    j = 1
    for ireg in 1:num_cellregions(grid)
        ispec = quantity.regionspec[ireg]
        if ispec > 0
            push!(projections, view(U[ispec, :], subgrids[j]))
            j = j + 1
        end
    end
    return projections
end

function views(U, quantity::ContinuousQuantity, subgrids, sys)
    grid = sys.grid
    projections = Vector[]
    j = 1
    for ireg in 1:num_cellregions(grid)
        push!(projections, view(U[quantity.ispec, :], subgrids[j]))
        j = j + 1
    end
    return projections
end

function views(U, quantity::InterfaceQuantity, bgrid, sys)
    return view(U[quantity.ispec, :], bgrid)
end

# just return the first which comes into mind.
# we need to ensure homogeneity of bregion-region structure.
# if that is true, this works.
function get_cellregion(sys, ibc)
    grid = sys.grid
    bfregions = grid[BFaceRegions]
    cregions = grid[CellRegions]
    for ibface in 1:num_bfaces(sys.grid)
        if bfregions[ibface] == ibc
            bfcells = grid[BFaceCells]
            return cregions[bfcells[ibface, 1]]
        end
    end
    return 0
end

"""
    boundary_dirichlet(system, quantity, ibc, value)

Set Dirichlet boundary `value` for `quantity` at boundary `ibc`.
"""
function boundary_dirichlet!(sys::AbstractSystem, q::DiscontinuousQuantity, ibc, v)
    return boundary_dirichlet!(
        sys,
        q.regionspec[get_cellregion(sys, ibc)],
        ibc,
        v
    )
end

function boundary_dirichlet!(sys::AbstractSystem, q::ContinuousQuantity, ibc, v)
    return boundary_dirichlet!(
        sys,
        q.ispec,
        ibc,
        v
    )
end

function boundary_neumann!(sys::AbstractSystem, q::DiscontinuousQuantity, ibc, v)
    return boundary_neumann!(
        sys,
        q.regionspec[get_cellregion(sys, ibc)],
        ibc,
        v
    )
end

function boundary_neumann!(sys::AbstractSystem, q::ContinuousQuantity, ibc, v)
    return boundary_neumann!(
        sys,
        q.ispec,
        ibc,
        v
    )
end

function boundary_robin!(sys::AbstractSystem, q::DiscontinuousQuantity, ibc, α, v)
    return boundary_robin!(
        sys,
        q.regionspec[get_cellregion(sys, ibc)],
        ibc,
        α,
        v
    )
end

function boundary_robin!(sys::AbstractSystem, q::ContinuousQuantity, ibc, α, v)
    return boundary_robin!(
        sys,
        q.ispec,
        ibc,
        α,
        v
    )
end

"""
    node[quantity]
    edge[quantity]
Return species number on [`AbstractNode`](@ref) or [`AbstractEdge`](@ref)
"""
Base.getindex(q::ContinuousQuantity, node::AbstractNode) = q.ispec
Base.getindex(q::InterfaceQuantity, node::AbstractNode) = q.ispec
Base.getindex(q::AbstractQuantity, edge::AbstractEdge) = q.ispec

Base.getindex(q::DiscontinuousQuantity, edge::Edge) = @inbounds q.regionspec[edge.region]
Base.getindex(q::DiscontinuousQuantity, edge::BEdge) = nothing
Base.getindex(q::DiscontinuousQuantity, node::Node) = @inbounds q.regionspec[node.region]

"""
    bnode[quantity]
Return species number of discontinuous quantity region `ireg`  adjacent
to  [`BNode`](@ref) for outer boundary nodes.
"""
Base.getindex(q::DiscontinuousQuantity, bnode::BNode) = @inbounds q.regionspec[bnode.cellregions[1]]

"""
    bnode[quantity,ireg]
Return species number of discontinuous quantity region `ireg`  adjacent
to  [`BNode`](@ref).
"""
Base.getindex(q::DiscontinuousQuantity, bnode::BNode, j) = @inbounds q.regionspec[bnode.cellregions[j]]

"""
    u[q,j]
Return value of quantity in unknowns on edge in flux callbacks.
"""
Base.getindex(u::AbstractEdgeData, q::AbstractQuantity, j) = @inbounds u[q[u.geom], j]

"""
    u[q]
Return value of quantity in unknowns on node in  node callbacks.
"""
Base.getindex(u::AbstractNodeData, q::AbstractQuantity) = @inbounds u[q[u.geom]]

"""
    f[q]=value
Set rhs value for quantity in callbacks
"""
Base.setindex!(f::AbstractNodeData, v, q::AbstractQuantity) = @inbounds f[q[f.geom]] = v

"""
    u[q,ireg]
Return value of discontinuous quantity in unknowns adjacent to unknowns on boundary node.
"""
Base.getindex(u::BNodeUnknowns, q::DiscontinuousQuantity, j) = @inbounds u[q[u.geom, j]]
Base.getindex(f::BNodeRHS, q::DiscontinuousQuantity, j) = @inbounds f[q[f.geom, j]]

"""
    f[q,ireg]=v
Set rhs value for discontinuous quantity in adjacent regions of  boundary node.
"""
Base.setindex!(f::BNodeRHS, v, q::DiscontinuousQuantity, j) = @inbounds f[q[f.geom, j]] = v

"""
    M[q,i]

Access columns  `M` using id of quantity `q`.
This is meant for nspecies x  nregions matrices e.g. defining parameters.
"""
Base.getindex(m::AbstractMatrix, q::AbstractQuantity, j) = m[q.id, j]

"""
    M[q,i]

Set element of `M` using id of quantity `q`.
This is meant for nspecies x  nregions matrices  e.g. defining parameters.
"""
Base.setindex!(m::AbstractMatrix, v, q::AbstractQuantity, j) = m[q.id, j] = v

"""
    A[q]

Access columns  of vectors `A` using id of quantity `q`.
This is meant for vectors indexed by species.
"""
Base.getindex(A::AbstractArray, q::AbstractQuantity) = A[q.id]

"""
    A[q]

Set element of `A` using id of quantity `q`
This is meant for vectors indexed by species.
"""
Base.setindex!(A::AbstractArray, v, q::AbstractQuantity) = A[q.id] = v

Base.getindex(I::SolutionIntegral, cspec::ContinuousQuantity, ireg) = I[cspec.id, ireg]
Base.getindex(I::SolutionIntegral, dspec::DiscontinuousQuantity, ireg) = I[dspec.regionspec[ireg], ireg]
function Base.getindex(I::SolutionIntegral, dspec::DiscontinuousQuantity, ::Colon)
    return [I[dspec.regionspec[ireg], ireg] for ireg in 1:size(I, 2)]
end
