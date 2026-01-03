##########################################################
"""
$(TYPEDEF)

Abstract type for physics.
"""
abstract type AbstractPhysics end

##########################################################
"""
$(TYPEDEF)

Abstract type for user data.

It is possible but not necessary to make user data a subtype of AbstractData
and get a  prettyprinting show method.
"""
abstract type AbstractData{Tv} end

function _showstruct(io::IO, this::AbstractData)
    myround(x; kwargs...) = round(Float64(value(x)); kwargs...)
    myround(s::Symbol; kwargs...) = s
    myround(i::Int; kwargs...) = i
    myround(b::Bool; kwargs...) = b
    println(io, typeof(this))
    for name in fieldnames(typeof(this))
        println(io, "$(lpad(name, 20)) = $(myround.(getfield(this, name), sigdigits = 5))")
    end
    return
end

"""
    copy!(to::AbstractData, from::AbstractData)

Copy [`AbstractData`](@ref).
"""
function Base.copy!(vdata::AbstractData{Tv}, udata::AbstractData{Tu}) where {Tv, Tu}
    vval(x::Any) = x
    vval(x::Tu) = Tv(x)
    for name in fieldnames(typeof(udata))
        setproperty!(vdata, name, vval(getproperty(udata, name)))
    end
    return vdata
end

"""
    $(TYPEDSIGNATURES)

Pretty print [`AbstractData`](@ref)
"""
Base.show(io::IO, ::MIME"text/plain", this::AbstractData) = _showstruct(io, this)

@doc """
    value(x)

Return the value of a dual number (for debugging in callback functions).
Re-exported from ForwardDiff.
""" value

#
# Dummy callbacks
#
function nofunc(args...) end


function default_storage(f, u, node, data)
    return f .= u
end


##########################################################
"""
````
struct Physics
````

Physics data record with the following fields:

$(TYPEDFIELDS)
"""
mutable struct Physics{
        Flux <: Function,
        Reaction <: Function,
        EdgeReaction <: Function,
        Storage <: Function,
        Source <: Function,
        BFlux <: Function,
        BReaction <: Function,
        BSource <: Function,
        BStorage <: Function,
        BOutflow <: Function,
        GenericOperator <: Function,
        GenericOperatorSparsity <: Function,
        Data,
    } <: AbstractPhysics
    """
    Flux between neighboring control volumes: `flux(f,u,edge,data)`
    should return in `f[i]` the flux of species i along the edge joining circumcenters
    of neighboring control volumes.  u is a  2D array such that for species i,
    `u[i,1]` and `u[i,2]` contain the unknown values at the corresponding ends of the edge.
    """
    flux::Flux

    """
    Storage term (term under time derivative): `storage(f,u,node,data)`

    It should return in `f[i]` the storage term for the i-th equation. `u[i]` contains the value of
    the i-th unknown.
    """
    storage::Storage

    """
    Reaction term: `reaction(f,u,node,data)`

    It should return in `f[i]` the reaction term for the i-th equation. `u[i]` contains the value of
    the i-th unknown.
    """
    reaction::Reaction

    """
    Edge reaction term: `edgereaction(f,u,edge,data)`

    It should return in `f[i]` the reaction term for the i-th equation. `u[i]` contains the value of
    the i-th unknown.  u is a  2D array such that for species i,
    `u[i,1]` and `u[i,2]` contain the unknown values at the corresponding ends of the edge.
    """
    edgereaction::EdgeReaction

    """
    Source term: `source(f,node,data)`.

    It should return the in `f[i]` the value of the source term for the i-th equation.
    """
    source::Source

    """
    Flux between neighboring control volumes on the boundary. Called on edges
    fully adjacent to the boundary: `bflux(f,u,bedge,data)
    """
    bflux::BFlux

    """
    Boundary reaction term: `breaction(f,u,node,data)`
    Similar to reaction, but restricted to the inner or outer boundaries.
    """
    breaction::BReaction

    """
    Boundary source term: `bsource(f,node,data)`.

    It should return in `f[i]` the value of the source term for the i-th equation.
    """
    bsource::BSource

    """
    Boundary storage term: `bstorage(f,u,node,data)`
    Similar to storage, but restricted to the inner or outer boundaries.
    """
    bstorage::BStorage

    """
    Outflow boundary term  `boutflow(f,u,edge,data)`
    This function is called for edges (including interior ones) which have at least one ode
    on one of the outflow boundaries. Within this function,
    [`outflownode`](@ref) and  [`isoutflownode`](@ref) can be used to identify
    that node.  There is some ambiguity in the case that both nodes are outflow
    nodes, in that case it is assumed that the contribution is zero.
    """
    boutflow::BOutflow

    """
    List (Vector) of boundary regions which carry outflow boundary conditions.
    Influences when `boutflow` is called.
    """
    outflowboundaries::Vector{Int}

    """
    Generic operator  `generic_operator(f,u,sys)`. 
    This operator acts on the full solution `u` of a system. Sparsity
    is detected automatically  unless `generic_operator_sparsity` is given.
    """
    generic_operator::GenericOperator

    """
    Function defining the sparsity structure of the generic operator.
    This should return the sparsity pattern of the `generic_operator`.
    """
    generic_operator_sparsity::GenericOperatorSparsity

    """
    User data (parameters).
    This allows to pass various parameters to the callback functions.
    """
    data::Data

    """
    Number of species including boundary species. Meaningless & deprecated.
    """
    num_species::Int8
end

##########################################################

"""
````
Physics(;num_species=0,
         data=nothing,
         flux,
         reaction,
         edgereaction,
         storage,
         source,
         breaction,
         bstorage,
         boutflow,
         outflowboundaries,
         generic,
         generic_sparsity
    )
````

Constructor for physics data. For the meaning of the optional keyword arguments, see [`VoronoiFVM.System(grid::ExtendableGrid; kwargs...)`](@ref).
"""
function Physics(;
        num_species = 0,
        data = nothing,
        flux::Function = nofunc,
        reaction::Function = nofunc,
        edgereaction::Function = nofunc,
        storage::Function = default_storage,
        source::Function = nofunc,
        bflux::Function = nofunc,
        breaction::Function = nofunc,
        bsource::Function = nofunc,
        bstorage::Function = nofunc,
        boutflow::Function = nofunc,
        outflowboundaries::Vector{Int} = Int[],
        generic::Function = nofunc,
        generic_sparsity::Function = nofunc,
        kwargs...
    )
    return Physics(
        flux,
        storage,
        reaction,
        edgereaction,
        source,
        bflux,
        breaction,
        bsource,
        bstorage,
        boutflow,
        outflowboundaries,
        generic,
        generic_sparsity,
        data,
        Int8(num_species)
    )
end

function Physics(physics::Physics, data)
    return Physics(
        physics.flux,
        physics.storage,
        physics.reaction,
        physics.edgereaction,
        physics.source,
        physics.bflux,
        physics.breaction,
        physics.bsource,
        physics.bstorage,
        physics.boutflow,
        physics.outflowboundaries,
        physics.generic,
        physics.generic_sparsity,
        data,
        physics.num_species
    )
end


"""
$(SIGNATURES)

Check if physics object has outflow callback
"""
hasoutflow(physics::Physics) = physics.boutflow != nofunc

"""
$(SIGNATURES)

Show physics object
"""
function Base.show(io::IO, physics::AbstractPhysics)
    str = "Physics("
    if !isnothing(physics.data)
        str = str * "data=$(typeof(physics.data)), "
    end

    # function addfunc(func, name)
    #     if func != nofunc
    #         str = str * ", $(name)=$(nameof(func))"
    #     end
    # end

    for name in fieldnames(typeof(physics))
        if (name != :num_species) && (name != :data)  && (name != :outflowboundaries)  && getfield(physics, name) != nofunc
            str = str * "$(name)=$(nameof(getfield(physics, name))), "
        end
    end
    if length(physics.outflowboundaries) > 0
        str = str * "outflowboundaries=$(physics.outflowboundaries)"
    end
    str = str * ")"
    return print(io, str)
end

"""
    $(TYPEDEF)

Abstract type for evaluator.
"""
abstract type AbstractEvaluator end

"""
    $(TYPEDEF)

Evaluator for functions from physics. Allows to call different types of physic functions (flux, reaction, source)
an provides  a common interface to different function formats (with data, without data etc.)

$(TYPEDFIELDS)
"""
struct ResEvaluator{Tv <: Number, Func <: Function, G} <: AbstractEvaluator
    """ wrapper function in Format ready for Diffetential equations"""
    fwrap::Func
    """ pre-allocated result """
    y::Vector{Tv}
    """ Geometry object # geometry (node, edge...)"""
    geom::G
    """ number of species """
    nspec::Int

    """ Is the function not nofunc """
    isnontrivial::Bool
end

"""
     ResEvaluator(physics,symb,uproto,geom,nspec)

Constructor for ResEvaluator
- `physics` Physics object
- `symb`: symbol naming one of the functions in physics to be wrapped.
- `uproto`: solution vector prototype,
- `geom`: node, edge...
- `nspec`: number of species
"""
function ResEvaluator(physics, data, symb::Symbol, uproto::Vector{Tv}, geom, nspec::Int) where {Tv}
    func = getproperty(physics, symb)

    # source functions need special handling here
    if symb == :source || symb == :bsource
        fwrap = function (y)
            y .= 0
            func(rhs(geom, y), geom, data)
            return nothing
        end
    else   # Normal functions wihth u as parameter
        fwrap = function (y, u)
            y .= 0
            ## for ii in ..  uu[geom.speclist[ii]]=u[ii]
            func(rhs(geom, y), unknowns(geom, u), geom, data)
            ## for ii in .. y[ii]=y[geom.speclist[ii]]
            return nothing
        end
    end
    isnontrivial = (func != nofunc)
    y = zeros(Tv, nspec)
    return ResEvaluator(fwrap, y, geom, nspec, isnontrivial)
end

"""
$(TYPEDSIGNATURES)

Call function in evaluator, store result in predefined memory.
"""
function evaluate!(e::ResEvaluator, u)
    e.isnontrivial ? e.fwrap(e.y, u) : nothing
    return nothing
end

"""
$(TYPEDSIGNATURES)

Call function in evaluator, store result in predefined memory.
"""
function evaluate!(e::ResEvaluator)
    e.isnontrivial ? e.fwrap(e.y) : nothing
    return nothing
end

"""
$(TYPEDSIGNATURES)

Retrieve evaluation result
"""
res(e::ResEvaluator) = e.y

"""
    $(TYPEDEF)

Evaluator for functions from physics and their Jacobians. Allows to call different types of physic functions (flux, reaction, source)
an provides  a common interface to different function formats (with data, without data etc.)

$(TYPEDFIELDS)
"""
struct ResJacEvaluator{Tv <: Number, Func <: Function, Cfg, Res, G} <: AbstractEvaluator
    """ wrapper function in Format ready for Differential equations"""
    fwrap::Func
    """ ForwardDiff.JacobianConfig """
    config::Cfg
    """ DiffResults.JacobianResult"""
    result::Res
    """ pre-allocated result """
    y::Vector{Tv}
    """ Geometry object # geometry (node, edge...)"""
    geom::G
    """ number of species """
    nspec::Int
    """ Is the function not nofunc """
    isnontrivial::Bool
end

"""
    $(SIGNATURES)

Constructor for ResJEvaluator
- `physics` Physics object
- `symb`: symbol naming one of the functions in physics to be wrapped.
- `uproto`: solution vector prototype,
- `geom`: node, edge...
- `nspec`: number of species
"""
function ResJacEvaluator(physics, data, symb::Symbol, uproto::Vector{Tv}, geom, nspec) where {Tv}
    func = getproperty(physics, symb)

    fwrap = function (y, u)
        y .= 0
        ## for ii in ..  uu[geom.speclist[ii]]=u[ii]
        func(rhs(geom, y), unknowns(geom, u), geom, data)
        ## for ii in .. y[ii]=y[geom.speclist[ii]]
        return nothing
    end

    isnontrivial = (func != nofunc)

    y = zeros(Tv, nspec)
    u = zeros(Tv, length(uproto))
    jac = zeros(Tv, nspec, length(u))
    result = DiffResults.DiffResult(u, jac)
    config = ForwardDiff.JacobianConfig(fwrap, y, u, ForwardDiff.Chunk(u, length(u)))
    return ResJacEvaluator(fwrap, config, result, y, geom, nspec, isnontrivial)
end

"""
$(TYPEDSIGNATURES)

Call function in evaluator, store result and jacobian in predefined memory.
"""
function evaluate!(e::ResJacEvaluator, u)
    if e.isnontrivial
        ForwardDiff.jacobian!(e.result, e.fwrap, e.y, u, e.config)
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Retrieve evaluation result
"""
res(e::ResJacEvaluator) = DiffResults.value(e.result)

"""
$(TYPEDSIGNATURES)

Retrieve Jacobian
"""
jac(e::ResJacEvaluator) = DiffResults.jacobian(e.result)

"""
$(TYPEDSIGNATURES)

Does calling the evaluator giva nontrivial (nonzero) result?
"""
isnontrivial(e::AbstractEvaluator) = e.isnontrivial

##########################################################

# "Generate" a flux function
function diffusion_flux(D::T) where {T}
    return (y, u, args...) -> y[1] = D(u[1, 1] + u[1, 2]) * (u[1, 1] - u[1, 2])
end

"""
     boundary_dirichlet!(y,u,bnode,ispec,ireg,val)

Set Dirichlet boundary condition for species ispec at boundary ibc.
"""
function boundary_dirichlet!(y, u, bnode::AbstractGeometryItem, ispec, ireg, val; penalty = bnode.Dirichlet)
    if bnode.region == ireg
        y[ispec] += penalty * (u[ispec] - val)
        # just for call during initialization, so we can convert from dual number
        bnode.dirichlet_value[ispec] = value(val)
    end
    return nothing
end

"""
     boundary_dirichlet!(y,u,bnode; kwargs...)

Keyword argument version:
- `species`: species number. Default: 1
- `region`: boundary region number. By default, all boundary regions.
- `value`: value 
"""
function boundary_dirichlet!(y, u, bnode::AbstractGeometryItem; species = 1, region = bnode.region, value = 0, penalty = bnode.Dirichlet)
    return boundary_dirichlet!(y, u, bnode, species, region, value; penalty)
end

"""
       ramp(t; kwargs...)
Ramp function for specifying time dependent boundary conditions

Keyword arguments:
- `dt`: Tuple: start and end time of ramp. Default: `(0,0.1)`
- `du`: Tuple: values at start and end time. Default: `(0,0)`
"""
function ramp(t; dt = (0, 0.1), du = (0, 0))
    (t, ubegin, uend, tbegin, tend) = promote(Float64(t), du[1], du[2], dt[1], dt[2])
    if t < tbegin
        return ubegin
    elseif t < tend
        return ubegin + (uend - ubegin) * (t - tbegin) / (tend - tbegin)
    else
        return uend
    end
end


"""
     boundary_neumann!(y,u,bnode,ispec,ireg,val)

Set Neumann boundary condition for species ispec at boundary ibc.
"""
boundary_neumann!(y, u, bnode::AbstractGeometryItem, ispec, ireg, val) = bnode.region == ireg ? y[ispec] -= val : nothing

"""
     boundary_neumann!(y,u,bnode; kwargs...)
Keyword argument version:
- `species`: species number. Default: 1
- `region`: boundary region number. By default, all boundary regions.
- `value`: value
"""
function boundary_neumann!(y, u, bnode::AbstractGeometryItem; species = 1, region = bnode.region, value = 0)
    return boundary_neumann!(y, u, bnode, species, region, value)
end


"""
     boundary_robin!(y,u,bnode,ispec,ireg,fac,val)

Set Robin boundary condition for species ispec at boundary ibc.
"""
boundary_robin!(y, u, bnode::AbstractGeometryItem, ispec, ireg, fac, val) = bnode.region == ireg ? y[ispec] += fac * u[ispec] - val : nothing

"""
     boundary_robin!(y,u,bnode, args...; kwargs...)
Keyword argument version:
- `species`: species number. Default: 1
- `region`: boundary region number. By default, all boundary regions.
- `factor`: factor
- `value`: value
"""
function boundary_robin!(y, u, bnode::AbstractGeometryItem; species = 1, region = bnode.region, factor = 0, value = 0)
    return boundary_robin!(y, u, bnode, species, region, factor, value)
end
