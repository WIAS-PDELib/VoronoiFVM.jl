struct SolutionIntegral{T} <: AbstractMatrix{T}
    value::Matrix{T}
end

Base.size(I::SolutionIntegral) = size(I.value)
Base.getindex(I::SolutionIntegral, ispec::Integer, ireg) = I.value[ispec, ireg]
Base.setindex!(I::SolutionIntegral, v, ispec::Integer, ireg) = I.value[ispec, ireg] = v

################################################################
"""
    integrate(system,F,U; boundary=false)

Integrate node function (same signature as reaction or storage)
 `F` of  solution vector region-wise over domain or boundary.
The result is an `nspec x nregion` matrix.
"""
function integrate(
        system::AbstractSystem{Tv, Tc, Ti, Tm}, F::Tf, U::AbstractMatrix{Tu};
        boundary = false, data = system.physics.data
    ) where {Tu, Tv, Tc, Ti, Tm, Tf}
    _complete!(system)
    grid = system.grid
    nspecies = num_species(system)
    res = zeros(Tu, nspecies)

    if boundary
        bnode = BNode(system)
        #!!!        bnode.time=time
        #!!!        bnode.embedparam=embedparam

        bfaceregions = grid[BFaceRegions]
        nbfaceregions = maximum(bfaceregions)
        integral = zeros(Tu, nspecies, nbfaceregions)

        for item in nodebatch(system.boundary_assembly_data)
            for ibnode in noderange(system.boundary_assembly_data, item)
                _fill!(bnode, system.boundary_assembly_data, ibnode, item)
                res .= zero(Tv)
                @views F(rhs(bnode, res), unknowns(bnode, U[:, bnode.index]), bnode, data)
                asm_res(idof, ispec) = integral[ispec, bnode.region] += bnode.fac * res[ispec]
                assemble_res(bnode, system, asm_res)
            end
        end
    else
        node = Node(system)
        #!!!        node.time=time
        #!!!        node.embedparam=embedparam
        cellregions = grid[CellRegions]
        ncellregions = maximum(cellregions)
        integral = zeros(Tu, nspecies, ncellregions)
        for item in nodebatch(system.assembly_data)
            for inode in noderange(system.assembly_data, item)
                _fill!(node, system.assembly_data, inode, item)
                res .= zero(Tv)
                @views F(rhs(node, res), unknowns(node, U[:, node.index]), node, data)
                asm_res(idof, ispec) = integral[ispec, node.region] += node.fac * res[ispec]
                assemble_res(node, system, asm_res)
            end
        end
    end

    return SolutionIntegral(integral)
end

"""
    integrate(system,U; boundary=false)

Integrate solution vector region-wise over domain or boundary.
The result is an `nspec x nregion` matrix.
"""
function integrate(system::AbstractSystem, U::AbstractMatrix; kwargs...)
    function f(f, u, node, data = nothing)
        return f .= u
    end
    return integrate(system, f, U; kwargs...)
end


"""
    edgeintegrate(system,F,U; boundary=false)

Integrate edge function (same signature as flux function)
 `F` of  solution vector region-wise over domain or boundary.
The result is an `nspec x nregion` matrix.
"""
function edgeintegrate(
        system::AbstractSystem{Tv, Tc, Ti, Tm}, F::Tf, U::AbstractMatrix{Tu};
        boundary = false, data = system.physics.data
    ) where {Tu, Tv, Tc, Ti, Tm, Tf}
    _complete!(system)
    grid = system.grid
    dim = dim_space(grid)
    nspecies = num_species(system)
    res = zeros(Tu, nspecies)
    nparams = system.num_parameters

    UKL = Array{Tv, 1}(undef, 2 * nspecies + nparams)

    if boundary
        error("missing implementation of boundary edge integrals")
    else
        edge = Edge(system)
        cellregions = grid[CellRegions]
        ncellregions = maximum(cellregions)
        integral = zeros(Tu, nspecies, ncellregions)
        for item in edgebatch(system.assembly_data)
            for iedge in edgerange(system.assembly_data, item)
                _fill!(edge, system.assembly_data, iedge, item)
                @views UKL[1:nspecies] .= U[:, edge.node[1]]
                @views UKL[(nspecies + 1):(2 * nspecies)] .= U[:, edge.node[2]]
                res .= zero(Tv)
                @views F(rhs(edge, res), unknowns(edge, UKL), edge, data)
                function asm_res(idofK, idofL, ispec)
                    h = meas(edge)
                    # This corresponds to the multiplication with the diamond volume.
                    return integral[ispec, edge.region] += h^2 * edge.fac * res[ispec] / dim
                end
                assemble_res(edge, system, asm_res)
            end
        end
    end
    return SolutionIntegral(integral)
end

############################################################################
"""
    nodeflux(system, F,U; data)

Reconstruction of an edge function as  vector function  on the nodes  of the
triangulation.  The result  can be seen as a  piecewiesw linear vector
function in the FEM space spanned by the discretization nodes exhibiting
the flux density.

The reconstruction is based on the  "magic formula"
R. Eymard, T. Gallouët, R. Herbin, IMA Journal of Numerical Analysis (2006)
26, 326−353, Lemma 2.4  (also: [hal.science/hal-00004840](https://hal.science/hal-00004840) ).

The return value is a `dim x nspec x nnodes` vector.

CAVEAT: there is a possible unsolved problem with the values at domain
corners in the code. Please see any potential boundary artifacts as a manifestation
of this issue and report them.
"""
function nodeflux(
        system::AbstractSystem{Tv, Tc, Ti, Tm}, F::Tf, U::AbstractArray{Tu, 2};
        data = system.physics.data
    ) where {Tu, Tv, Tc, Ti, Tm, Tf}
    _complete!(system)
    grid = system.grid
    dim = dim_space(grid)
    nnodes = num_nodes(grid)
    nspecies = num_species(system)
    nodeflux = zeros(Tu, dim, nspecies, nnodes)
    edgeflux = zeros(Tu, nspecies)
    xsigma = grid[VoronoiFaceCenters]
    coord = grid[Coordinates]
    nodevol = zeros(Tv, nnodes)
    cellnodes = grid[CellNodes]
    physics = system.physics
    edge = Edge(system)
    node = Node(system)
    res = zeros(Tu, nspecies)

    # !!! TODO Parameter handling here
    UKL = Array{Tu, 1}(undef, 2 * nspecies)

    for item in edgebatch(system.assembly_data)
        for iedge in edgerange(system.assembly_data, item)
            _fill!(edge, system.assembly_data, iedge, item)
            K = edge.node[1]
            L = edge.node[2]
            @views UKL[1:nspecies] .= U[:, edge.node[1]]
            @views UKL[(nspecies + 1):(2 * nspecies)] .= U[:, edge.node[2]]
            res .= zero(Tv)
            @views F(rhs(edge, res), unknowns(edge, UKL), edge, data)

            function asm_res(idofK, idfoL, ispec)
                @views nodeflux[:, ispec, K] .+= edge.fac * res[ispec] * (xsigma[:, edge.index] - coord[:, K])
                return @views nodeflux[:, ispec, L] .-= edge.fac * res[ispec] * (xsigma[:, edge.index] - coord[:, L])
            end
            assemble_res(edge, system, asm_res)
        end
    end

    for item in nodebatch(system.assembly_data)
        for inode in noderange(system.assembly_data, item)
            _fill!(node, system.assembly_data, inode, item)
            nodevol[node.index] += node.fac
        end
    end

    for inode in 1:nnodes
        @views nodeflux[:, :, inode] /= nodevol[inode]
    end
    return nodeflux
end

"""
    nodeflux(system, U; data)

Reconstruction of the edge flux as vector function  on the nodes  of the
triangulation. See [`nodeflux(system,F,U;data)`](@ref).

 The flux of species i
can  e.g. plotted via GridVisualize.vectorplot.

Example:
```julia
    ispec=3
    vis=GridVisualizer(Plotter=Plotter)
    scalarplot!(vis,grid,solution[ispec,:],clear=true,colormap=:summer)
    vectorplot!(vis,grid,nf[:,ispec,:],clear=false)
    reveal(vis)
```


"""
function nodeflux(
        system::AbstractSystem{Tv, Tc, Ti, Tm}, U::AbstractArray{Tu, 2};
        data = system.physics.data
    ) where {Tu, Tv, Tc, Ti, Tm}
    return nodeflux(system, system.physics.flux, U)
end

#####################################################################################
"""
    $(SIGNATURES)

Calculate Euklidean norm of the degree of freedom vector.
"""
function LinearAlgebra.norm(system::AbstractSystem, u, p::Number = 2) end

function LinearAlgebra.norm(system::DenseSystem, u, p::Number = 2)
    _initialize_inactive_dof!(u, system)
    _complete!(system)
    return norm(u, p)
end

LinearAlgebra.norm(system::SparseSystem, u::SparseSolutionArray, p::Number = 2) = LinearAlgebra.norm(u.u.nzval, p)

"""
    $(SIGNATURES)

Calculate weighted discrete ``L^p`` norm of a solution vector.
"""
function lpnorm(sys, u, p, species_weights = ones(num_species(sys)))
    _complete!(sys)
    nspec = num_species(sys)
    II = integrate(sys, (y, u, node, data = nothing) -> y .= u .^ p, u)
    return (sum([sum(II[i, :]) for i in 1:nspec] .* species_weights))^(1 / p)
end

"""
    $(SIGNATURES)

Calculate weighted discrete ``L^2(\\Omega)`` norm of a solution vector.
"""
function l2norm(sys, u, species_weights = ones(num_species(sys)))
    _complete!(sys)
    return lpnorm(sys, u, 2, species_weights)
end

"""
    $(SIGNATURES)

Calculate weighted discrete ``W^{1,p}(\\Omega)`` seminorm of a solution vector.
"""
function w1pseminorm(sys, u, p, species_weights = ones(num_species(sys)))
    _complete!(sys)
    nspec = num_species(sys)
    dim = dim_space(sys.grid)
    function f(y, u, edge, data = nothing)
        h = meas(edge)
        for ispec in 1:nspec
            y[ispec] = dim * ((u[ispec, 1] - u[ispec, 2]) / h)^p
        end
        return
    end
    II = edgeintegrate(sys, f, u)
    return (sum([sum(II[i, :]) for i in 1:nspec] .* species_weights))^(1 / p)
end

"""
    $(SIGNATURES)

Calculate weighted discrete ``H^1(\\Omega)`` seminorm of a solution vector.
"""
function h1seminorm(sys, u, species_weights = ones(num_species(sys)))
    _complete!(sys)
    return w1pseminorm(sys, u, 2, species_weights)
end

"""
    $(SIGNATURES)

Calculate weighted discrete ``W^{1,p}(\\Omega)`` norm of a solution vector.
"""
function w1pnorm(sys, u, p, species_weights = ones(num_species(sys)))
    _complete!(sys)
    return lpnorm(sys, u, p, species_weights) + w1pseminorm(sys, u, p, species_weights)
end

"""
    $(SIGNATURES)

Calculate weighted discrete ``H^1(\\Omega)`` norm of a solution vector.
"""
function h1norm(sys, u, species_weights = ones(num_species(sys)))
    _complete!(sys)
    return w1pnorm(sys, u, 2, species_weights)
end

function _bochnernorm(sys, u::TransientSolution, p, species_weights, spatialnorm::F) where {F}
    _complete!(sys)
    n = length(u.t)
    nrm = spatialnorm(sys, u.u[1], p, species_weights)^p * (u.t[2] - u.t[1]) / 2
    for i in 2:(n - 1)
        nrm += spatialnorm(sys, u.u[i], p, species_weights)^p * (u.t[i + 1] - u.t[i - 1]) / 2
    end
    nrm += spatialnorm(sys, u.u[n], p, species_weights)^p * (u.t[end] - u.t[end - 1]) / 2
    return nrm^(1 / p)
end

"""
    $(SIGNATURES)

Calculate weighted discrete ``L^p([0,T];W^{1,p}(\\Omega))`` norm of a transient solution.
"""
function lpw1pnorm(sys, u::TransientSolution, p, species_weights = ones(num_species(sys)))
    return _bochnernorm(sys, u, p, species_weights, w1pnorm)
end

"""
    $(SIGNATURES)

Calculate weighted discrete ``L^p([0,T];W^{1,p}(\\Omega))`` seminorm of a transient solution.
"""
function lpw1pseminorm(sys, u::TransientSolution, p, species_weights = ones(num_species(sys)))
    return _bochnernorm(sys, u, p, species_weights, w1pseminorm)
end

"""
    $(SIGNATURES)

Calculate weighted discrete ``L^2([0,T];H^1(\\Omega))`` seminorm of a transient solution.
"""
function l2h1seminorm(sys, u::TransientSolution, species_weights = ones(num_species(sys)))
    return lpw1pseminorm(sys, u, 2, species_weights)
end

"""
    $(SIGNATURES)

Calculate weighted discrete ``L^2([0,T];H^1(\\Omega))`` norm of a transient solution.
"""
function l2h1norm(sys, u::TransientSolution, species_weights = ones(num_species(sys)))
    return lpw1pnorm(sys, u, 2, species_weights)
end

"""
    $(SIGNATURES)

Calculate volumes of Voronoi cells.
"""
function nodevolumes(system)
    _complete!(system)
    if isnothing(system.assembly_data)
        update_grid!(system)
    end
    node = Node(system)
    nodevol = zeros(num_nodes(system.grid))
    for item in nodebatch(system.assembly_data)
        for inode in noderange(system.assembly_data, item)
            _fill!(node, system.assembly_data, inode, item)
            nodevol[node.index] += node.fac
        end
    end
    return nodevol
end

"""
    nondelaunay(grid;tol=1.0e-14, verbose=true)

Return non-Delaunay edges. Returns a vector of tuples:
Each tuple consists of `(node1, node2, edge factor, region)`

If the vector has length 0, the grid is boundary conforming Delaunay
with respect to each cell region. This means that up to `tol`, all
edge form factors are nonnegative.
"""
function nondelaunay(grid; tol = 1.0e-16, verbose = false)
    atol = abs(tol)
    sys = System(grid; unknown_storage = :dense, species = [1], assembly = :edgewise)
    update_grid!(sys)
    (; colptr, rowval, nzval) = sys.assembly_data.edgefactors
    ndelaunay = Tuple{Int, Int, Int, Float64}[]
    enodes = grid[EdgeNodes]
    coord = grid[Coordinates]
    for i in 1:(length(colptr) - 1)
        delaunay = true
        for k in colptr[i]:(colptr[i + 1] - 1)
            if nzval[k] < -atol
                if verbose
                    @warn "Non-Delaunay edge: $(coord[:, enodes[1, i]]), $(coord[:, enodes[2, i]]), region=$(rowval[k])"
                end
                push!(ndelaunay, (enodes[1, i], enodes[2, i], rowval[k], nzval[k]))
            end
        end
    end
    @info "nondelaunay:  $(length(ndelaunay)) non-Delaunay edges detected"
    return ndelaunay
end
