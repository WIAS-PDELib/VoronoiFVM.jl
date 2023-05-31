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
The result is an `nspec x nregion` vector.
"""
function integrate(system::AbstractSystem{Tv, Tc, Ti, Tm}, F::Function, U::AbstractMatrix{Tu};
                   boundary = false) where {Tu, Tv, Tc, Ti, Tm}
    grid = system.grid
    data = system.physics.data
    nspecies = num_species(system)
    res = zeros(Tu, nspecies)

    if boundary
        bnode = BNode(system)
        bnodeparams = (bnode,)
        if isdata(data)
            bnodeparams = (bnode, data)
        end
        #!!!        bnode.time=time
        #!!!        bnode.embedparam=embedparam

        bfaceregions = grid[BFaceRegions]
        nbfaceregions = maximum(bfaceregions)
        integral = zeros(Tu, nspecies, nbfaceregions)

        for item in nodebatch(system.boundary_assembly_data)
            for ibnode in noderange(system.boundary_assembly_data,item)
                _fill!(bnode,system.boundary_assembly_data,ibnode,item)
                res .= zero(Tv)
                @views F(rhs(bnode, res), unknowns(bnode, U[:, bnode.index]), bnodeparams...)
                asm_res(idof, ispec) = integral[ispec, bnode.region] += bnode.fac * res[ispec]
                assemble_res(bnode, system, asm_res)
            end
        end
    else
        node = Node(system)
        nodeparams = (node,)
        if isdata(data)
            nodeparams = (node, data)
        end
        #!!!        node.time=time
        #!!!        node.embedparam=embedparam
        cellregions = grid[CellRegions]
        ncellregions = maximum(cellregions)
        integral = zeros(Tu, nspecies, ncellregions)
        for item in nodebatch(system.assembly_data)
            for inode in noderange(system.assembly_data,item)
                _fill!(node,system.assembly_data,inode,item)
                res .= zero(Tv)
                @views F(rhs(node, res), unknowns(node, U[:, node.index]), nodeparams...)
                asm_res(idof, ispec) = integral[ispec, node.region] += node.fac * res[ispec]
                assemble_res(node, system, asm_res)
            end
        end
    end

    return SolutionIntegral(integral)
end

"""
    edgeintegrate(system,F,U; boundary=false)    

Integrate edge function (same signature as flux function)
 `F` of  solution vector region-wise over domain or boundary.
The result is an `nspec x nregion` vector.
"""
function edgeintegrate(system::AbstractSystem{Tv, Tc, Ti, Tm}, F::Function, U::AbstractMatrix{Tu};
                       boundary = false) where {Tu, Tv, Tc, Ti, Tm}
    grid = system.grid
    dim=dim_space(grid)
    data = system.physics.data
    nspecies = num_species(system)
    res = zeros(Tu, nspecies)
    nparams = system.num_parameters

    UKL = Array{Tv, 1}(undef, 2 * nspecies + nparams)

    if boundary
        error("missing implementation of boundary edge integrals")
    else
        edge = Edge(system)
        edgeparams = (edge,)
        if isdata(data)
            edgeparams = (edge, data)
        end
        cellregions = grid[CellRegions]
        ncellregions = maximum(cellregions)
        integral = zeros(Tu, nspecies, ncellregions)
        for item in edgebatch(system.assembly_data)
            for iedge in edgerange(system.assembly_data,item)
                _fill!(edge,system.assembly_data,iedge,item) 
                @views UKL[1:nspecies] .= U[:, edge.node[1]]
                @views UKL[(nspecies+1):(2*nspecies)] .= U[:, edge.node[2]]
                res .= zero(Tv)
                @views F(rhs(edge, res), unknowns(edge,UKL), edgeparams...)
                function asm_res(idofK, idofL, ispec)
                    h=meas(edge)
                    # This corresponds to the multiplication with the diamond volume.
                    integral[ispec,edge.region] += h^2*edge.fac*res[ispec]/dim
                end
                assemble_res(edge, system, asm_res)
            end
        end
    end
    return SolutionIntegral(integral)
end


############################################################################
"""
$(SIGNATURES)

Reconstruction of  edge flux as  vector function  on the nodes  of the
triangulation.  The result  can be seen as a  piecewiesw linear vector
function in the FEM space spanned by the discretization nodes exhibiting
the flux density.  

The reconstruction is based on the  "magic formula"
R. Eymard, T. Gallouet, R. Herbin, IMA Journal of Numerical Analysis (2006)
26, 326−353, Lemma 2.4 .

The return value is a `dim x nspec x nnodes` vector. The flux of species i
can  e.g. plotted via GridVisualize.vectorplot.

Example:
```julia
    ispec=3
    vis=GridVisualizer(Plotter=Plotter)
    scalarplot!(vis,grid,solution[ispec,:],clear=true,colormap=:summer)
    vectorplot!(vis,grid,nf[:,ispec,:],clear=false)
    reveal(vis)
```

CAVEAT: there is a possible unsolved problem with the values at domain 
corners in the code. Please see any potential boundary artifacts as a manifestation
of this issue and report them.
"""
function nodeflux(system::AbstractSystem{Tv, Tc, Ti, Tm}, U::AbstractArray{Tu, 2}) where {Tu, Tv, Tc, Ti, Tm}
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

    
    # !!! TODO Parameter handling here
    UKL = Array{Tu, 1}(undef, 2 * nspecies)
    flux_eval = ResEvaluator(system.physics, :flux, UKL, edge, nspecies)

    for item in edgebatch(system.assembly_data)
        for iedge in edgerange(system.assembly_data,item)
            _fill!(edge,system.assembly_data,iedge,item) 
            K = edge.node[1]
            L = edge.node[2]
            @views UKL[1:nspecies] .= U[:, edge.node[1]]
            @views UKL[(nspecies + 1):(2 * nspecies)] .= U[:, edge.node[2]]
            evaluate!(flux_eval, UKL)
            edgeflux = res(flux_eval)
            function asm_res(idofK, idfoL, ispec)
                @views nodeflux[:, ispec, K] .+= edge.fac * edgeflux[ispec] * (xsigma[:, edge.index] - coord[:, K])
                @views nodeflux[:, ispec, L] .-= edge.fac * edgeflux[ispec] * (xsigma[:, edge.index] - coord[:, L])
            end
            assemble_res(edge, system, asm_res)
        end
    end
    
    for item in nodebatch(system.assembly_data)
        for inode in noderange(system.assembly_data,item)
            _fill!(node,system.assembly_data,inode,item)        
            nodevol[node.index] += node.fac
        end
    end
    
    for inode = 1:nnodes
        @views nodeflux[:, :, inode] /= nodevol[inode]
    end
    nodeflux
end


#####################################################################################
"""
    $(SIGNATURES)

Calculate Euklidean norm of the degree of freedom vector.
"""
function LinearAlgebra.norm(system::AbstractSystem, u, p::Number=2) end

function LinearAlgebra.norm(system::DenseSystem, u, p::Number=2)
    _initialize_inactive_dof!(u, system)
    norm(u, p)
end


LinearAlgebra.norm(system::SparseSystem, u::SparseSolutionArray, p::Number=2) = LinearAlgebra.norm(u.node_dof.nzval, p)

"""
    $(SIGNATURES)

Calculate weighted discrete ``L^p`` norm of a solution vector.
"""
function lpnorm(sys,u,p,species_weights=ones(num_species(sys)))
    nspec=num_species(sys)
    II=integrate(sys, (y,u,node,data=nothing)->y.=u.^p,u)
    (sum([sum(II[i,:]) for i=1:nspec].*species_weights))^(1/p)
end

"""
    $(SIGNATURES)

Calculate weigthed discrete ``L^2(\\Omega)`` norm of a solution vector. 
"""
function l2norm(sys,u,species_weights=ones(num_species(sys)))
    lpnorm(sys,u,2,species_weights)
end

"""
    $(SIGNATURES)

Calculate weighted discrete ``W^{1,p}(\\Omega)`` seminorm of a solution vector.
"""
function w1pseminorm(sys,u,p,species_weights=ones(num_species(sys)))
    nspec=num_species(sys)
    dim=dim_space(sys.grid)
    function f(y,u,edge,data=nothing)
        h=meas(edge)
        for ispec=1:nspec
            y[ispec]=dim*((u[ispec,1]-u[ispec,2])/h)^p
        end
    end
    II=edgeintegrate(sys,f,u)
    (sum([sum(II[i,:]) for i=1:nspec].*species_weights))^(1/p)
end

"""
    $(SIGNATURES)

Calculate weighted discrete ``H^1(\\Omega)`` seminorm of a solution vector.
"""
function h1seminorm(sys,u,species_weights=ones(num_species(sys)))
    w1pseminorm(sys,u,2,species_weights)
end

"""
    $(SIGNATURES)

Calculate weighted discrete ``W^{1,p}(\\Omega)`` norm of a solution vector.
"""
function w1pnorm(sys,u,p,species_weights=ones(num_species(sys)))
    lpnorm(sys,u,p,species_weights)+w1pseminorm(sys,u,p,species_weights)
end

"""
    $(SIGNATURES)

Calculate weighted discrete ``H^1(\\Omega)`` norm of a solution vector.
"""
function h1norm(sys,u,species_weights=ones(num_species(sys)))
    w1pnorm(sys,u,2,species_weights)
end

function _bochnernorm(sys,u::TransientSolution,p,species_weights,spatialnorm::F) where F
    n=length(u.t)
    nrm=spatialnorm(sys,u.u[1],p,species_weights)^p*(u.t[2]-u.t[1])/2
    for i=2:n-1
        nrm+=spatialnorm(sys,u.u[i],p,species_weights)^p*(u.t[i+1]-u.t[i-1])/2
    end
    nrm+=spatialnorm(sys,u.u[n],p,species_weights)^p*(u.t[end]-u.t[end-1])/2
    nrm^(1/p)
end

"""
    $(SIGNATURES)

Calculate weighted discrete ``L^p([0,T];W^{1,p}(\\Omega))`` norm of a transient solution.
"""
function lpw1pnorm(sys,u::TransientSolution,p,species_weights=ones(num_species(sys)))
    _bochnernorm(sys,u,p,species_weights,w1pnorm)
end


"""
    $(SIGNATURES)

Calculate weighted discrete ``L^p([0,T];W^{1,p}(\\Omega))`` seminorm of a transient solution.
"""
function lpw1pseminorm(sys,u::TransientSolution,p,species_weights=ones(num_species(sys)))
    _bochnernorm(sys,u,p,species_weights,w1pseminorm)
end



"""
    $(SIGNATURES)

Calculate weighted discrete ``L^2([0,T];H^1(\\Omega))`` seminorm of a transient solution.
"""
function l2h1seminorm(sys,u::TransientSolution,species_weights=ones(num_species(sys)))
    lpw1pseminorm(sys,u,2,species_weights)
end


"""
    $(SIGNATURES)

Calculate weighted discrete ``L^2([0,T];H^1(\\Omega))`` norm of a transient solution.
"""
function l2h1norm(sys,u::TransientSolution,species_weights=ones(num_species(sys)))
    lpw1pnorm(sys,u,2,species_weights)
end


