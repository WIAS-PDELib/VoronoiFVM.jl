################################################
"""
$(TYPEDEF)

Data structure containing DenseSystem used to calculate
test functions for boundary flux calculations.


$(TYPEDFIELDS)
"""
mutable struct TestFunctionFactory
    """
    Original system
    """
    system::AbstractSystem

    """
    Test function system
    """
    tfsystem::DenseSystem

    """
    Solver control
    """
    control::SolverControl
end

################################################
"""
$(TYPEDSIGNATURES)

Constructor for TestFunctionFactory from System
"""
function TestFunctionFactory(system::AbstractSystem; control = SolverControl())
    physics = Physics(; flux = function (f, u, edge)
                          f[1] = u[1] - u[2]
                      end,
                      storage = function (f, u, node)
                          f[1] = u[1]
                      end)
    tfsystem = System(system.grid, physics; unknown_storage = :dense)
    enable_species!(tfsystem, 1, [i for i = 1:num_cellregions(system.grid)])
    return TestFunctionFactory(system, tfsystem, control)
end

############################################################################
"""
$(TYPEDSIGNATURES)

Create testfunction which has Dirichlet zero boundary conditions  for boundary
regions in bc0 and Dirichlet one boundary conditions  for boundary
regions in bc1.
"""
function testfunction(factory::TestFunctionFactory, bc0, bc1)
    u = unknowns(factory.tfsystem)
    f = unknowns(factory.tfsystem)
    u .= 0
    f .= 0

    factory.tfsystem.boundary_factors .= 0
    factory.tfsystem.boundary_values .= 0

    for i = 1:length(bc1)
        factory.tfsystem.boundary_factors[1, bc1[i]] = Dirichlet
        factory.tfsystem.boundary_values[1, bc1[i]] = -1
    end

    for i = 1:length(bc0)
        factory.tfsystem.boundary_factors[1, bc0[i]] = Dirichlet
        factory.tfsystem.boundary_values[1, bc0[i]] = 0
    end
    _complete!(factory.tfsystem)
    eval_and_assemble(factory.tfsystem, u, u, f, Inf, Inf, 0.0, zeros(0))

    _initialize!(u, factory.tfsystem)

    method_linear = factory.control.method_linear
    if isnothing(method_linear)
        method_linear = UMFPACKFactorization()
    end

    p = LinearProblem(factory.tfsystem.matrix, vec(f))
    sol = solve(p, method_linear)
    sol.u
end

############################################################################

"""
$(SIGNATURES)

Calculate test function integral for transient solution.
"""
function integrate(system::AbstractSystem, tf, U::AbstractMatrix{Tv},
                   Uold::AbstractMatrix{Tv}, tstep; params = Tv[]) where {Tv}
    grid = system.grid
    nspecies = num_species(system)
    integral = zeros(Tv, nspecies)
    tstepinv = 1.0 / tstep

    nparams = system.num_parameters
    @assert nparams == length(params)

    # !!! params etc 
    physics = system.physics
    node = Node(system, 0.0, 1.0, params)
    bnode = BNode(system, 0.0, 1.0, params)
    edge = Edge(system, 0.0, 1.0, params)
    bedge = Edge(system, 0.0, 1.0, params)

    UKL = Array{Tv, 1}(undef, 2 * nspecies + nparams)
    UK = Array{Tv, 1}(undef, nspecies + nparams)
    UKold = Array{Tv, 1}(undef, nspecies + nparams)

    if nparams > 0
        UK[(nspecies + 1):end] .= params
        UKold[(nspecies + 1):end] .= params
        UKL[(2 * nspecies + 1):end] .= params
    end

    src_eval = ResEvaluator(physics, :source, UK, node, nspecies + nparams)
    rea_eval = ResEvaluator(physics, :reaction, UK, node, nspecies + nparams)
    erea_eval = ResEvaluator(physics, :edgereaction, UK, edge, nspecies + nparams)
    stor_eval = ResEvaluator(physics, :storage, UK, node, nspecies + nparams)
    storold_eval = ResEvaluator(physics, :storage, UKold, node, nspecies + nparams)
    flux_eval = ResEvaluator(physics, :flux, UKL, edge, nspecies + nparams)

    geom = grid[CellGeometries][1]

    for icell = 1:num_cells(grid)
        for iedge = 1:num_edges(geom)
            _fill!(edge, iedge, icell)

            @views UKL[1:nspecies] .= U[:, edge.node[1]]
            @views UKL[(nspecies + 1):(2 * nspecies)] .= U[:, edge.node[2]]

            evaluate!(flux_eval, UKL)
            flux = res(flux_eval)

            function asm_res(idofK, idofL, ispec)
                integral[ispec] += system.celledgefactors[iedge, icell] * flux[ispec] * (tf[edge.node[1]] - tf[edge.node[2]])
            end
            assemble_res(edge, system, asm_res)

            if isnontrivial(erea_eval)
                evaluate!(erea_eval, UKL)
                erea = res(erea_eval)

                function easm_res(idofK, idofL, ispec)
                    integral[ispec] += system.celledgefactors[iedge, icell] * erea[ispec] * (tf[edge.node[1]] + tf[edge.node[2]])
                end
                assemble_res(edge, system, easm_res)
            end
        end

        for inode = 1:num_nodes(geom)
            _fill!(node, inode, icell)
            for ispec = 1:nspecies
                UK[ispec] = U[ispec, node.index]
                UKold[ispec] = Uold[ispec, node.index]
            end

            evaluate!(rea_eval, UK)
            rea = res(rea_eval)
            evaluate!(stor_eval, UK)
            stor = res(stor_eval)
            evaluate!(storold_eval, UKold)
            storold = res(storold_eval)
            evaluate!(src_eval)
            src = res(src_eval)

            function asm_res(idof, ispec)
                integral[ispec] += system.cellnodefactors[inode, icell] *
                                   (rea[ispec] - src[ispec] + (stor[ispec] - storold[ispec]) * tstepinv) * tf[node.index]
            end
            assemble_res(node, system, asm_res)
        end
    end
    return integral
end

############################################################################
"""
$(SIGNATURES)

Calculate test function integral for steady state solution.
"""
function integrate(system::AbstractSystem, tf::Vector{Tv}, U::AbstractMatrix{Tu}; kwargs...) where {Tu, Tv}
    integrate(system, tf, U, U, Inf; kwargs...)
end

############################################################################
"""
$(SIGNATURES)

Steady state part of test function integral.
"""
function integrate_stdy(system::AbstractSystem, tf::Vector{Tv}, U::AbstractArray{Tu, 2}) where {Tu, Tv}
    grid = system.grid
    nspecies = num_species(system)
    integral = zeros(Tu, nspecies)

    physics = system.physics
    node = Node(system)
    bnode = BNode(system)
    edge = Edge(system)
    bedge = BEdge(system)

    UKL = Array{Tu, 1}(undef, 2 * nspecies)
    UK = Array{Tu, 1}(undef, nspecies)
    geom = grid[CellGeometries][1]

    src_eval = ResEvaluator(physics, :source, UK, node, nspecies)
    rea_eval = ResEvaluator(physics, :reaction, UK, node, nspecies)
    erea_eval = ResEvaluator(physics, :edgereaction, UK, node, nspecies)
    flux_eval = ResEvaluator(physics, :flux, UKL, edge, nspecies)

    for icell = 1:num_cells(grid)
        for iedge = 1:num_edges(geom)
            _fill!(edge, iedge, icell)

            @views UKL[1:nspecies] .= U[:, edge.node[1]]
            @views UKL[(nspecies + 1):(2 * nspecies)] .= U[:, edge.node[2]]
            evaluate!(flux_eval, UKL)
            flux = res(flux_eval)

            function asm_res(idofK, idofL, ispec)
                integral[ispec] += system.celledgefactors[iedge, icell] * flux[ispec] * (tf[edge.node[1]] - tf[edge.node[2]])
            end
            assemble_res(edge, system, asm_res)

            if isnontrivial(erea_eval)
                evaluate!(erea_eval, UKL)
                erea = res(erea_eval)

                function easm_res(idofK, idofL, ispec)
                    integral[ispec] += system.celledgefactors[iedge, icell] * erea[ispec] * (tf[edge.node[1]] + tf[edge.node[2]])
                end
                assemble_res(edge, system, easm_res)
            end
        end

        for inode = 1:num_nodes(geom)
            _fill!(node, inode, icell)

            @views UK .= U[:, node.index]

            evaluate!(rea_eval, UK)
            rea = res(rea_eval)
            evaluate!(src_eval)
            src = res(src_eval)

            function asm_res(idof, ispec)
                integral[ispec] += system.cellnodefactors[inode, icell] * (rea[ispec] - src[ispec]) * tf[node.index]
            end
            assemble_res(node, system, asm_res)
        end
    end
    return integral
end

############################################################################
"""
$(SIGNATURES)

Calculate transient part of test function integral.
"""
function integrate_tran(system::AbstractSystem, tf::Vector{Tv}, U::AbstractArray{Tu, 2}) where {Tu, Tv}
    grid = system.grid
    nspecies = num_species(system)
    integral = zeros(Tu, nspecies)

    physics = system.physics
    node = Node(system)
    bnode = BNode(system)
    edge = Edge(system)
    bedge = BEdge(system)
    # !!! Parameters

    UK = Array{Tu, 1}(undef, nspecies)
    geom = grid[CellGeometries][1]
    csys = grid[CoordinateSystem]
    stor_eval = ResEvaluator(physics, :storage, UK, node, nspecies)

    for icell = 1:num_cells(grid)
        for inode = 1:num_nodes(geom)
            _fill!(node, inode, icell)
            @views UK .= U[:, node.index]
            evaluate!(stor_eval, UK)
            stor = res(stor_eval)

            asm_res(idof, ispec) = integral[ispec] += system.cellnodefactors[inode, icell] * stor[ispec] * tf[node.index]
            assemble_res(node, system, asm_res)
        end
    end
    return integral
end

# function checkdelaunay(grid)
#     nreg=num_cellregions(grid)
#     D=ones(nreg)

#     physics=Physics( 
#         flux=function(f,u,edge)
#         f[1]=Du[1]-u[2]
#         end,
#         storage=function(f,u,node)
#         f[1]=0
#         end
#     )
#     tfsystem=System(system.grid,physics,unknown_storage=:dense)
#     enable_species!(tfsystem,1,[i for i=1:num_cellregions(system.grid)])
#     return TestFunctionFactory{Tv}(system,tfsystem)
# end
