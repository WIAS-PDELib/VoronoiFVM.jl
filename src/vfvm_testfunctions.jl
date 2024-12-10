################################################
"""
$(TYPEDEF)

Data structure containing DenseSystem used to calculate
test functions for boundary flux calculations.

Type parameters:
- `Tu`: value type of test functions
- `Tv`: Default value type of system
$(TYPEDFIELDS)
"""
mutable struct TestFunctionFactory{Tu, Tv}
    """
    Original system
    """
    system::AbstractSystem{Tv}

    """
    Test function system state
    """
    state::SystemState{Tu}

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
function TestFunctionFactory(system::AbstractSystem{Tv}; control = SolverControl()) where {Tv}
    physics = Physics(;
        flux = function (f, u, edge, data)
            return f[1] = u[1] - u[2]
        end,
        storage = function (f, u, node, data)
            return f[1] = u[1]
        end
    )
    tfsystem = System(system.grid, physics; unknown_storage = :dense)
    enable_species!(tfsystem, 1, [i for i in 1:num_cellregions(system.grid)])
    state = SystemState(tfsystem)
    return TestFunctionFactory(system, state, control)
end

############################################################################
"""
$(TYPEDSIGNATURES)

Create testfunction which has Dirichlet zero boundary conditions  for boundary
regions in bc0 and Dirichlet one boundary conditions  for boundary
regions in bc1.
"""
function testfunction(factory::TestFunctionFactory{Tv}, bc0, bc1) where {Tv}
    u = unknowns(factory.state.system)
    f = unknowns(factory.state.system)
    u .= 0
    f .= 0

    factory.state.system.boundary_factors .= 0
    factory.state.system.boundary_values .= 0

    for i in 1:length(bc1)
        factory.state.system.boundary_factors[1, bc1[i]] = Dirichlet(Tv)
        factory.state.system.boundary_values[1, bc1[i]] = -1
    end

    for i in 1:length(bc0)
        factory.state.system.boundary_factors[1, bc0[i]] = Dirichlet(Tv)
        factory.state.system.boundary_values[1, bc0[i]] = 0
    end

    eval_and_assemble(factory.state.system, u, u, f, factory.state.matrix, factory.state.dudp, Inf, Inf, 0.0, nothing, zeros(0))

    _initialize!(u, factory.state.system, nothing)

    method_linear = factory.control.method_linear
    if isnothing(method_linear)
        method_linear = UMFPACKFactorization()
    end

    p = LinearProblem(SparseMatrixCSC(factory.state.matrix), dofs(f))
    sol = solve(p, method_linear)
    return sol.u
end

############################################################################

"""
$(SIGNATURES)

Calculate test function integral for transient solution.
"""
function integrate(
        system::AbstractSystem, tf, U::AbstractMatrix{Tv},
        Uold::AbstractMatrix{Tv}, tstep; params = Tv[], data = system.physics.data
    ) where {Tv}
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

    src_eval = ResEvaluator(physics, data, :source, UK, node, nspecies + nparams)
    rea_eval = ResEvaluator(physics, data, :reaction, UK, node, nspecies + nparams)
    erea_eval = ResEvaluator(physics, data, :edgereaction, UK, edge, nspecies + nparams)
    stor_eval = ResEvaluator(physics, data, :storage, UK, node, nspecies + nparams)
    storold_eval = ResEvaluator(physics, data, :storage, UKold, node, nspecies + nparams)
    flux_eval = ResEvaluator(physics, data, :flux, UKL, edge, nspecies + nparams)

    for item in nodebatch(system.assembly_data)
        for inode in noderange(system.assembly_data, item)
            _fill!(node, system.assembly_data, inode, item)
            for ispec in 1:nspecies
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
                return integral[ispec] += node.fac *
                    (rea[ispec] - src[ispec] + (stor[ispec] - storold[ispec]) * tstepinv) * tf[node.index]
            end
            assemble_res(node, system, asm_res)
        end
    end

    for item in edgebatch(system.assembly_data)
        for iedge in edgerange(system.assembly_data, item)
            _fill!(edge, system.assembly_data, iedge, item)
            @views UKL[1:nspecies] .= U[:, edge.node[1]]
            @views UKL[(nspecies + 1):(2 * nspecies)] .= U[:, edge.node[2]]

            evaluate!(flux_eval, UKL)
            flux = res(flux_eval)

            function asm_res(idofK, idofL, ispec)
                return integral[ispec] += edge.fac * flux[ispec] * (tf[edge.node[1]] - tf[edge.node[2]])
            end
            assemble_res(edge, system, asm_res)

            if isnontrivial(erea_eval)
                evaluate!(erea_eval, UKL)
                erea = res(erea_eval)

                function easm_res(idofK, idofL, ispec)
                    return integral[ispec] += edge.fac * erea[ispec] * (tf[edge.node[1]] + tf[edge.node[2]])
                end
                assemble_res(edge, system, easm_res)
            end
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
    return integrate(system, tf, U, U, Inf; kwargs...)
end

############################################################################
"""
$(SIGNATURES)

Steady state part of test function integral.
"""
function integrate_stdy(system::AbstractSystem, tf::Vector{Tv}, U::AbstractArray{Tu, 2}; data = system.physics.data) where {Tu, Tv}
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

    src_eval = ResEvaluator(physics, data, :source, UK, node, nspecies)
    rea_eval = ResEvaluator(physics, data, :reaction, UK, node, nspecies)
    erea_eval = ResEvaluator(physics, data, :edgereaction, UK, node, nspecies)
    flux_eval = ResEvaluator(physics, data, :flux, UKL, edge, nspecies)

    for item in nodebatch(system.assembly_data)
        for inode in noderange(system.assembly_data, item)
            _fill!(node, system.assembly_data, inode, item)
            @views UK .= U[:, node.index]

            evaluate!(rea_eval, UK)
            rea = res(rea_eval)
            evaluate!(src_eval)
            src = res(src_eval)

            function asm_res(idof, ispec)
                return integral[ispec] += node.fac * (rea[ispec] - src[ispec]) * tf[node.index]
            end
            assemble_res(node, system, asm_res)
        end
    end

    for item in edgebatch(system.assembly_data)
        for iedge in edgerange(system.assembly_data, item)
            _fill!(edge, system.assembly_data, iedge, item)
            @views UKL[1:nspecies] .= U[:, edge.node[1]]
            @views UKL[(nspecies + 1):(2 * nspecies)] .= U[:, edge.node[2]]
            evaluate!(flux_eval, UKL)
            flux = res(flux_eval)

            function asm_res(idofK, idofL, ispec)
                return integral[ispec] += edge.fac * flux[ispec] * (tf[edge.node[1]] - tf[edge.node[2]])
            end
            assemble_res(edge, system, asm_res)

            if isnontrivial(erea_eval)
                evaluate!(erea_eval, UKL)
                erea = res(erea_eval)

                function easm_res(idofK, idofL, ispec)
                    return integral[ispec] += edge.fac * erea[ispec] * (tf[edge.node[1]] + tf[edge.node[2]])
                end
                assemble_res(edge, system, easm_res)
            end
        end
    end

    return integral
end

############################################################################
"""
$(SIGNATURES)

Calculate transient part of test function integral.
"""
function integrate_tran(system::AbstractSystem, tf::Vector{Tv}, U::AbstractArray{Tu, 2}; data = system.physics.data) where {Tu, Tv}
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
    stor_eval = ResEvaluator(physics, data, :storage, UK, node, nspecies)

    for item in nodebatch(system.assembly_data)
        for inode in noderange(system.assembly_data, item)
            _fill!(node, system.assembly_data, inode, item)
            @views UK .= U[:, node.index]
            evaluate!(stor_eval, UK)
            stor = res(stor_eval)
            asm_res(idof, ispec) = integral[ispec] += node.fac * stor[ispec] * tf[node.index]
            assemble_res(node, system, asm_res)
        end
    end

    return integral
end
