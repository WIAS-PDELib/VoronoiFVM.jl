################################################
"""
$(TYPEDEF)

Data structure containing DenseSystem used to calculate
test functions for boundary flux calculations.

Type parameters:
- `Tu`: value type of test functions
- `Tv`: Default value type of system

Fields:
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
    TestFunctionFactory(system; control= SolverControl())

Constructor for TestFunctionFactory from system.
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
    testfunction(factory::TestFunctionFactory, bc0, bc1)

Create testfunction which has Dirichlet zero boundary conditions  for boundary
regions listed in `bc0` and Dirichlet one boundary conditions  for boundary
regions listed in `bc1`. 
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

    eval_and_assemble(
        factory.state.system, u, u, f,
        factory.state.matrix, factory.state.generic_matrix, factory.state.dudp,
        Inf, Inf, 0.0, nothing, zeros(0)
    )

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
     integrate(system, T, U::AbstractMatrix)

Calculate test function integral for a steady state solution ``u``  ``∫_Γ T \\vec j(u) ⋅ \\vec n dω ≈ I_j(T,u)-I_r(T,u)``,
see the definition of [test function integral contributions](@ref discrete_appr).

The result is a `nspec` vector giving one value of the integral for each species. 
"""
function integrate(
        system::AbstractSystem,
        tf::Vector{Tv},
        U::AbstractMatrix{Tu};
        kwargs...
    ) where {Tu, Tv}
    return integrate(system, tf, U, U, Inf; kwargs...)
end

############################################################################

"""
    integrate(system,T, U::AbstractTransientSolution; rate=true, params, data)

Calculate test function integrals for the transient solution  ``∫_Γ T \\vec j(u^n) ⋅ \\vec n dω ≈ I_j(T,u^n)-I_r(T,u^n)-I_{s_t}(T,u^{n-1}, u^n)``
for each time interval ``(t^{n-1}, t^n)``,see the definition of [test function integral contributions](@ref discrete_appr).

Keyword arguments:
- `params`: vector of parameters used to calculate `U`. Default: `[]`.
- `data`: user data   used to calculate `U`. Default: `data(system)`
- `rate`: If `rate=true` (default), calculate the flow rate (per second) 
   through the corresponding boundary. Otherwise, calculate the absolute 
   amount per time inteval. 

The result is a `nspec x (nsteps-1)` DiffEqArray.
"""
function integrate(
        sys::AbstractSystem,
        tf::Vector,
        U::AbstractTransientSolution;
        rate = true,
        kwargs...
    )
    nsteps = length(U.t) - 1
    integral = [
        VoronoiFVM.integrate(
                sys,
                tf,
                U.u[istep + 1],
                U.u[istep],
                U.t[istep + 1] - U.t[istep];
                kwargs...
            ) / (rate ? U.t[istep + 1] - U.t[istep] : 1)
            for istep in 1:nsteps
    ]
    return DiffEqArray(integral, U.t[2:end])
end


"""
  integrate(system, T, U::AbstractMatrix, Uold::AbstractMatrix, tstep; kwargs...)

Calculate test function integrals for the implicit Euler time step results  ``∫_Γ T \\vec j(u) ⋅ \\vec n dω ≈ I_j(T,u)-I_r(T,u)-I_{s_t}(T,u_{old}, u)``,
see the definition of [test function integral contributions](@ref discrete_appr).

Keyword arguments:
- `params`: vector of parameters used to calculate `U`. Default: `[]`.
- `data`: user data   used to calculate `U`. Default: `data(system)`

The result is a `nspec` vector giving one value of the integral for each species. 
"""
function integrate(
        system::AbstractSystem,
        tf,
        U::AbstractMatrix{Tv},
        Uold::AbstractMatrix{Tv},
        tstep;
        params = Tv[],
        data = system.physics.data
    ) where {Tv}

    integral1 = integrate_nodebatch(system, tf, U, Uold, tstep; params, data)
    integral2 = integrate_edgebatch(system, tf, U, Uold, tstep; params, data)

    return integral1 .+ integral2
end


"""
     integrate_flux_time_derivative(system, T, U::AbstractMatrix, Uold::AbstractMatrix, tstep; kwargs...)

Calculate test function integral for the time time derivative of the current density. More precisely, this method
computes the current ``\\int_{\\Omega} \\nabla T \\cdot \\partial_t \\vec j dx`` using a test function approach.
This method can be used for to calculate the displacement current for  Poisson-Nernst Planck or the van Roosbroeck models.
See [Farrell et al.,  Numerical Methods for Drift-DIffusion Models, WIAS Preprint No. 2263, (2016)](http://dx.doi.org/10.20347/WIAS.PREPRINT.2263).
"""
function integrate_flux_time_derivative(
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
    edge = Edge(system, 0.0, 1.0, params)
    bedge = Edge(system, 0.0, 1.0, params)

    UKL = Array{Tv, 1}(undef, 2 * nspecies + nparams)
    UK = Array{Tv, 1}(undef, nspecies + nparams)
    UKLold = Array{Tv, 1}(undef, 2 * nspecies + nparams)
    UKold = Array{Tv, 1}(undef, nspecies + nparams)

    if nparams > 0
        UKL[(2 * nspecies + 1):end] .= params
        UKLold[(2 * nspecies + 1):end] .= params
    end

    erea_eval = ResEvaluator(physics, data, :edgereaction, UK, edge, nspecies + nparams)
    ereaold_eval = ResEvaluator(physics, data, :edgereaction, UKold, edge, nspecies + nparams)
    flux_eval = ResEvaluator(physics, data, :flux, UKL, edge, nspecies + nparams)
    fluxold_eval = ResEvaluator(physics, data, :flux, UKLold, edge, nspecies + nparams)

    for item in edgebatch(system.assembly_data)
        for iedge in edgerange(system.assembly_data, item)
            _fill!(edge, system.assembly_data, iedge, item)
            @views UKL[1:nspecies] .= U[:, edge.node[1]]
            @views UKL[(nspecies + 1):(2 * nspecies)] .= U[:, edge.node[2]]

            @views UKLold[1:nspecies] .= Uold[:, edge.node[1]]
            @views UKLold[(nspecies + 1):(2 * nspecies)] .= Uold[:, edge.node[2]]

            evaluate!(flux_eval, UKL)
            flux = res(flux_eval)

            evaluate!(fluxold_eval, UKLold)
            fluxold = res(fluxold_eval)

            function asm_res(idofK, idofL, ispec)
                return integral[ispec] += edge.fac * (flux[ispec] - fluxold[ispec]) * tstepinv * (tf[edge.node[1]] - tf[edge.node[2]])
            end
            assemble_res(edge, system, asm_res)

            if isnontrivial(erea_eval)
                evaluate!(erea_eval, UKL)
                erea = res(erea_eval)

                evaluate!(ereaold_eval, UKLold)
                ereaold = res(ereaold_eval)

                function easm_res(idofK, idofL, ispec)
                    return integral[ispec] += edge.fac * (erea[ispec] - ereaold[ispec]) * tstepinv * (tf[edge.node[1]] + tf[edge.node[2]])
                end
                assemble_res(edge, system, easm_res)
            end

        end
    end

    return integral

end


"""
    integrate_nodebatch(system, T, U, Uold, tstep; kwargs...)

Calculate  node contribution ``I_{s_t}(T, u_{old}, u) + I_r(T,u)``  to test function
integral   ``∫_Γ T \\vec j(u) ⋅ \\vec n dω`` for a given timestep solution.
See [the time step `integrate` method](@ref VoronoiFVM.integrate(::VoronoiFVM.AbstractSystem,::Any, ::AbstractMatrix{Tv},::AbstractMatrix{Tv}, ::Any; kwargs...) where {Tv}).
"""
function integrate_nodebatch(
        system::AbstractSystem,
        tf,
        U::AbstractMatrix{Tv},
        Uold::AbstractMatrix{Tv},
        tstep;
        params = Tv[],
        data = system.physics.data
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

    UK = Array{Tv, 1}(undef, nspecies + nparams)
    UKold = Array{Tv, 1}(undef, nspecies + nparams)

    if nparams > 0
        UK[(nspecies + 1):end] .= params
        UKold[(nspecies + 1):end] .= params
    end

    src_eval = ResEvaluator(physics, data, :source, UK, node, nspecies + nparams)
    rea_eval = ResEvaluator(physics, data, :reaction, UK, node, nspecies + nparams)
    stor_eval = ResEvaluator(physics, data, :storage, UK, node, nspecies + nparams)
    storold_eval = ResEvaluator(physics, data, :storage, UKold, node, nspecies + nparams)

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
    return integral
end


"""
    integrate_edgebatch(system, T, U, Uold, tstep; kwargs...)

Calculate  edge (flux) contribution ``I_j(T,u)``  to test function
integral   ``∫_Γ T \\vec j(u) ⋅ \\vec n dω`` for a given timestep solution.
See [the time step `integrate` method](@ref VoronoiFVM.integrate(::VoronoiFVM.AbstractSystem,::Any, ::AbstractMatrix{Tv},::AbstractMatrix{Tv}, ::Any; kwargs...) where {Tv}).
"""
function integrate_edgebatch(
        system::AbstractSystem,
        tf,
        U::AbstractMatrix{Tv},
        Uold::AbstractMatrix{Tv},
        tstep;
        params = Tv[],
        data = system.physics.data
    ) where {Tv}

    grid = system.grid
    nspecies = num_species(system)
    integral = zeros(Tv, nspecies)
    tstepinv = 1.0 / tstep
    nparams = system.num_parameters
    @assert nparams == length(params)

    # !!! params etc
    physics = system.physics
    edge = Edge(system, 0.0, 1.0, params)
    bedge = Edge(system, 0.0, 1.0, params)

    UK = Array{Tv, 1}(undef, nspecies + nparams)
    UKL = Array{Tv, 1}(undef, 2 * nspecies + nparams)

    if nparams > 0
        UK[(nspecies + 1):end] .= params
        UKL[(2 * nspecies + 1):end] .= params
    end

    erea_eval = ResEvaluator(physics, data, :edgereaction, UK, edge, nspecies + nparams)
    flux_eval = ResEvaluator(physics, data, :flux, UKL, edge, nspecies + nparams)

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
    integrate_stdy(system, T, U; kwargs...)

Calculate  steady state contribution ``I_j(T,u) - I_r(T,u)``  to test function
integral   ``∫_Γ T \\vec j(u) ⋅ \\vec n dω`` for a given timestep solution.
See [the time step `integrate` method](@ref VoronoiFVM.integrate(::VoronoiFVM.AbstractSystem,::Any, ::AbstractMatrix{Tv},::AbstractMatrix{Tv}, ::Any; kwargs...) where {Tv}).

Used for impedance calculations.
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
    integrate_tran(system, T, U; kwargs...)

Calculate  storage term contribution to test function integral. Used for impedance calculations.
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
