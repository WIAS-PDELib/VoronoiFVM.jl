"""
    $(SIGNATURES)

Evaluate functiaon and Jacobian at u if they have not been evaluated before at u.
See https://github.com/SciML/DifferentialEquations.jl/issues/521 for discussion of another way to do this.

Evaluation needs to be done with tstep=Inf, like for a stationary problem, because the ODE solvers
handle the time derivative. 
"""
function _eval_res_jac!(state, u, t)
    uhash = hash(u)
    if uhash != state.uhash
        ur = reshape(u, state.system)
        eval_and_assemble(state.system, ur, ur, state.residual, state.matrix, state.dudp, value(t), Inf, 0.0, state.system.physics.data, state.params)
        state.uhash = uhash
        state.history.nd += 1
    end
    return nothing
end

"""
$(SIGNATURES)

Interpret the  discrete problem as an ODE/DAE problem. Provide the 
rhs function for [`SciMLBase.ODEFunction`](@ref).
"""
function eval_rhs!(du, u, state, t)
    _eval_res_jac!(state, u, t)
    du .= -dofs(state.residual)
    state.history.nf += 1
    return nothing
end

"""
$(SIGNATURES)

Interpret the  discrete problem as an ODE/DAE problem. Provide the 
jacobi matrix calculation function for [`SciMLBase.ODEFunction`](@ref)
"""
function eval_jacobian!(J, u, state, t)
    _eval_res_jac!(state, u, t)
    # Need to implement broadcast for ExtendableSparse.
    J .= -state.matrix.cscmatrix
    state.history.njac += 1
    return nothing
end

"""
$(SIGNATURES)

Calculate the mass matrix for use with [`SciMLBase.ODEFunction`](@ref).
Return a Diagonal matrix if it occurs to be diagonal, otherwise return a SparseMatrixCSC.
"""
function mass_matrix(state::SystemState{Tv, TMatrix, TSolArray, TData}) where {Tv, TMatrix, TSolArray, TData}
    physics = state.system.physics
    data = physics.data
    node = Node(state.system)
    bnode = BNode(state.system)
    nspecies = num_species(state.system)
    ndof = num_dof(state.system)
    data = state.system.physics.data

    stor_eval = ResJacEvaluator(physics, data, :storage, zeros(Tv, nspecies), node, nspecies)
    bstor_eval = ResJacEvaluator(physics, data, :bstorage, zeros(Tv, nspecies), bnode, nspecies)

    U = unknowns(state.system; inival = 0)
    M = similar(state.matrix)

    asm_res(idof, ispec) = nothing
    asm_param(idof, ispec, iparam) = nothing

    for item in nodebatch(state.system.assembly_data)
        for inode in noderange(state.system.assembly_data, item)
            _fill!(node, state.system.assembly_data, inode, item)
            @views evaluate!(stor_eval, U[:, node.index])
            jac_stor = jac(stor_eval)
            asm_jac(idof, jdof, ispec, jspec) = _addnz(M, idof, jdof, jac_stor[ispec, jspec], node.fac)
            assemble_res_jac(node, state.system, asm_res, asm_jac, asm_param)
        end
    end

    if isnontrivial(bstor_eval)
        for item in nodebatch(state.system.boundary_assembly_data)
            for ibnode in noderange(state.system.boundary_assembly_data, item)
                _fill!(bnode, state.system.boundary_assembly_data, ibnode, item)
                K = bnode.index
                @views evaluate!(bstor_eval, U[:, K])
                jac_bstor = jac(bstor_eval)
                asm_jac(idof, jdof, ispec, jspec) = _addnz(M, idof, jdof, jac_bstor[ispec, jspec], bnode.fac)
                assemble_res_jac(bnode, state.system, asm_res, asm_jac, asm_param)
            end
        end
    end
    Mcsc = SparseMatrixCSC(M)
    return isdiag(Mcsc) ? Diagonal([Mcsc[i, i] for i in 1:ndof]) : Mcsc
end

"""
$(SIGNATURES)
Prepare system for use with VoronoiFVMDiffEq.

- `jacval`: value at which to evaluate jacobian to obtatin prototype
- `tjac`: time moment for jacobian
 
Returns a prototype for the jacobian.
"""
function prepare_diffeq!(state, jacval, tjac)
    state.history = DiffEqHistory()
    _complete!(state.system)
    _eval_res_jac!(state, jacval, tjac)
    flush!(state.matrix)
    return state.matrix.cscmatrix
end

###################################################################################################
# API

"""
     ODEFunction(state; jacval = unknowns(sys, 0), tjac = 0 )
    
Create an [ODEFunction](https://diffeq.sciml.ai/stable/basics/overview/#Defining-Problems).
For more documentation, see [`SciMLBase.ODEFunction(state::VoronoiFVM.SystemState; kwargs...)`](@ref)
`jacval` and `tjac` are passed to [`prepare_diffeq!`](@ref) and used there to calculate the Jacobian prototype.

Defined in VoronoiFVM.jl.
"""
function SciMLBase.ODEFunction(state::VoronoiFVM.SystemState; jacval = unknowns(sys, 0), tjac = 0)
    return SciMLBase.ODEFunction(
        eval_rhs!;
        jac = eval_jacobian!,
        jac_prototype = prepare_diffeq!(state, dofs(jacval), tjac),
        mass_matrix = mass_matrix(state)
    )
end

"""
     ODEFunction(system; jacval=unknowns(system,inival=0),tjac=0)
    
Create an [ODEFunction](https://diffeq.sciml.ai/stable/basics/overview/#Defining-Problems)
in [mass matrix form](https://diffeq.sciml.ai/stable/solvers/dae_solve/#OrdinaryDiffEq.jl-(Mass-Matrix))
to be handled by ODE solvers from [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

Parameters:
- `system`: A [`VoronoiFVM.System`](https://WIAS-PDELib.github.io/VoronoiFVM.jl/stable/system/#VoronoiFVM.System-Tuple{ExtendableGrid})
- `jacval` (optional): Initial value. Default is a zero vector. Consider to  pass a stationary solution at time `tjac`.
- `tjac` (optional): tjac, Default: 0

The `jacval` and `tjac` are passed  for a first evaluation of the Jacobian, allowing to detect
the sparsity pattern which is passed to the solver.

Defined in VoronoiFVM.jl.
"""
SciMLBase.ODEFunction(sys::VoronoiFVM.System; kwargs...) = SciMLBase.ODEFunction(SystemState(sys); kwargs...)

"""
    ODEProblem(state,inival,tspan,callback=SciMLBase.CallbackSet())

Create an [ODEProblem](https://diffeq.sciml.ai/stable/basics/overview/#Defining-Problems) from  a
system state. See [`SciMLBase.ODEProblem(sys::VoronoiFVM.System, inival, tspan;kwargs...)`](@ref)
for more documentation.

Defined in VoronoiFVM.jl.
"""
function SciMLBase.ODEProblem(
        state::VoronoiFVM.SystemState, inival, tspan;
        params = state.params, callback = SciMLBase.CallbackSet()
    )
    state.params .= params
    odefunction = SciMLBase.ODEFunction(state; jacval = dofs(inival), tjac = tspan[1])
    return SciMLBase.ODEProblem(odefunction, dofs(inival), tspan, state, callback)
end

"""
    ODEProblem(system,inival,tspan,callback=SciMLBase.CallbackSet())
    
Create an [ODEProblem](https://diffeq.sciml.ai/stable/basics/overview/#Defining-Problems)
in [mass matrix form](https://diffeq.sciml.ai/stable/solvers/dae_solve/#OrdinaryDiffEq.jl-(Mass-Matrix))
which can  be handled by ODE solvers from [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

Parameters:
- `system`: A [`VoronoiFVM.System`](https://WIAS-PDELib.github.io/VoronoiFVM.jl/stable/system/#VoronoiFVM.System-Tuple{ExtendableGrid})
- `inival`: Initial value. Consider to  pass a stationary solution at `tspan[1]`.
- `tspan`: Time interval 
- `callback` : (optional) [callback](https://diffeq.sciml.ai/stable/features/callback_functions/#Using-Callbacks) for ODE solver 

The method returns an [ODEProblem](https://diffeq.sciml.ai/stable/basics/overview/#Defining-Problems) which can be solved
by [solve()](https://diffeq.sciml.ai/stable/basics/common_solver_opts/).

Defined in VoronoiFVM.jl.
"""
function SciMLBase.ODEProblem(
        sys::VoronoiFVM.System, inival, tspan;
        params = zeros(sys.num_parameters), kwargs...
    )
    return SciMLBase.ODEProblem(SystemState(sys), inival, tspan; params, kwargs...)
end

"""
    reshape(ode_solution, system; times=nothing, state=nothing)
Create a [`TransientSolution`](@ref) from the output of the ode solver which
reflects the species structure of the system ignored by the ODE solver.
Howvever the interpolation behind `reshaped_sol(t)` will be linear and ignores the possibility
of higher order interpolations with `ode_sol`.

If `times` is specified, the (possibly higher order) interpolated solution at the given moments of time will be returned.

Defined in VoronoiFVM.jl.
"""
function Base.reshape(sol::AbstractDiffEqArray, sys::VoronoiFVM.AbstractSystem; times = nothing, state = nothing)
    if isnothing(times)
        tsol = TransientSolution([reshape(sol.u[i], sys) for i in 1:length(sol.u)], sol.t)
    else
        isol = sol(times)
        tsol = TransientSolution([reshape(isol[i], sys) for t in 1:length(times)], times)
    end
    if !isnothing(state)
        tsol.history = state.history
    end
    return tsol
end
