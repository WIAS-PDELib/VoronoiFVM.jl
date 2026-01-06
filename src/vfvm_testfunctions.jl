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
    testfunction(factory::TestFunctionFactory, B_0, B_1)

Create a test function ``T`` defined in the domain ``Ω`` given by the system
behind the factory object. Assume  that the  boundary ``Γ=∂Ω=⋃_{i∈ B}Γ_i`` is subdivided into non-overlapping parts
where the set of boundary regions ``B`` is subdivided into ``B=B_N∪ B_0 ∪ B_1``.
Create ``T`` by solving
```math
   -Δ T =0 \\quad \\text{in}\\; Ω
```
such that
- ``∇  T ⋅  \\vec n|_{Γ_i}=  0`` for ``i∈ B_N``
- ``T|_{Γ_i} = 0`` for  ``i∈ B_0``,
- ``T|_{Γ_i} = 1`` for  ``i∈ B_1``.

Returns a vector representing T.
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

Calculate test function integral for a steady state solution ``u``  
``∫_Γ T \\vec j(u) ⋅ \\vec n dω ≈ I_{flux}(j, T,u)-I_{func}(r,T,u) + I_{src}(f,T) ``,
see the definition of [test function integral contributions](@ref discrete_appr).

The result is a `nspec` vector giving one value of the integral for each species. 
"""
function integrate(
        system::AbstractSystem,
        T::Vector{Tv},
        U::AbstractMatrix{Tu};
        kwargs...
    ) where {Tu, Tv}
    I_react = integrate_TxFunc(system, T, system.physics.reaction, U; kwargs...)
    I_src = integrate_TxSrc(system, T, system.physics.source; kwargs...)
    I_flux = integrate_∇TxFlux(system, T, system.physics.flux, U; kwargs...)

    I = I_flux + I_react - I_src

    if system.physics.edgereaction != nofunc
        I += integrate_TxEdgefunc(system, T, system.physics.edgereaction, U; kwargs...)
    end

    return I
end

############################################################################

"""
    integrate(system,T, U::AbstractTransientSolution; rate=true, params, data)

Calculate test function integrals for the transient solution  
``∫_Γ T \\vec j(u^n) ⋅ \\vec n dω ≈ I_{flux}(j,T,u^n)-I_{func}(r,T,u^n)-\\frac{I_{func}(s,T,u^n)-I_{func}(s,T,u^{n-1})}{t^n-t^{n-1}} + I_{src}(f,T)``
for each time interval ``(t^{n-1}, t^n)``, see the definition of [test function integral contributions](@ref discrete_appr).

Keyword arguments:
- `params`: vector of parameters used to calculate `U`. Default: `[]`.
- `data`: user data   used to calculate `U`. Default: `data(system)`
- `rate`: If `rate=true` (default), calculate the flow rate (per second) 
   through the corresponding boundary. Otherwise, calculate the absolute 
   amount per time interval. 

The result is a `nspec x (nsteps-1)` DiffEqArray.
"""
function integrate(
        sys::AbstractSystem,
        T::Vector,
        U::AbstractTransientSolution;
        rate = true,
        kwargs...
    )
    nsteps = length(U.t) - 1
    integral = [
        integrate(
                sys,
                T,
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
     integrate(system, T, U::AbstractMatrix, Uold::AbstractMatrix, Δt; kwargs...)

Calculate test function integrals for the implicit Euler time step results ``u^n≡`` `U` and
``u^{n-1}≡`` `Uold` as 
``∫_Γ T \\vec j(u^n) ⋅ \\vec n dω ≈ I_{flux}(j,T,u^n)-I_{func}(r,T,u^n)-\\frac{I_{func}(s,T,u^n)-I_{func}(s,T,u^{n-1})}{Δt} + I_{src}(f,T)``,
see the definition of [test function integral contributions](@ref discrete_appr).

Keyword arguments:
- `params`: vector of parameters used to calculate `U`. Default: `[]`.
- `data`: user data   used to calculate `U`. Default: `data(system)`

The result is a `nspec` vector giving one value of the integral for each species. 
"""
function integrate(
        system::AbstractSystem,
        T,
        U::AbstractMatrix{Tv},
        Uold::AbstractMatrix{Tv},
        tstep;
        kwargs...
    ) where {Tv}

    I_react = integrate_TxFunc(system, T, system.physics.reaction, U; kwargs...)
    I_stor = integrate_TxFunc(system, T, system.physics.storage, U; kwargs...)
    I_oldstor = integrate_TxFunc(system, T, system.physics.storage, Uold; kwargs...)
    I_src = integrate_TxSrc(system, T, system.physics.source; kwargs...)
    I_flux = integrate_∇TxFlux(system, T, system.physics.flux, U; kwargs...)
    I = I_flux + I_react - I_src + (I_stor - I_oldstor) / tstep
    if system.physics.edgereaction != nofunc
        I += integrate_TxEdgefunc(system, T, system.physics.edgereaction, U; kwargs...)
    end
    return I
end

"""
    integrate_TxFunc(system, T, f!, U; kwargs...)

Calculate ``I_{func}(f!,T, U)=∫_Ω T⋅ f!(U) dω`` for a test function `T` and unknown vector `U`. 
The function `f!` shall have the same signature as a storage or reaction function.
See the definition of [test function integral contributions](@ref discrete_appr).
"""
function integrate_TxFunc(
        system::AbstractSystem, T, func!::Tfunc,
        U::AbstractMatrix{Tv};
        params = Tv[], data = system.physics.data
    ) where {Tfunc, Tv}
    grid = system.grid
    nspecies = num_species(system)
    integral = zeros(Tv, nspecies)
    nparams = system.num_parameters
    @assert nparams == length(params)

    physics = system.physics
    node = Node(system, 0.0, 1.0)

    UK = Array{Tv, 1}(undef, nspecies + nparams)
    YK = Array{Tv, 1}(undef, nspecies)

    if nparams > 0
        UK[(nspecies + 1):end] .= params
    end
    VK = unknowns(node, UK)

    for item in nodebatch(system.assembly_data)
        for inode in noderange(system.assembly_data, item)
            _fill!(node, system.assembly_data, inode, item)
            for ispec in 1:nspecies
                UK[ispec] = U[ispec, node.index]
            end
            YK .= zero(Tv)
            func!(YK, VK, node, data)
            for ispec in 1:nspecies
                integral[ispec] += node.fac * YK[ispec] * T[node.index]
            end
        end
    end
    return integral
end

"""
    integrate_TxSrc(system, T, f!; kwargs...)

Calculate ``I_{src}(f!,T)= ∫_Ω T⋅ f! dω`` for a test function `T. 
The function `f!` shall have the same signature as a storage or reaction function.
"""
function integrate_TxSrc(
        system::AbstractSystem{Tv}, T, src!::Tsrc;
        params = Tv[], data = system.physics.data
    ) where {Tsrc, Tv}
    grid = system.grid
    nspecies = num_species(system)
    integral = zeros(Tv, nspecies)
    nparams = system.num_parameters
    @assert nparams == length(params)

    physics = system.physics
    node = Node(system, 0.0, 1.0)

    YK = Array{Tv, 1}(undef, nspecies)


    for item in nodebatch(system.assembly_data)
        for inode in noderange(system.assembly_data, item)
            _fill!(node, system.assembly_data, inode, item)
            YK .= zero(Tv)
            src!(YK, node, data)
            for ispec in 1:nspecies
                integral[ispec] += node.fac * YK[ispec] * T[node.index]
            end
        end
    end
    return integral
end

"""
      integrate_∇TxFlux(system, T, f!, U; kwargs...)

Calculate ``I_{flux}(f!,T, U)=∫_Ω ∇T⋅ f!(U) dω`` for a test function `T` and unknown vector `U`. 
The function `f!` shall have the same signature as a flux function.
"""
function integrate_∇TxFlux(
        system::AbstractSystem,
        T,
        flux!::Tflux,
        U::AbstractMatrix{Tv};
        params = Tv[],
        data = system.physics.data
    ) where {Tflux, Tv}
    grid = system.grid
    nspecies = num_species(system)
    integral = zeros(Tv, nspecies)
    nparams = system.num_parameters
    @assert nparams == length(params)

    edge = Edge(system, 0.0, 1.0)
    UKL = Array{Tv, 1}(undef, 2 * nspecies + nparams)
    YK = Array{Tv, 1}(undef, nspecies)

    if nparams > 0
        UKL[(2 * nspecies + 1):end] .= params
    end

    VKL = unknowns(edge, UKL)

    for item in edgebatch(system.assembly_data)
        for iedge in edgerange(system.assembly_data, item)
            _fill!(edge, system.assembly_data, iedge, item)
            @views UKL[1:nspecies] .= U[:, edge.node[1]]
            @views UKL[(nspecies + 1):(2 * nspecies)] .= U[:, edge.node[2]]
            YK .= zero(Tv)
            flux!(YK, VKL, edge, data)
            for ispec in 1:nspecies
                integral[ispec] += edge.fac * YK[ispec] * (T[edge.node[1]] - T[edge.node[2]])
            end
        end
    end
    return integral
end

"""
      integrate_∇TxFlux(system, T, U; kwargs...)

Calculate ``I_{flux}(f!,T, U)=∫_Ω ∇T⋅ f!(U) dω`` for a test function `T` and unknown vector `U`,
where `flux!` is the system flux.
"""
integrate_∇TxFlux(system, T, U; kwargs...) = integrate_∇TxFlux(system, T, system.physics.flux, U; kwargs...)

"""
      integrate_TxEdgeFunc(system, T, f!, U; kwargs...)

Calculate ``∫_Ω T⋅ f!(U) dω`` for a test function `T` and unknown vector `U`. 
The function `f!` shall have the same signature as a flux function, but is assumed
to describe a reaction term given on edges.
"""
function integrate_TxEdgefunc(
        system::AbstractSystem, T, func!::Tfunc, U::AbstractMatrix{Tv};
        params = Tv[], data = system.physics.data
    ) where {Tfunc, Tv}
    grid = system.grid
    nspecies = num_species(system)
    integral = zeros(Tv, nspecies)
    nparams = system.num_parameters
    @assert nparams == length(params)

    edge = Edge(system, 0.0, 1.0)
    UKL = Array{Tv, 1}(undef, 2 * nspecies + nparams)
    YK = Array{Tv, 1}(undef, nspecies)

    VKL = unknowns(edge, UKL)
    if nparams > 0
        UKL[(2 * nspecies + 1):end] .= params
    end

    for item in edgebatch(system.assembly_data)
        for iedge in edgerange(system.assembly_data, item)
            _fill!(edge, system.assembly_data, iedge, item)
            @views UKL[1:nspecies] .= U[:, edge.node[1]]
            @views UKL[(nspecies + 1):(2 * nspecies)] .= U[:, edge.node[2]]
            YK .= zero(Tv)
            func!(YK, VKL, edge, data)
            for ispec in 1:nspecies
                integral[ispec] += edge.fac * YK[ispec] * (T[edge.node[1]] + T[edge.node[2]])
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
function integrate_stdy(system::AbstractSystem, T::Vector{Tv}, U::AbstractArray{Tu, 2}; kwargs...) where {Tu, Tv}

    I_react = integrate_TxFunc(system, T, system.physics.reaction, U; kwargs...)
    I_src = integrate_TxSrc(system, T, system.physics.source; kwargs...)
    I_flux = integrate_∇TxFlux(system, T, system.physics.flux, U; kwargs...)

    I = I_flux + I_react - I_src

    if system.physics.edgereaction != nofunc
        I += integrate_TxEdgefunc(system, T, system.physics.edgereaction, U; kwargs...)
    end
    return I
end

############################################################################
"""
    integrate_tran(system, T, U; kwargs...)

Calculate  storage term contribution to test function integral. Used for impedance calculations.
"""
function integrate_tran(system::AbstractSystem, T::Vector{Tv}, U::AbstractArray{Tu, 2}; kwargs...) where {Tu, Tv}
    return I_stor = integrate_TxFunc(system, T, system.physics.storage, U; kwargs...)
end
