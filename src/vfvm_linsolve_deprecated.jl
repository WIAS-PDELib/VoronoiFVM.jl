################################################
# deprecated "homegrown API" to work around the need to specify preconditioners
# instead of preconditioner constructors to LinearSolve.jl pre 2.36.
# TODO: remove with v3.0

abstract type AbstractStrategy end

"""
    VoronoiFVM.LinearSolverStrategy

An linear solver strategy provides the possibility to construct [`SolverControl`](@ref) objects as follows:
```
    SolverControl(strategy,sys;kwargs...)
```,
e.g.
```
    SolverControl(GMRESIteration(UMFPackFactorization(), EquationBlock()),sys; kwargs...)
```

A linear solver strategy combines a Krylov method  with a preconditioner
which by default is calculated from the linearization of the initial value of the
Newton iteration. For coupled systems, a blocking strategy can be chosen. The [`EquationBlock`](@ref) strategy
calculates preconditioners or LU factorization  separately for each species equation and combines them
to a block Jacobi preconditioner.  The [`PointBlock`](@ref) strategy treats the linear system as consisting
of `nspecies x nspecies` blocks. 

Which is the best strategy, depends on the space dimension. The following is a rule of thumb for choosing strategies
- For 1D problems use direct solvers
- For 2D stationary problems, use direct solvers, for transient problems consider iterative solvers which 
  can take advantage of the diagonal dominance of the implicit timestep problem
- For 3D problems avoid direct solvers


Currently available strategies are:
- [`DirectSolver`](@ref)
- [`CGIteration`](@ref)
- [`BICGstabIteration`](@ref)
- [`GMRESIteration`](@ref)

Notable LU Factorizations/direct solvers are:
- [`UMFPACKFactorization`](https://docs.sciml.ai/LinearSolve/stable/solvers/solvers/#SuiteSparse.jl)  (`using LinearSolve`)
- [`KLUFactorization`](https://docs.sciml.ai/LinearSolve/stable/solvers/solvers/#SuiteSparse.jl) (`using LinearSolve`)
- [`SparspakFactorization`](https://docs.sciml.ai/LinearSolve/stable/solvers/solvers/#Sparspak.jl)  (`using LinearSolve`), [`SparspakLU`](https://WIAS-PDELib.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.SparspakLU) (`using ExtendableSparse,Sparspak`)
- [`MKLPardisoLU`](https://WIAS-PDELib.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.MKLPardisoLU) (`using ExtendableSparse, Pardiso`), openmp parallel
- [`AMGSolver`](https://j-fu.github.io/AMGCLWrap.jl/stable/solvers/#AMGCLWrap.AMGSolver) (`using AMGCLWrap`), openmp parallel
- [`RLXSolver`](https://j-fu.github.io/AMGCLWrap.jl/stable/solvers/#AMGCLWrap.RLXSolver) (`using AMGCLWrap`), openmp parallel


Notable incomplete factorizations/preconditioners
- Incomplete LU factorizations written in Julia:
    - [`ILUZeroPreconditioner`](https://WIAS-PDELib.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.ILUZeroPreconditioner)
    - [`ILUTPrecondidtioner`](https://WIAS-PDELib.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.ILUTPreconditioner) (`using ExtendableSparse, IncompleteLU`)
- Algebraic multigrid written in Julia: (`using ExtendableSparse, AlgebraicMultigrid`)
    - [`RS_AMGPreconditioner`](https://WIAS-PDELib.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.RS_AMGPreconditioner)
    - [`SA_AMGPreconditioner`](https://WIAS-PDELib.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.SA_AMGPreconditioner)
- AMGCL based preconditioners (`using ExtendableSparse, AMGCLWrap`), openmp parallel
    - [`AMGCL_AMGPreconditioner`](https://WIAS-PDELib.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.AMGCL_AMGPreconditioner)
    - [`AMGCL_RLXPreconditioner`](https://WIAS-PDELib.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.AMGCL_RLXPreconditioner)

Blocking strategies are:
- [`NoBlock`](@ref)
- [`EquationBlock`](@ref)
- [`PointBlock`](@ref)

"""
abstract type LinearSolverStrategy <: AbstractStrategy end

"""
    VoronoiFVM.BlockStrategy

Abstract supertype for various block preconditioning strategies.
"""
abstract type BlockStrategy <: AbstractStrategy end

"""
    NoBlock()

No blocking.
"""
struct NoBlock <: BlockStrategy end

"""
    EquationBlock()

Equation-wise blocking. Can be combined with any preconditioner resp. factorization including direct solvers.
"""
struct EquationBlock <: BlockStrategy end

"""
    PointBlock()

Point-wise blocking. Currently only together with ILUZeroFactorization.
This requires a system with `unknown_storage=:dense`.
"""
struct PointBlock <: BlockStrategy end

"""
    DirectSolver(;factorization=UMFPACKFactorization())

LU Factorization solver.
"""
Base.@kwdef struct DirectSolver <: LinearSolverStrategy
    factorization::FactorizationStrategy = UMFPACKFactorization()
    blocking::NoBlock = NoBlock() # prevent ambiguity in constructor definition
end

DirectSolver(factorization::FactorizationStrategy; kwargs...) = DirectSolver(; factorization, kwargs...)

function VoronoiFVM.SolverControl(strat::DirectSolver, sys; kwargs...)
    return SolverControl(; method_linear = strat.factorization, kwargs...)
end

"""
    GMRESIteration(;factorization=ILUZeroFactorization(), memory=20, restart=true)
    GMRESIteration(factorization; memory=20, restart=true)
    
GMRES Iteration from Krylov.jl via LinearSolve.jl.
"""
Base.@kwdef struct GMRESIteration <: LinearSolverStrategy
    memory::Int = 20
    restart::Bool = true
    factorization::FactorizationStrategy = UMFPACKFactorization()
    blocking::BlockStrategy = NoBlock()
end

function GMRESIteration(factorization::FactorizationStrategy, blocking = NoBlock(); kwargs...)
    return GMRESIteration(; factorization, blocking, kwargs...)
end

function VoronoiFVM.SolverControl(strat::GMRESIteration, sys; kwargs...)
    return SolverControl(;
        method_linear = KrylovJL_GMRES(;
            gmres_restart = strat.memory,
            restart = strat.restart
        ),
        precon_linear = factorizationstrategy(strat.factorization, strat.blocking, sys),
        kwargs...,
    )
end

"""
    CGIteration(;factorization=UMFPACKFactorization())
    CGIteration(factorization)
    
CG Iteration from Krylov.jl via LinearSolve.jl.
"""
Base.@kwdef struct CGIteration <: LinearSolverStrategy
    factorization::FactorizationStrategy = UMFPACKFactorization()
    blocking::BlockStrategy = NoBlock()
end

function CGIteration(factorization::FactorizationStrategy, blocking = NoBlock(); kwargs...)
    return CGIteration(; factorization, blocking, kwargs...)
end

function VoronoiFVM.SolverControl(strat::CGIteration, sys; kwargs...)
    return SolverControl(;
        method_linear = KrylovJL_CG(),
        precon_linear = factorizationstrategy(strat.factorization, strat.blocking, sys),
        kwargs...,
    )
end

"""
    BICGstabIteration(;factorization=UMFPACKFactorization())
    BICGstabIteration(factorization)
    
BICGstab Iteration from Krylov.jl via LinearSolve.jl.
"""
Base.@kwdef struct BICGstabIteration <: LinearSolverStrategy
    factorization::FactorizationStrategy = UMFPACKFactorization()
    blocking::BlockStrategy = NoBlock()
end

function BICGstabIteration(factorization::FactorizationStrategy, blocking = NoBlock(); kwargs...)
    return BICGstabIteration(; factorization, blocking, kwargs...)
end

function VoronoiFVM.SolverControl(strat::BICGstabIteration, sys; kwargs...)
    return SolverControl(;
        method_linear = KrylovJL_BICGSTAB(),
        precon_linear = factorizationstrategy(strat.factorization, strat.blocking, sys),
        kwargs...,
    )
end

"""
    factorizationstrategy(preconditioner, blockstratrgy, system)

Create a factorizations strategy from preconditioner and block information
"""
factorizationstrategy(p::FactorizationStrategy, ::NoBlock, sys) = p

function factorizationstrategy(strat::FactorizationStrategy, ::EquationBlock, sys)
    return BlockPreconditioner(;
        partitioning = partitioning(sys, Equationwise()),
        factorization = factorizationstrategy(strat, NoBlock(), sys),
    )
end

function factorizationstrategy(::ExtendableSparse.ILUZeroPreconditioner, ::PointBlock, sys)
    !isdensesystem(sys) ?
        error("Point block preconditioner needs dense system") : nothing
    return PointBlockILUZeroPreconditioner(; blocksize = num_species(sys))
end

VoronoiFVM.SolverControl(::AbstractStrategy, sys; kwargs...) = SolverControl(; kwargs...)
VoronoiFVM.SolverControl(::Nothing, sys; kwargs...) = SolverControl(; kwargs...)

"""
    Make preconditioner constructors from methods
"""

function (method::LinearSolve.AbstractFactorization)(A)
    pr = LinearProblem(A, zeros(eltype(A), size(A, 1)))
    return init(pr, method)
end

function (method::LinearSolve.SciMLLinearSolveAlgorithm)(A)
    pr = LinearProblem(SparseMatrixCSC(A), zeros(eltype(A), size(A, 1)))
    return init(pr, method)
end

function (f::ExtendableSparse.AbstractFactorization)(A)
    return factorize!(f, A)
end

function LinearAlgebra.ldiv!(u, cache::LinearSolve.LinearCache, b)
    cache.b = b
    sol = solve!(cache)
    return copyto!(u, sol.u)
end
