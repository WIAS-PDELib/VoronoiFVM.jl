# Solvers

The package comes with a [built-in solve method](@ref "Built-in solver") which solves 
stationary problems, homotopy embedding problems and transient problems 
via the implicit Euler method. In particular, the transient solver allows to use
nonlinear storage terms.

Alternatively, [OrdinaryDiffEq.jl based solvers](@ref diffeq) can be used 
for transient problems.


## Built-in solver
This solver and its default parameters are tuned for robustness,
possibly at the expense of solution speed. Careful tuning of the parameters, or -- in the case of transient problems --
the choice of a [OrdinaryDiffEq.jl based solver](@ref diffeq) can significantly improve the performance.

Overview:
- [Solve method](@ref "Solve method")
- [Solver control](@ref "Solver control")
- [System state](@ref "System state")
- [Linear solver control](@ref "Linear solver control")
- [History handling](@ref "History handling")
- [Matrix extraction](@ref "Matrix extraction")

### Solve method
```@docs
VoronoiFVM.solve(system::VoronoiFVM.AbstractSystem; kwargs...)
``` 
### Solver control
```@docs 
SolverControl
fixed_timesteps!
```
### System state

```@docs
VoronoiFVM.SystemState
VoronoiFVM.SystemState(::Type, system::VoronoiFVM.AbstractSystem; data)
VoronoiFVM.SystemState(system::VoronoiFVM.AbstractSystem; data)
VoronoiFVM.solve!(state::VoronoiFVM.SystemState; kwargs...)
Base.similar(state::VoronoiFVM.SystemState; kwargs...)
```

### Linear solver control
Linear systems are solved using LinearSolve.jl. 
Linear solve compatible solver strategies (factorizations, iterative solvers) can be
specified wia  `method_linear` keyword argument to LinearSolve (equivalent to the `method_linear` 
entry of [`SolverControl`](@ref).

Currently supported possibilities are documented in the
[documentation of ExtendableSparse.jl](https://wias-pdelib.github.io/ExtendableSparse.jl/stable/linearsolve/#Solving-with-LinearSolve.jl).

VoronoiFVM.jl provides partitioning methods for block preconditioners.
```@docs
Equationwise
partitioning
```


### History handling
If `log` is set to true in `solve`, the history of newton iterations and  time/embedding
steps is recorded and returned as `history(solution)`

```@docs
NewtonSolverHistory
TransientSolverHistory
VoronoiFVM.DiffEqHistory
Base.summary(::NewtonSolverHistory)
Base.summary(::TransientSolverHistory)
details
history

history_details
history_summary
```



### Matrix extraction
For testing and teaching purposes, one can obtain residual and linearization at a given vector of unknowns

```@docs
evaluate_residual_and_jacobian
```

## [OrdinaryDiffEq.jl transient solver](@id diffeq)

For transient problems, as an alternative to the [built-in implicit Euler method](@ref "Built-in solver"), (stiff) ODE solvers from 
[OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl)  can be used.

The interface just provides two methods: creation of an [ODEProblem](https://diffeq.sciml.ai/stable/basics/overview/#Defining-Problems) from a [`VoronoiFVM.System`](@ref) and a reshape method
which turns the output of the ode solver into a [`TransientSolution`](@ref).

The basic usage pattern is as follows: use [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) and replace the call to the built-in solver
```julia
sol=solve(sys; times=(t0,t1), inival=inival)
```
by
```julia
problem = ODEProblem(sys,inival,(t0,t1))
odesol = solve(problem, solver)
sol=reshape(odesol,sys)
```
Here, `solver` is some  ODE/DAE solver from [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl).
It is preferable to choose methods able to handle [stiff problems](https://diffeq.sciml.ai/stable/solvers/ode_solve/#Stiff-Problems).
Moreover, often, discretized PDE systems (e.g. containing elliptic equations) are differential agebraic equation (DAE) systems 
which should be solved by [DAE solvers](https://diffeq.sciml.ai/stable/solvers/dae_solve/).
Some choices to start with are Rosenbrock methods like 
[Rosenbrock23](https://docs.sciml.ai/DiffEqDocs/stable/solvers/dae_solve/#Rosenbrock-W-Methods)
and multistep methods like [QNDF](https://docs.sciml.ai/DiffEqDocs/stable/solvers/dae_solve/#Multistep-Methods)
and [FBDF](https://docs.sciml.ai/DiffEqDocs/stable/solvers/dae_solve/#Multistep-Methods).

If the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)
package is loaded, the `solver` parameter can be omitted, and some default is chosen.

The solution `odesol` returned by `solve` conforms to the [ArrayInterface](https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/#Array-Interface)
but "forgot" the VoronoiFVM species structure. Using 

Accessing `odesol(t)` will return an interpolated solution vector giving
the value of the solution at moment `t`. Using [`reshape(::AbstractVector, ::VoronoiFVM.AbstractSystem)`](@ref) on `(odesol(t),system)` it can be turned into into a
sparse or dense array reflecting the species structure of `system`. The order of the [interpolation](https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/#Interpolations-and-Calculating-Derivatives)
depends on the ODE solver.

Using [`reshape(::AbstractDiffEqArray,::VoronoiFVM.AbstractSystem)`](@ref) on `(odesol, system)` returns a [`TransientSolution`](@ref) knowing
the species structure.

```@docs
SciMLBase.ODEProblem
reshape(::AbstractDiffEqArray,::VoronoiFVM.AbstractSystem)
SciMLBase.ODEFunction
```



## Legacy API
During the development of the code, a number of API variants have been developed which 
are supported for backward compatibility.

```@docs
NewtonControl
``` 

