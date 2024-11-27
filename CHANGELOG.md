# Changes

## v3.0 Planned  (pending some improvements in LinearSolve)
  - Remove VoronoiFVM solver strategies
  - Stop re-exporting ForwardDiff.value
  - Try to remove type piracies
  - Remove `params` from edge, node structs (apparently never used)

## v2.6.1 November 27, 2024
  - Allow to calculate node flux reconstruction from general function
  - Runic formatting 
  - pre-commit checks
  
## v2.5.0 November 17, 2024
  - update show methods for physics, grid
  
## v2.4.0 November 11, 2024
  - Use `precs` based linear solver API, see https://github.com/SciML/LinearSolve.jl/pull/514 with  ExtendableSparse 1.6
  - Deprecate VoronoiFVM solver strategies

## v2.3.0 November 5, 2024
  - Allow to use result of time embedding a steady state for impedance
  - Move to WIAS-PDELib org

## v2.2.0 October 30, 2024
  - Add `params` to SystemState, allow to pass params to ODEProblem
  - Fix use of end results of time evolution as steady state for impedance calculations
  - new internal solutionarray method. 

## October 28, 2024

Moved repository from https://github.com/j-fu/VoronoiFVM.jl to https://github.com/WIAS-PDELib/VoronoiFVM.jl. [WIAS-PDELib](https://github.com/WIAS-PDELib/) is a github organization created to collectively manage the Julia packages developed under the lead of the [WIAS Numerical Mathematics and Scientific Computing](https://wias-berlin.de/research/rgs/fg3)  research group. According to the [github docs on repository transfer](https://docs.github.com/en/repositories/creating-and-managing-repositories/transferring-a-repository#whats-transferred-with-a-repository), all links to the previous repository location are automatically redirected to the new location, and all relationships with forks stay intact.

## v2.1.0 October 23, 2024
  - Features:
    - Add VoronoiFVMExtendableFEMBaseExt  extension allowing to use velocity fields calculated  with [ExtendableFEM](https://github.com/chmerdon/ExtendableFEM.jl) in convection-diffusion problems

## v2.0.2 October 12, 2024
  - Bugfixes
    - replace internally used `values` and `vec` methods by new `dofs`
    - fix #119, update Example115, Example120, test sparse unknown storage & ODE solver
    - directly re-export ForwardDiff.value, dont type-pirate a new method

## v2.0.1 October 7, 2024
  - Bugfixes
    - fix case where the last element of a time solution is used as initial value for a next `solve`

## v2.0.0 September 12, 2024

### Breaking:
   - Removed deprecated solve methods:
````
    solve(inival, system::VoronoiFVM.AbstractSystem, lambdas;kwargs...)
    solve(inival, system::VoronoiFVM.AbstractSystem; kwargs...)
    solve!(solution,inival, system::VoronoiFVM.AbstractSystem; kwargs...)
````
   - History is now a field of the solution but not system, to be accessed with  ``history(solution)``, instead of ``history(system)`` both for stationary and transient solutions
   - [`VoronoiFVM.Physics`](@ref) callbacks (`flux`, `storage`, etc.) need a `data` argument, even if `data` is not provided during ``Physics`` or ``System`` construction

### Added
  - Introduced ``SystemState`` which contains entries (matrix, residuals, data) which before were
    part of ``System``. One can have different states with different data for the same system.
  - Introduced `solve!(state::SystemState; kwargs...)` method
  - Stationary solutions now contain a history entry and are subtypes of VoronoiFVM.AbstractSolutionArray
  - Changelog now in package root

### Modified
  - Sparse solver default for Float64 in all space dimensions is now  UMFPACKFactorization. For all other
    number types it is SparspakFactorization. Before it was KLU for 1D, Sparspak for 2D an UMFPACK for 3D.

## v1.25.0, Sept 9, 2024
- [`VoronoiFVM.Physics`](@ref) callbacks (`flux`, `storage`, etc.) without `data` argument are now deprecated

## v1.24.0 August 20, 2024
- Bernoulli function using expm1 (thanks @tcaduser !)

## v1.23.0 August 2, 2024
- Add `integrate(sys,solution)` method

## v1.22.0 July 17, 2024
- Multithreaded assembly

## v1.21.0 July 1, 2024
- Introduce nondelaunay method for checking Delaunay properties of grid

## v1.20.0 April 26 2024
- Drop Julia 1.6 support, improve package quality (partial Aqua tests, explicit imports)

## v1.19.0 Feb 01 2024
- Enable equation block preconditioning for sparse unknown storage

## v1.18.0 Feb 01 2024
- Re-shoring of OrdinaryDiffEq interface, no need of VoronoiFVMDiffEq.jl anymore
  It appeared that it is sufficient to depend on SciMLBase for this, and
  all tests can be done with OrdinaryDiffEq.jl

## v1.17.1 Jan 30, 2024
- Bugfix for boundary node factors
- Bugfix with types for RecursiveArrayTools

## v1.16.0 Dec 15, 2023
- Bugfix for assembly of outflow bc
- Bugfix for matrixtype=:banded
- Updated [`plothistory`](@ref):

## v1.15.0 Dec 1, 2023
- Adjusted time/embedding stepping scheme, added `num_final_steps` to  [`VoronoiFVM.SolverControl`](@ref). This may lead to
  slightly different results when solving time dependent problems. Some unit test values have been adapted. Before,
  accidentally, very small time steps at the end of an evolution were possible.

## v1.14.0 Nov 27, 2023
- Add `Δu_max_factor` and `Δt_decrease` to [`VoronoiFVM.SolverControl`](@ref).
  Before (with values 2.0 and 0.5, respectively) the were fixed. New `Δu_max_factor=1.2` default.
- Add history to [`VoronoiFVM.TransientSolution`](@ref), prepare deprecation of `system.history`
- Add [`plothistory`](@ref) method

## v1.13.0 July 24, 2023
- Add [`nodevolumes`](@ref) method

## v1.12.0 July 22, 2023
- Add functionality for outflow boundary conditions

## v1.11.0 July 17, 2023
- Add [`calc_divergences`](@ref) method for checking velocity field divergences
- Fix form factor calculation and velocity projecion for unstructructured grids and cylindrical symmetry

## v1.10.0 July 11, 2023
- Use AbstractTransientSolution in gridvisualize stuff

## v1.9.0 June 27, 2023
- With `control.handle_exceptions=true`, in case of a failing step,
  time stepping and embedding now returns the solution calculated so far instead of
  throwing an error

## v1.8.0 June 20, 2023
- LinearSolve 2.x

## v1.7.0 May 17, 2023
- Discrete Sobolev norms

## v1.6.0 May 12, 2023
- Rework linear solver strategies - prevent combinatorial explosion.
  E.g. `gmres_iluzero` is now `GMRESIteration(ILUZeroPreconditioner())` etc.

## v1.5.0 May 9, 2023
- Introduced solver strategies like `gmres_iluzero()`, `direct_umfpack()` etc.
  See documentation of the module `VoronoiFVM.SolverStrategies`. More to come.
- edgewise assembly - faster in particular for 3D
- Plan: edgewise assembly seems to be non-breaking,
  if this is confirmed, will be made default in 1.6 or (if it appears to be breaking
  for some) in 2.0.

## v1.4.0 May 3, 2023
- equation-block preconditioning support with the help of ExtendableSparse.jl

## v1.3.0 April 13, 2023
- inplace_linsolve! for dense linear system solution in flux functions
- Mixture flow example 510

## v1.2.0 March 17, 2023
- Initialization of quantities, create unknowns using Base.map

## v1.1.0 Feb 22, 2023
- Edge reactions, Joule heating

## v1.0.0 Feb 22, 2023
- full LinearSolve compatibility

## v0.19.0 Jan 31, 2023
This is a breaking release. Implementations using default solver settings should continue to work without modifications,
albeit possibly showing  deprecation and  allocation warnings.
Really breaking is control of iterative linear solvers and allocation checks.

### Breaking
- Make `solve` a method of `CommonSolve.solve` (and re-export it).
- Rely on `LinearSolve.jl` for linear system solution including control of iterative solvers.
- New verbosity handling. `verbose` can now be a Bool or a String of flag characters, allowing for control of different output categories. I would love to do this via  logging, but there is still a [long way to go](https://github.com/JuliaLang/julia/issues/33418) IMHO
- Allocation check is active by default with warnings instead of throwing an error. These warnings can be muted by passing a `verbose` string without 'a'. This is now the only control in this respect. All `check_allocs` methods/kwargs, control via environment variables have been removed.
- Deprecation warnings can be switched off by passing a `verbose` string without 'd'.
- Improve iteration logging etc., allow for logging of linear iterations ('l' flag character)


### Deprecations
- Deprecated all `VoronoiFVM.solve` methods with signatures others than `solve(system; kwargs...)` which renders them incompatible to the philosophy of `CommonSolve. Updated examples accordingly.
- Deprecated the following entries of SolverControl/solve kwargs:
   `:tol_absolute` => `:abstol`,
    `:tol_relative` => `:reltol`,
    `:damp` => `:damp_initial`,
    `:damp_grow` => `:damp_growth`
    `:max_iterations` => `e:maxiters`
    `:tol_linear` => `:reltol_linear`
    `:max_lureuse` =>
- `NewtonControl`


## v0.18.8 - 0.18.10  Dec 11, 2022 - Jan 5, 2023
- Internal restructuring: remove `@create_physics_wrappers` macro, reduce boilerplate in assembly, wrap repeating pattenrns
  into functions.  The price in the moment is a bit of a slowdown of assembly.
- Fix parameter dependency handling (now we can get parameter derivative without solving in dual numbers; see
  the runh() example in Example430. However in the moment the advantatge is not very clear, so this  is
  on hold...
## v0.18.7  Dec 7, 2022
- bump gridvisualize compat
## v0.18.6  Dec 3, 2022
- enable non-diagonal mass matrices for VoronoiFVMDiffEq
## v0.18.5  Nov 30, 2022
- ready for Julia 1.9, re-enable CI on nighly
## v0.18.4  Nov 29, 2022
- Add API methods used by VoronoiFVMDiffEq.jl
## v0.18.3  Oct 18 2022
- Removed some allocations
## v0.18.2  Oct 13 2022
- Emerging capability of differencing wrt. parameters (experimental, see example 430)
- Allow iterative methods from Krylov.jl
- Proper Dirichlet initialization with bcondition
- Allow for more general matrix structures (banded, tridiagonal, multidiagonal)

## v0.18.1 September 25 2022
- Working spherical symmetry case
## v0.18 September 12 2022
- Remove DifferentialEquations interface (move this to a glue package)
The current method of activating it through require is too brittle with
respect to versioning.

## v0.17.1 August 20 2022
- Fix DifferentialEquations interface, start transition to LinearSolve

## v0.17.0 July 1 2022
- ensure not to assemble data for species where they are not enabled
  This change should be breaking only for incorrect code where physics
  callbacks write into degrees of freedom which are not enabled

## v0.16.5 June 30, 2022
- add `iteration` to solver options, allow to choose :cg, :bicgstab.
- allow setting penalty with boundary_dirichlet!

## v0.16.4 May 25, 2022
- fix x-t plots

## v0.16.3 March 18, 2022
- Linearization API
- relax some type constraints

## v0.16.2 Feb 18, 2022
- ExtendableGrids 0.9

## v0.16.1 Feb 17, 2022
- fix quantity postprocessing
- define unknown access for abstract vectors instead of vectors
- pass rhs/unknowns wrappers in postprocessing methods
- integrals as a wrapper type with proper quantity handling

## v0.16.0 Feb 13, 2022
- Expose ODEProblem (and possibly ODEFunc) from VoronoiFVM.System.
- __Breaking__: Remove `solve` wrapper for DifferentialEquations.solve, instead recommend to call that directly
- __Breaking__: Handle DifferentialEquations.jl via Requires.jl.

## v0.15.1 Jan 15, 2022
- Documentation fixes
- Fix OrdinaryDiffEq interface
- added example for current calculation with Quantities
- Fixed type instabilities in quantities interface

## v0.15.0 Jan 1, 2022
- __Breaking__: History is not anymore returned by `solve`, instead it can be accessed via [`history`](@ref) after the solution.
- Cleaned API:
  - [`VoronoiFVM.solve(system::VoronoiFVM.AbstractSystem; kwargs...)`](@ref) is now the main method of `solve` which allows to
    access stationary, transient, embedding and DifferentialEquations based solvers.
  - Joint implementation for implicit Euler timestepping and parameter embedding
  - Handle more kwargs via SolverControl (e.g. log)
  - Use Parameters.jl in struct definition
  - Add history types [`NewtonSolverHistory`](@ref), [`TransientSolverHistory`](@ref)
  - `detailed` and `summary` methods for both history types
  - Nonlinear solver example notebook (under development): [nonlinear-solvers.jl](@ref nonlinear-solvers)
- OrdinaryDiffEq solver now in CI
- scalarplot for 1D transient solutions
- Sparsity detection via Symbolics.jl instead of the sunsetted SparsityDetection.jl

## v0.14.0 Dec 24, 2021
Backward compatible, hopefully nonbreaking API simplification
- Boundary conditions are now specified in breaction.
  Advantages:
   - easy x/t dependency
   - unified (upcoming) interface for parameters
   - unified handling of standard and nonstandard boundary conditions
   - simpler documentation
- Made NewtonControl alias of SolverControl, continue to work with SolverControl
- System constructor now directly takes physics callback functions, no need anymore to work with extra physics struct
- solve() now takes "SolverControl" parameters as kwargs,no need anymore to work with extra NewtonControl/SolverControl struct
- Notebooks as part of documentation and CI
- See also the pluto notebook [api-update.jl](@ref api-update)

## v0.13.2 Oct 29, 2021
- Bernoulli function overhaul
## v0.13.1
- sorted things with ExtendableGrids
- nodal flux reconstruction (e.g. for visualization)
## v0.13.0, Oct 13, 2021
- various bug fixes, explicit numbering of edge nodes

## v0.12.3, July 7, 2021
- Add quantity id
- Document quantities

## v0.12.2, July 7, 2021
Introduce the notion of quantities which can be continuous or discontinuous at interfaces.

- Quantity handling is implemented on top of species handling
- Unknowns u and rhs y now passed to callbacks as wrapper types,
  and can be indexed by quantity or by species numbers.
  Moreover, this will allow to abstract parameters, gradients etc. in future versions.

## v0.12.0, July 2 2021
-  By default, the `u` parameter in flux callbacks is now a  `nspec x 2` array
- `unknowns(edge,u)`, `viewK, viewL` are obsolete, they still work for backward compatibility
- `physics.num_species` is now meaningless, num_species is automatically detected.
- `SparseSystem` and `DenseSystem` are now type aliases of a parametrized type instead of two
   independent subtypes of `System`

## v0.11.8
- increase chunk size threshold to match argument length in calls to vector_mode_jacobian

## v0.11.7
- First attempts on surface flux

## v0.11.1, April 13, 2021
- Assembly loops cleaned from type instabilities
- Optionally check for allocations due to type instabilities introduced in physics callbacks.
  See [`check_allocs!`](@ref) for more information.

## v0.11, April 12, 2021
- Depending on Julia 1.5 now
- Lineaer solvers now based on factorize! from ExtendableSparse 0.5
- Documentation overhaul
- Re-checking impedance calculation

## v0.10.13 April 1, 2021
- Outflow boundary conditions

## v0.10.8 March 22, 2021
- TransientSolution structure, transient solve
- Solve compatible with DifferentialEquations.jl

## v0.10.3 Feb 11, 2021
- Introduce non-mutating solve
- Optionally record history if log kw is true in solve.

## v0.10.0 Jan 9 2021
- Moving visualization to the package [GridVisualize.jl](https://github.com/j-fu/GridVisualize.jl), emerging from
  the visualization methods in ExtendableGrids


## v0.9.0 Dec 21 2020
- Add the possibility to interface with [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/)
- __Breaking__: The API change to passing the unknowns to the
  an edge callback as a matrix turned out to be a dead end in the
  strategic sense. In order to extend functionality, we need to be able
  to pass more data to  which we can apply differetiation. Particular plans
  involve bifurcation parameters and reconstructed gradients.
  So we return to the viewK/viewL pattern we had before. However,
  these are now aliases:
  - `viewK(edge,u)=unknowns(edge,u,1)`
  - `viewL(edge,u)=unknowns(edge,u,2)`
  In order to ease refactoring in the case where models have been
  implemented with Matrix access to the unknowns, `unknowns(edge,u)`
  returns a matrix of the edge unknowns.
  For refactoring, just rewrite e.g.
```julia
    function flux(y,u,edge)
        for ispec=1:nspec
            y[ispec]=u[ispec,1]-u[ispec,2]
        end
    end
```
  to
```julia
    function flux(y,u0,edge)
        u=unknowns(edge,u0)
        for ispec=1:nspec
            y[ispec]=u[ispec,1]-u[ispec,2]
        end
    end
```


## v0.8.5 Sep 1 2020
- allow any object in Physics.data  (thanks Jan Weidner)
- add generic operator for non-canonical problem structures

## v0.8.4 July 25 2020
- Update ExtendableGrids + ExtendableSparse

## v0.8.3 June 25 2020
- Replace splatting by dispatch on availability of data record

## v0.8.2 May 15 2020
- Form factors are now pre-calculated and stored
- Introduced update_grid! for triggering re-calculation
  if coordinates have changed

## v0.8.1 May 2 2020
- Introduce evolve! : time solver with automatic timestep control

## v0.8 Apr 28, 2020
- Replaced VoronoiFVM grid module by  [ExtendableGrids.jl](https://github.com/j-fu/ExtendableGrids.jl)
- Moved grid generation, modification, plotting  over to ExtendableGrids
- Necessary changes in codes using VoronoiFVM:
  - Replace `grid.coord` by `coord` obtained via `coord=coordinates(grid)` or  `coord=grid[Coordinates]` after importing ExtendableGrids
  - Replace `VoronoiFVM.plot` by `ExtendableGrids.plot`.
  - In the plot method, Plotter is now a keyword argument
  - `VoronoiFVM.Grid()` now returns a ExtendableGrids.ExtendableGrid,
    in fact it is just an alias to [ExtendableGrids.simplexgrid](https://j-fu.github.io/ExtendableGrids.jl/stable/search/?q=simplexgrid)
  - For using any methods on grids like cellmask! one needs to use `ExtendableGrids`
  - Subgrids now are of the same type `ExtendableGrids`,  views are currently defined for vectors only.

## v0.7 Feb 28 2020
- API modification:
  - __Breaking__:
    - `data` parameter passed to physics callbacks only if `Physics` object is created with `data` parameter.

      This makes the API more consistent in the case that parameters are just taken from the closure (the scope where the physics functions are defined)
      and no data `object` has been created.
    - Replace `node.coord[i]` by  `node[i]`.
    - Replace `edge.coordK[i]` by  `edge[i,1]`.
    - Replace `edge.coordL[i]` by  `edge[i,2]`.

      This now directly accesses the coordinate field of the grid and avoids copying of the coordinates
  - Backward compatible:
    - No need for `viewK and viewL` in edge callbacks (they also make trouble with allocations...)
        - Replace `uk[i]` by `u[i,1]`
        - Replace `ul[i]` by `u[i,2]`
     - Replace `VoronoiFVM.DenseSystem(...)`  by `VoronoiFVM.System(..., unknown_storage=:dense)`
     - Replace `VoronoiFVM.SparseSystem(...)` by `VoronoiFVM.System(..., unknown_storage=:sparse)`

- No allocations anymore in assembly loop:
  - Replaced `ElasticArray` in `Grid` by normal one - this was the largest regression
  - Return `nothing` from mutating methods to avoid some allocations
  - Indexing in `formfactors.jl` with `Int`

## v0.6.5 Jan 25 2020
- use updateindex! for matrix, depend on ExtendableSparse 0.2.6

## v0.6.4 2020-01-20
- Rearranged + commented boundary assembly loop
- Reworked + renamed some examples
- Document that unknowns doesn't initialize values.

## v0.6.3 2019-12-21
remove xcolptrs call
Update dependency on ExtendableSparse

## v0.6.2 2019-12-20
Updated dependency list (Triangulate ^0.4.0)

## v0.6.1, 2019-12-17
- return "plotted" for being able to  place colormap
- require Triangulate >= 0.3.0

## v0.6.0, Dec 15 2019
- Removed Triangle submodule, depend on new Triangulate.jl Triangle wrapper
- link to source code in examples
- boundary_dirichlet! etc methods for setting boundary conditions
## v0.5.6 Dec 5 2019
- Bug fixes
- check triangle input for min 3 points
- check triangle edgelist for C_NULL
- voronoi plot
## v0.5.5 Dec 4 2019
- (Temporary) Copy of TriangleRaw as Triangle submodule. To be replaced by dependency on evisioned package
## v0.5.4 Dec 3 2019
- Re-enabled ElasticArrays in grid structure (for the time being)
- Added potkink example: this adds an inner boundary
## v0.5.3 Dec 1 2019
- triangle in optional submodule
- Modified API for plotting
   - Removed formal dependency on Plots and PyPlot
   - Use Plotter module as first parameter to plot methods  - replaces fvmplot
     and fvmpyplot functions. Use `VoronoiFVM.plot(PyPlot,...)` resp.  `VoronoiFVM.plot(Plots,...)`
   - No more complaints when package is used in environment with plots or pyplot installed
- Modified API for impedance

## v0.5.2 Nov 19, 2019
- Reorganized grid stuff
- Included triangle (after Ideas from TriangleMesh.jl)

## v0.5.1 Nov 13, 2019
- Fixed performance regression: AbstractArrays for Grid components were slow.
- Added handling of cylindrical coordinates

## V0.5, November 10, 2019
- Velocity projections
- Added edge handling to grid struct

## V0.4.2, November 6, 2019
- Replaced PyPlot by Plots
- Better and more examples

## V0.4, July 12, 2019
- Registered with Julia ecosystem
- Enhance Newton solver by embedding, exception handling
- Replace SparseMatrixCSC with ExtendableSparseMatrix
- fixed allocation issues in assembly
- assured that users get allocation stuff right via
  typed functions in physics structure
- more julianic API

## V0.3, April 9 2019
- Renamed from TwoPointFluxFVM to  VoronoiFVM
- Complete rewrite of assembly allowing sparse or dense matrix
  to store degree of freedom information
    - Solution is a nnodes x nspecies sparse or dense matrix
    - The wonderful array interface of Julia still provides slicing
      etc in order to access species without need to write
      any bulk_solution stuff or whatever when using the sparse variant
- Re-export value() for debugging in physics functions
- Test function handling for flux calculation
- First working steps to impedance handling
- Abolished Graph in favor of  Grid, Graph was premature optimization...

## V0.2, Feb 20, 2019

- Changed signature of all callback functions:
  This also allows to pass user defined arrays etc. to the callback functions.
  In particular, velocity vectors can be passed this way.

  - Besides of `flux!()`, they now all have `node::VoronoiFVM.Node`
    as a second argument.

  - `flux!()` has `edge::VoronoiFVM.Edge` as a second argument

  - the `x` argument in `source!()` is omitted, the same data
     are now found in `node.coord`


- New method `edgelength(edge::VoronoiFVM.Edge)`

## V0.1, Dec. 2018

- Initial release
