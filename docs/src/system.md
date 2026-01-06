# System

The computational grid required is assumed to correspond to a domain
``\Omega=\cup_{r=1}^{n_\Omega} \Omega_r`` 

Grids for VoronoiFVM are managed by the packages [ExtendableGrids.jl](https://github.com/WIAS-PDELib/ExtendableGrids.jl)
and [SimplexGridFactory.jl](https://github.com/WIAS-PDELib/SimplexGridFactory.jl)


with boundary  ``\partial\Omega=\Gamma=\cup_{b=1}^{n_\Gamma} \Gamma_b``.

The subdomains ``\Omega_r`` are called "regions" and the boundary subdomains ``\Gamma_b`` are called "boundary regions".

On this complex of domains "lives"  a number of species which are either attached to a number of regions or to a number of boundary regions.

All these data, the matrix for the linear system and other things are hold together by a struct `VoronoiFVM.System`. 
This type is not exported to avoid name clashes.


## System constructors

```@docs
VoronoiFVM.System(grid::ExtendableGrid; kwargs...)
VoronoiFVM.System(X::AbstractVector; kwargs...)
VoronoiFVM.System(X::AbstractVector,Y::AbstractVector; kwargs...)
VoronoiFVM.System(X::AbstractVector,Y::AbstractVector,Z::AbstractVector; kwargs...)
update_grid!
physics!
```

## Adding species by species numbers
```@docs
enable_species!(system::VoronoiFVM.AbstractSystem,ispec::Integer, regions::AbstractVector)
enable_species!(system::VoronoiFVM.AbstractSystem; kwargs...)
enable_boundary_species!
VoronoiFVM.is_boundary_species
VoronoiFVM.is_bulk_species
```


## Allocation warnings

The code checks for allocations in the assembly loop. 
Care has been taken to ensure that allocations in the assembly loop don't emerge
from VoronoiFVM.jl code.

If allocations occur in the assembly  loop, they happen in the physics
callbacks.  The corresponding warnings can bee switched off by passing
a  verbosity strings  without  'a'  to the  solver.   If  no data  are
allocated in the physics callbacks, these allocations are probably due to 
type instabilities in physics callbacks.  Type instabilities
can be debugged via the `@time`  macro applied to expressions in a
physics callback.

The following  cases provide some ideas  where to look for  reasons of
the problem and possible remedies:

Case 1: a parameter changes its value, and Julia is not sure about the type.
```julia
eps=1.0

flux(f,u,edge)
    f[1]=eps*(u[1,1]-[1,2])
end
... solve etc ...
eps=2.0
```
Remedy: use a type annotation `eps::Float64=...` to signalize your intent to Julia.
This behaviour is explained in the [Julia documentation](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-captured).



Case 2: variables in the closure have the same name as a variable
introduced in a callback.
```julia
flux(f,u,edge)
    x=(u[1,1]-[1,2])
    f[1]=x
end

... create etc ...

x=solve(...)
```
Remedy: rename e.g. `x=solve()` to `sol=solve()`



## Various tools

```@docs
num_dof
num_species
VoronoiFVM.data
VoronoiFVM.unknowns(system::VoronoiFVM.AbstractSystem; kwargs...)
VoronoiFVM.unknowns(Tu::Type, system::VoronoiFVM.AbstractSystem; kwargs...)
Base.map
Base.map!
VoronoiFVM.isunknownsof
Base.reshape(::AbstractVector, ::VoronoiFVM.AbstractSystem)
```

## Types

```@docs
VoronoiFVM.AbstractSystem
VoronoiFVM.System{Tv,Ti, Tm, TSpecMat<:AbstractMatrix, TSolArray<:AbstractMatrix}
```


## Legacy API
```@docs
boundary_dirichlet!(system::VoronoiFVM.AbstractSystem{Tv}, ispec, ibc, v) where {Tv}
boundary_dirichlet!(system::VoronoiFVM.AbstractSystem; kwargs...)
boundary_neumann!(system::VoronoiFVM.AbstractSystem, ispec, ibc, v)
boundary_neumann!(system::VoronoiFVM.AbstractSystem; kwargs...)
boundary_robin!(system::VoronoiFVM.AbstractSystem, ispec, ibc,alpha, v)
boundary_robin!(system::VoronoiFVM.AbstractSystem; kwargs...)
VoronoiFVM.DenseSystem
VoronoiFVM.SparseSystem
```
