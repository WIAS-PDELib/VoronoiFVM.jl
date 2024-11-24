# Postprocessing

## Plotting

Plotting can be performed using the package [GridVisualize.jl](https://github.com/WIAS-PDELib/GridVisualize.jl).
This package extends the API with a couple of methods for systems:
    
```@docs
GridVisualize.gridplot
GridVisualize.gridplot!
GridVisualize.scalarplot
GridVisualize.scalarplot!
VoronoiFVM.plothistory
```
## Grid verification

```@docs
VoronoiFVM.nondelaunay
```

## Norms & volumes
```@docs
LinearAlgebra.norm
lpnorm
l2norm
w1pseminorm
h1seminorm
w1pnorm
h1norm
lpw1pseminorm
l2h1seminorm
lpw1pnorm
l2h1norm
nodevolumes
```

## Solution integrals
```@docs
VoronoiFVM.integrate(system::VoronoiFVM.AbstractSystem, tf::Vector, U::AbstractMatrix; kwargs...)
VoronoiFVM.integrate(system::VoronoiFVM.AbstractSystem{Tv, Tc, Ti, Tm}, F::Tf, U::AbstractMatrix{Tu}; boundary = false, data = system.physics.data) where {Tu, Tv, Tc, Ti, Tm, Tf}
VoronoiFVM.integrate(system::VoronoiFVM.AbstractSystem, U::AbstractMatrix; kwargs...)
VoronoiFVM.integrate(system::VoronoiFVM.AbstractSystem, F, U::AbstractMatrix; boundary, data)
VoronoiFVM.edgeintegrate
```

## Nodal flux reconstruction
```@docs
nodeflux
```

## Boundary flux calculation
```@autodocs
Modules = [VoronoiFVM]
Pages = ["vfvm_testfunctions.jl"]
```

## Impedance calculatiom
Impedance calculation can be seen as a postprocessing step
after the solution of the unexcited stationary system.


```@autodocs
Modules = [VoronoiFVM]
Pages = ["vfvm_impedance.jl"]
```

