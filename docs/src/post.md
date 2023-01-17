# Postprocessing

## Plotting

Plotting can be performed using the package [GridVisualize.jl](https://github.com/j-fu/GridVisualize.jl).
This package extends the API with a couple of methods:

```@docs
    GridVisualize.gridplot
    GridVisualize.gridplot!
    GridVisualize.scalarplot
    GridVisualize.scalarplot!
```


## Solution integrals

```@autodocs
Modules = [VoronoiFVM]
Pages = ["vfvm_intergrals.jl"]
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

