# Internal API


Besides of the interface methods for `VoronoiFVMDiffEq`, 
these are not exported and therefore should not be used outside of the package



## Wrapping evaluators for physics callbacks
```@docs 
VoronoiFVM.hasdata
VoronoiFVM.AbstractEvaluator
VoronoiFVM.ResEvaluator
VoronoiFVM.ResEvaluator(physics::Any,symb::Symbol,uproto::Vector{Tv},geom::Any,nspec::Int) where Tv
VoronoiFVM.evaluate!(e::VoronoiFVM.ResEvaluator)
VoronoiFVM.evaluate!(e::VoronoiFVM.ResEvaluator, u::Any)
VoronoiFVM.res(e::VoronoiFVM.ResEvaluator)
VoronoiFVM.ResJacEvaluator
VoronoiFVM.ResJacEvaluator(physics::Any,symb::Symbol,uproto::Vector{Tv},geom::Any,nspec::Int) where Tv
VoronoiFVM.evaluate!(e::VoronoiFVM.ResJacEvaluator, u::Any)
VoronoiFVM.res(e::VoronoiFVM.ResJacEvaluator)
VoronoiFVM.jac(e::VoronoiFVM.ResJacEvaluator)
VoronoiFVM.isnontrivial
```


## Global node and edge assembly loops
```@docs 
VoronoiFVM.AbstractAssemblyData
VoronoiFVM.CellwiseAssemblyData
VoronoiFVM.EdgewiseAssemblyData
VoronoiFVM.nodebatch
VoronoiFVM.noderange
VoronoiFVM.edgebatch
VoronoiFVM.edgerange
VoronoiFVM._fill!
```



## Local node and edge assembly loops

Local assembly methods organize the assembly of data to those degrees of freedom (dofs) which are defined for a given node or edge.
E.g. for an node residual for `nspec` defined species, only those entries need to be assembled into the global residual vector which correspond to actually defined degrees of freedom. 

Similarly for  `nspec x nspec` node Jacobian, an for the `nparam x nspec` parameter derivatives.

These local assembly methods organize the correct loops and call back to the concrete assembly methods passed to them.
These receive global degrees of freedom and the local species numbers to be handled. The callbacks can be used as well for other purposes than assembly

```@docs 
VoronoiFVM.assemble_res_jac
VoronoiFVM.assemble_res
```


## Degree of Freedom management



We distinguish
- active degrees of freedom: these are the actual degrees of freedom 
- degrees of freedom (dof)  potential degrees of freedom - the may be active dofs or dummy ones
  With sparse arrays there are no dummy ones, with dense arrays dummy are maske in the node_dof field
- species: each degree of freedom has associated the species it represents and the node index where it is localized  




```@docs 
VoronoiFVM.isnodespecies
VoronoiFVM.isregionspecies
VoronoiFVM.firstnodedof
VoronoiFVM.lastnodedof
VoronoiFVM.getspecies
VoronoiFVM.getnodedof
```


## Abstract geometry data
```@docs
VoronoiFVM.AbstractGeometryItem
VoronoiFVM.AbstractNode
VoronoiFVM.AbstractNodeData
VoronoiFVM.AbstractEdge
VoronoiFVM.AbstractEdgeData
```

## Global assembly & helpers

```@docs 
VoronoiFVM.eval_and_assemble
VoronoiFVM._addnz
VoronoiFVM._add
```

## Interface methods for `VoronoiFVMDiffEq.jl`
```@docs
VoronoiFVM._eval_res_jac!
eval_rhs!
eval_jacobian!
mass_matrix
prepare_diffeq!
```
