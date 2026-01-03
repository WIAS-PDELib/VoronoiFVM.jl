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
VoronoiFVM.integrate(::VoronoiFVM.AbstractSystem, ::Any, ::AbstractMatrix{Tu}; boundary, data) where Tu
VoronoiFVM.integrate(::VoronoiFVM.AbstractSystem, ::Any, ::VoronoiFVM.AbstractTransientSolution; boundary, data)
VoronoiFVM.integrate(::VoronoiFVM.AbstractSystem, ::AbstractMatrix)
VoronoiFVM.edgeintegrate
```

## Boundary fluxes

### Continuum level motivation 
For  the calculation  of boundary  fluxes, we  use a method inspired by H. Gajewski, WIAS Report No 6, 1993:
"The regularity of the (weak) solutions of the drift-diffusion equations, which are guaranteed by existence statements,
is not sufficient to justify a naive calculation of the currents as boundary integrals. 
Rather, it has also proven to be practical and necessary to reduce contact currents to volume integrals in accordance with the definition of weak solutions using suitable test functions."

We start with  the  continuous  case for a problem
```math
    \partial_t s(u) + \nabla\cdot \vec j(u)  + r(u) - f  =0

```
defined in a domain ``\Omega`` with a boundary ``Γ=⋃_{i\in B}Γ_i`` subdivided into non-overlapping parts
where the set of boundary regions ``B`` is subdivided into ``B=B_N\cup B_0 \cup B_1``.
Assume that ``\vec  j\cdot \vec n=0`` on  ``\Gamma_i`` for ``i\in B_N``.

Let  ``T(x)`` be  a smooth test function  such that
- ``\nabla  T \cdot  \vec  n|_{Γ_i}=  0`` for ``i\in B_N``,
- ``T|_{Γ_i} = 0`` for  ``i\in B_0``,
- ``T|_{Γ_i} = 1`` for  ``i\in B_1``.
Here, we obtain it by solving the Laplace equation 
```math
   -\Delta T =0 \quad \text{in}\; \Omega
```

To obtain as a quantity of interest the  flux ``Q`` through the boundaries ``\Gamma_i`` for ``i\in B_1`` of, calculate:
```math
  \begin{aligned}
    Q=Q(t)&=\sum_{i \in B_1}\int\limits_{\Gamma_i} T\vec j(u)\cdot\vec n \,d\gamma\\
    &=\sum_{i \in B_1}\int\limits_{\Gamma_i} T\vec j(u)\cdot\vec n \,d\gamma 
      + \sum_{i \in B_N}\int\limits_{\Gamma_i} T\vec j(u)\cdot\vec n \,d\gamma 
      + \sum_{i \in B_0}\int\limits_{\Gamma_l} T\vec j(u)\cdot\vec n \,d\gamma\\
    &=\int\limits_{\Gamma} T\vec j(u)\cdot\vec n \,d\gamma\\
    &=\int\limits_{\Omega}\nabla\cdot(T\vec j(u))\,d\omega\\
    &=\int\limits_{\Omega}\nabla T \cdot \vec{j(u)} \,d\omega + \int\limits_{\Omega}T \nabla\cdot \vec{j(u)} \,d\omega\\
    &=\int\limits_{\Omega}\nabla T \cdot \vec{j}(u) \,d\omega - \int\limits_{\Omega}T r(u) \,d\omega +  \int\limits_{\Omega}T f \,d\omega - \int\limits_{\Omega}T \partial_ts(u) \,d\omega
  \end{aligned}
```

### [Discrete approximations of the integrals](@id discrete_appr)
The discete versions of these integrals for evaluation at ``t=t^n`` are as follows:
```math
  \begin{aligned}
    \int\limits_{\Omega}T r(u) \,d\omega &\approx \sum_{k\in N} T_kr(u^n_k)|\omega_k| =: I_{func}(r, T,u^n)\\
    \int\limits_{\Omega}T f \,d\omega &\approx \sum_{k\in N} T_k f_k|\omega_k| =: I_{src}(f, T)\\
    \int\limits_{\Omega}T \partial_ts(u) \,d\omega &\approx \frac{ I_{func}(s,T,u^n) -I_{func}(s,T,u^{n-1})}{t^n-t^{n-1}}\\
    \int\limits_{\Omega}\nabla T \cdot \vec{j}(u^n) d\omega &\approx
    \sum\limits_{\stackrel{k,l}{\partial \omega_k \cup \partial \omega_l\neq\emptyset}} (T_{k}-T_{l})\frac{\sigma_{kl}}{h_{kl}}g(u^n_k, u^n_l) =: I_{flux}(j,T,u)
  \end{aligned}
```

### General flux calculation API

```@docs
TestFunctionFactory
TestFunctionFactory(::VoronoiFVM.AbstractSystem; kwargs...)
testfunction
VoronoiFVM.integrate(::VoronoiFVM.AbstractSystem,::Vector{Tv},::AbstractMatrix{Tu}; kwargs...) where{Tu, Tv}
VoronoiFVM.integrate(::VoronoiFVM.AbstractSystem,::Vector,::VoronoiFVM.AbstractTransientSolution; kwargs...)
```

### Special purpose cases
```@docs
VoronoiFVM.integrate(::VoronoiFVM.AbstractSystem,::Any, ::AbstractMatrix{Tv},::AbstractMatrix{Tv}, ::Any; kwargs...) where {Tv}
VoronoiFVM.integrate_TxFunc
VoronoiFVM.integrate_TxSrc
VoronoiFVM.integrate_∇TxFlux
VoronoiFVM.integrate_TxEdgefunc
VoronoiFVM.integrate_stdy
VoronoiFVM.integrate_tran
```


## Nodal flux reconstruction
```@docs
nodeflux
```


## Impedance calculatiom
Impedance calculation can be seen as a postprocessing step
after the solution of the unexcited stationary system.


```@autodocs
Modules = [VoronoiFVM]
Pages = ["vfvm_impedance.jl"]
```
