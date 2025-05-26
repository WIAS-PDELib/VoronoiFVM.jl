# # 510: Mixture
# ([source code](@__SOURCE_URL__))
#=

Test mixture diffusion flux. The problem is here that in the flux function we need to
solve a linear system of equations which calculates the fluxes from the gradients.

``u_i`` are the species partial pressures, ``\vec N_i`` are the species fluxes.
``D_i^K`` are the Knudsen diffusion coefficients, and ``D^B_{ij}`` are the binary diffusion coefficients.
```math
  -\nabla \cdot \vec N_i =0 \quad (i=1\dots n)\\
  \frac{\vec N_i}{D^K_i} + \sum_{j\neq i}\frac{u_j \vec N_i - u_i \vec N_j}{D^B_{ij}} = -\vec \nabla u_i \quad (i=1\dots n)
```
From this representation, we can derive the matrix ``M=(m_{ij})`` with
```math
m_{ii}= \frac{1}{D^K_i} + \sum_{j\neq i} \frac{u_j}{D_ij}\\
m_{ij}= -\sum_{j\neq i} \frac{u_i}{D_ij}
```
such that 
```math
	M\begin{pmatrix}
\vec N_1\\
\vdots\\
\vec N_n
\end{pmatrix}
=
\begin{pmatrix}
\vec \nabla u_1\\
\vdots\\
\vec \nabla u_n
\end{pmatrix}
```
In the two point flux finite volume discretization, this results into a corresponding linear system which calculates the discrete edge fluxes from the discrete gradients.  Here we demonstrate how to implement this in a fast, (heap) allocation free way.

For this purpose, intermediate arrays have to be used. They need to have the same element type as the unknowns passed to the flux function
(which could be Float64 or some dual number). 

To do so without (heap) allocations can be achieved at least in three ways tested in this example:
- Stack allocation within the flux function using [`StrideArrays`](https://github.com/JuliaSIMD/StrideArrays.jl)`.StrideArray`, with the need to have static (compile time) information about the size of the local arrays to be allocated via e.g. a global constant, or, as demonstrated here, a type parameter. As [documented in  StrideArrays.jl](https://juliasimd.github.io/StrideArrays.jl/stable/stack_allocation/), use `@gc_preserve` when passing a `StrideArray` as a function parameter. See also [this Discourse thread](https://discourse.julialang.org/t/what-is-stridearrays-jl/97146).
- Stack allocation within the flux function using [`StaticArrays`](https://github.com/JuliaArrays/StaticArrays.jl)`.MArray`, with the need to have static (compile time) information about the size of the local arrays to be allocated via e.g. a global constant, or, as demonstrated here, a type parameter. However, this may run into  [this issue](https://github.com/JuliaArrays/StaticArrays.jl/issues/874), requiring `@inbounds` e.g. with reverse order loops.
- Preallocation using [`PreallocationTools`](https://github.com/SciML/PreallocationTools.jl)`.DiffCache`. While this avoids the need to pass the size via a compile time constant, one has to ensure that each running thread has its own cache. Here this is achieved by providing a cache for each partition.

=#

module Example510_Mixture

using Printf
using VoronoiFVM

using ExtendableGrids
using GridVisualize
using LinearAlgebra
using AMGCLWrap
using Random
using StrideArraysCore: @gc_preserve, StrideArray, StaticInt, PtrArray
using LinearSolve, ExtendableSparse
using ExtendableSparse: ILUZeroPreconBuilder
using StaticArrays
using ExtendableSparse
using PreallocationTools
using Metis

## Userdata structure for passing number of species as parameter known at compile time.
## Buffers are stack allocated 
Base.@kwdef struct MyDataStaticSizeInfo{NSpec}
    DBinary::Symmetric{Float64, Matrix{Float64}} = Symmetric(fill(0.1, NSpec, NSpec))
    DKnudsen::Vector{Float64} = ones(NSpec)
    diribc::Vector{Int} = [1, 2]
end
nspec(::MyDataStaticSizeInfo{NSpec}) where {NSpec} = NSpec
MyDataStaticSizeInfo(nspec;  kwargs...) = MyDataStaticSizeInfo{nspec}(; kwargs...)

## Flux with stack allocated buffers using StrideArray
function flux_strided(f, u, edge, data)
    T = eltype(u)
    M = StrideArray{T}(undef, StaticInt(nspec(data)), StaticInt(nspec(data)))
    au = StrideArray{T}(undef, StaticInt(nspec(data)))
    du = StrideArray{T}(undef, StaticInt(nspec(data)))
    ipiv = StrideArray{Int}(undef, StaticInt(nspec(data)))

    for ispec in 1:nspec(data)
        M[ispec, ispec] = 1.0 / data.DKnudsen[ispec]
        du[ispec] = u[ispec, 1] - u[ispec, 2]
        au[ispec] = 0.5 * (u[ispec, 1] + u[ispec, 2])
    end

    for ispec in 1:nspec(data)
        for jspec in 1:nspec(data)
            if ispec != jspec
                M[ispec, ispec] += au[jspec] / data.DBinary[ispec, jspec]
                M[ispec, jspec] = -au[ispec] / data.DBinary[ispec, jspec]
            end
        end
    end

    ## Pivoting linear system solution via RecursiveFactorizations.jl (see vfvm_functions.jl)
    inplace_linsolve!(M, du, ipiv)

    for ispec in 1:nspec(data)
        f[ispec] = du[ispec]
    end
    return
end

## Flux with stack allocated buffers using MArray
function flux_marray(f, u, edge, data)
    T = eltype(u)
    n = nspec(data)

    M = MMatrix{nspec(data), nspec(data), T}(undef)
    au = MVector{nspec(data), T}(undef)
    du = MVector{nspec(data), T}(undef)
    ipiv = MVector{nspec(data), Int}(undef)

    for ispec in 1:nspec(data)
        M[ispec, ispec] = 1.0 / data.DKnudsen[ispec]
        du[ispec] = u[ispec, 1] - u[ispec, 2]
        au[ispec] = 0.5 * (u[ispec, 1] + u[ispec, 2])
    end

    for ispec in 1:nspec(data)
        for jspec in 1:nspec(data)
            if ispec != jspec
                M[ispec, ispec] += au[jspec] / data.DBinary[ispec, jspec]
                M[ispec, jspec] = -au[ispec] / data.DBinary[ispec, jspec]
            end
        end
    end

    ## Pivoting linear system solution via RecursiveFactorizations.jl (see vfvm_functions.jl)
    inplace_linsolve!(M, du, ipiv)

    for ispec in 1:nspec(data)
        f[ispec] = du[ispec]
    end
    return nothing
end

## Userdata structure for passing number of species as  a field in the structure, with 
## multithreading-aware pre-allocated buffers
Base.@kwdef struct MyDataPrealloc
    nspec::Int = 5
    npart::Int = 1
    DBinary::Symmetric{Float64, Matrix{Float64}} = Symmetric(fill(0.1, nspec, nspec))
    DKnudsen::Vector{Float64} = ones(nspec)
    diribc::Vector{Int} = [1, 2]
    M::Vector{DiffCache{Matrix{Float64}, Vector{Float64}}} = [DiffCache(ones(nspec, nspec)) for i in 1:npart]
    au::Vector{DiffCache{Vector{Float64}, Vector{Float64}}} = [DiffCache(ones(nspec)) for i in 1:npart]
    du::Vector{DiffCache{Vector{Float64}, Vector{Float64}}} = [DiffCache(ones(nspec)) for i in 1:npart]
    ipiv::Vector{Vector{Int}} = [zeros(Int, nspec) for i in 1:npart]
end
nspec(data::MyDataPrealloc)  = data.nspec


## Flux using pre-allocated buffers
function flux_diffcache(f, u, edge, data)
    T = eltype(u)
    n = data.nspec
    ipart = partition(edge)
    M = get_tmp(data.M[ipart], u)
    au = get_tmp(data.au[ipart], u)
    du = get_tmp(data.du[ipart], M)
    ipiv = data.ipiv[ipart]

    for ispec in 1:nspec(data)
        M[ispec, ispec] = 1.0 / data.DKnudsen[ispec]
        du[ispec] = u[ispec, 1] - u[ispec, 2]
        au[ispec] = 0.5 * (u[ispec, 1] + u[ispec, 2])
    end
    for ispec in 1:nspec(data)
        for jspec in 1:nspec(data)
            if ispec != jspec
                M[ispec, ispec] += au[jspec] / data.DBinary[ispec, jspec]
                M[ispec, jspec] = -au[ispec] / data.DBinary[ispec, jspec]
            end
        end
    end

    ## Pivoting linear system solution via RecursiveFactorizations.jl (see vfvm_functions.jl)
    inplace_linsolve!(M, du, ipiv)

    for ispec in 1:nspec(data)
        f[ispec] = du[ispec]
    end

    return nothing
end


function bcondition(f, u, node, data)
    for species in 1:nspec(data)
        boundary_dirichlet!(
            f, u, node; species, region = data.diribc[1],
            value = species % 2
        )
        boundary_dirichlet!(
            f, u, node; species, region = data.diribc[2],
            value = 1 - species % 2
        )
    end
    return nothing
end

function main(;
        n = 11, nspec = 5,
        dim = 2,
        Plotter = nothing,
        verbose = "",
        unknown_storage = :dense,
        flux = :flux_strided,
        strategy = nothing,
        assembly = :cellwise,
        npart = 1
    )
    h = 1.0 / convert(Float64, n - 1)
    X = collect(0.0:h:1.0)
    DBinary = Symmetric(fill(0.1, nspec, nspec))
    for ispec in 1:nspec
        DBinary[ispec, ispec] = 0
    end

    DKnudsen = fill(1.0, nspec)

    if dim == 1
        grid = simplexgrid(X)
        diribc = [1, 2]
    elseif dim == 2
        grid = simplexgrid(X, X)
        diribc = [4, 2]
    else
        grid = simplexgrid(X, X, X)
        diribc = [4, 2]
    end

    if npart > 1
        grid = partition(grid, PlainMetisPartitioning(; npart), nodes = true, edges = true)
    end

    function storage(f, u, node, data)
        f .= u
        return nothing
    end

    if flux == :flux_strided
        _flux = flux_strided
        data = MyDataStaticSizeInfo(nspec; DBinary, DKnudsen, diribc)
    elseif flux == :flux_diffcache
        _flux = flux_diffcache
        data = MyDataPrealloc(;nspec, npart = num_partitions(grid), DBinary, DKnudsen, diribc)
    else
        _flux = flux_marray
        data = MyDataStaticSizeInfo(nspec; DBinary, DKnudsen, diribc)
    end

    sys = VoronoiFVM.System(grid; flux = _flux, storage, bcondition, species = 1:nspec, data, assembly, unknown_storage)

    if !isnothing(strategy) && hasproperty(strategy, :precs)
        if isa(strategy.precs, BlockPreconBuilder)
            strategy.precs.partitioning = A -> partitioning(sys, Equationwise())
        end
        if isa(strategy.precs, ILUZeroPreconBuilder) && strategy.precs.blocksize != 1
            strategy.precs.blocksize = nspec
        end
    end
    control = SolverControl(method_linear = strategy)
    control.maxiters = 500
    u = solve(sys; verbose, control, log = true)
    return norm(u)
end

using Test
function runtests()
    strategies = [
        (method = UMFPACKFactorization(), dims = (1, 2, 3)),
        (method = KrylovJL_GMRES(precs = LinearSolvePreconBuilder(UMFPACKFactorization())), dims = (1, 2, 3)),
        (method = KrylovJL_GMRES(precs = BlockPreconBuilder(precs = LinearSolvePreconBuilder(UMFPACKFactorization()))), dims = (1, 2, 3)),
        (method = KrylovJL_GMRES(precs = BlockPreconBuilder(precs = AMGPreconBuilder())), dims = (2, 3)),
        (method = KrylovJL_BICGSTAB(precs = BlockPreconBuilder(precs = AMGPreconBuilder())), dims = (2,)),
        (method = KrylovJL_GMRES(precs = BlockPreconBuilder(precs = ILUZeroPreconBuilder())), dims = (2, 3)),
        (method = KrylovJL_GMRES(precs = ILUZeroPreconBuilder(blocksize = 5)), dims = (1, 2)),
    ]

    dimtestvals = [4.788926530387466, 15.883072449873742, 52.67819183426213]
    for dim in [1,2,3]
        for assembly in [:edgewise, :cellwise]
            for flux in [:flux_marray, :flux_strided, :flux_diffcache]
                for strategy in strategies
                    if dim in strategy.dims
                        result = main(; dim, assembly, flux, strategy = strategy.method) ≈ dimtestvals[dim]
                        if !result
                            @show dim, assembly, flux, strategy
                        end
                        @test result
                    end
                end
            end
        end
    end

    for dim in [2]
        for assembly in [:edgewise, :cellwise]
            for flux in [:flux_marray, :flux_strided, :flux_diffcache]
                result = main(; dim, n=100, assembly, flux, npart=20)
                @test  result ≈ 141.54097792523987
            end
        end
    end

    return nothing
end

end
