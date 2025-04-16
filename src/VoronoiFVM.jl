"""
    VoronoiFVM

$(read(joinpath(@__DIR__, "..", "README.md"), String))
"""
module VoronoiFVM

using BandedMatrices: BandedMatrices, BandedMatrix, Zeros
import Colors
using CommonSolve: CommonSolve, solve, solve!
using DiffResults: DiffResults
using DocStringExtensions: DocStringExtensions, SIGNATURES, TYPEDEF,
    TYPEDFIELDS, TYPEDSIGNATURES
using ExtendableGrids: ExtendableGrids, BEdgeNodes, BFaceCells, BFaceEdges,
    BFaceGeometries, BFaceNodes, BFaceNormals, BFaceRegions,
    Cartesian1D, Cartesian2D, Cartesian3D, CellEdges,
    CellGeometries, CellNodes, CellRegions,
    CoordinateSystem, Coordinates, Cylindrical2D, Edge1D,
    EdgeCells, EdgeNodes, ExtendableGrid,
    Polar1D, Spherical1D, Tetrahedron3D,
    Triangle2D, Vertex0D, VoronoiFaceCenters, coord_type,
    dim_space, index_type, local_celledgenodes, num_bfaces,
    num_cells, num_edges, num_nodes, num_cellregions, num_bfaceregions, num_targets,
    simplexgrid, subgrid, tricircumcenter!,
    num_partitions, pcolor_partitions, pcolors, num_pcolors,
    PColorPartitions, PartitionCells, PartitionBFaces, PartitionNodes, PartitionEdges

using DifferentiationInterface: DifferentiationInterface, AutoSparse, AutoForwardDiff, prepare_jacobian

using ExtendableSparse: ExtendableSparse, BlockPreconditioner,
    ExtendableSparseMatrix,
    ExtendableSparseMatrixCSC,
    MTExtendableSparseMatrixCSC,
    AbstractExtendableSparseMatrixCSC,
    PointBlockILUZeroPreconditioner, factorize!, flush!,
    nnz, rawupdateindex!, updateindex!, nnznew

using ForwardDiff: ForwardDiff, value
using GridVisualize: GridVisualize, GridVisualizer
using InteractiveUtils: InteractiveUtils
using JLD2: JLD2, jldopen
using LinearAlgebra: LinearAlgebra, Diagonal, I, Tridiagonal, isdiag, ldiv!, norm
using LinearSolve: LinearSolve, KrylovJL_BICGSTAB,
    KrylovJL_CG, KrylovJL_GMRES, LinearProblem,
    SparspakFactorization, UMFPACKFactorization, init, reinit!
using Printf: Printf, @printf, @sprintf
using Random: Random, AbstractRNG
using RecursiveArrayTools: RecursiveArrayTools, AbstractDiffEqArray, DiffEqArray
import RecursiveFactorization
using SciMLBase: SciMLBase
using SparseArrays: SparseArrays, SparseMatrixCSC, dropzeros!, nonzeros,
    nzrange, spzeros, issparse
using SparseConnectivityTracer: SparseConnectivityTracer, TracerSparsityDetector
using SparseMatrixColorings: GreedyColoringAlgorithm, sparsity_pattern
using StaticArrays: StaticArrays, @MVector, @SArray, @SMatrix
using Statistics: Statistics, mean
using TextWrap: print_wrapped

"""
   $(TYPEDEF)

Abstract type for geometry items (node,bnode,edge, bedge)
"""
abstract type AbstractGeometryItem{Tc <: Number, Tp <: Number, Ti <: Integer} end
export AbstractGeometryItem

"""
   $(TYPEDEF)

Abstract type for stationary solution. Subtype of `AbstractArray`.
"""
abstract type AbstractSolutionArray{T, N} <: AbstractArray{T, N} end
Base.getindex(a::AbstractSolutionArray, i::Int, j::Int) = getindex(a.u, i, j)
Base.setindex!(a::AbstractSolutionArray, v, i::Int, j::Int) = setindex!(a.u, v, i, j)
Base.size(a::AbstractSolutionArray) = size(a.u)
solutionarray(a::AbstractSolutionArray) = a

export AbstractSolutionArray

include("vfvm_physics.jl")
# see https://discourse.julialang.org/t/is-compat-jl-worth-it-for-the-public-keyword/119041/34
VERSION >= v"1.11.0-DEV.469" && eval(Meta.parse("public Physics, AbstractPhysics, AbstractData"))

include("vfvm_functions.jl")
export fbernoulli
export fbernoulli_pm
export inplace_linsolve!


include("vfvm_history.jl")
export NewtonSolverHistory, TransientSolverHistory, details

include("vfvm_densesolution.jl")
include("vfvm_sparsesolution.jl")
VERSION >= v"1.11.0-DEV.469" && eval(Meta.parse("public dofs"))
export num_dof
export dof
export getdof
export setdof!
export unknown_indices, SparseSolutionIndices

include("vfvm_transientsolution.jl")
export TransientSolution

include("vfvm_xgrid.jl")
export cartesian!, circular_symmetric!, spherical_symmetric!
export coordinates

"""
$(TYPEDEF)

Abstract type for finite volume system structure.
"""
abstract type AbstractSystem{Tv <: Number, Tc <: Number, Ti <: Integer, Tm <: Integer} end
include("vfvm_geometryitems.jl")
include("vfvm_assemblydata.jl")
include("vfvm_system.jl")
include("vfvm_state.jl")
export unknowns
export num_species
export enable_species!
export enable_boundary_species!
export update_grid!
export boundary_dirichlet!
export boundary_neumann!
export boundary_robin!
export ramp
export value
export physics!
export history, history_summary, history_details
export evaluate_residual_and_jacobian
export edgelength
export viewK, viewL, data
export hasoutflownode, isoutflownode, outflownode
export parameters

VERSION >= v"1.11.0-DEV.469" && eval(Meta.parse("public System, AbstractSystem, SystemState"))

# export to be deprecated
export partitioning, Equationwise

include("vfvm_logging_exceptions.jl")
VERSION >= v"1.11.0-DEV.469" && eval(Meta.parse("public print_output!, log_output!"))

include("vfvm_formfactors.jl")
export meas, project
export unknown_indices
export edgevelocities, bfacevelocities, bfacenodefactors
export time, region, embedparam
export calc_divergences

include("vfvm_solvercontrol.jl")
export fixed_timesteps!, NewtonControl, SolverControl
include("vfvm_linsolve_deprecated.jl")
include("vfvm_linsolve.jl")
export DirectSolver, GMRESIteration, CGIteration, BICGstabIteration, NoBlock, EquationBlock, PointBlock

include("vfvm_assembly.jl")
include("vfvm_solver.jl")
export solve!, solve

include("vfvm_postprocess.jl")
export nodeflux
export integrate
export l2norm, lpnorm
export w1pseminorm, h1seminorm
export w1pnorm, h1norm
export lpw1pnorm, l2h1norm
export lpw1pseminorm, l2h1seminorm
export nodevolumes
export nondelaunay

include("vfvm_testfunctions.jl")
export testfunction
export TestFunctionFactory

include("vfvm_quantities.jl")
export ContinuousQuantity
export DiscontinuousQuantity
export InterfaceQuantity
export subgrids, views

include("vfvm_impedance.jl")
export impedance, freqdomain_impedance
export measurement_derivative

include("vfvm_diffeq_interface.jl")
export eval_rhs!, eval_jacobian!, mass_matrix, prepare_diffeq!

include("gridvisualize.jl")
export plothistory
#include("precompile.jl")

end
