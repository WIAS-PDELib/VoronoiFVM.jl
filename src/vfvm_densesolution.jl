"""
    $(TYPEDEF)

Dense storage of solution. Subtype of [`AbstractSolutionArray`](@ref)


Fields:

$(TYPEDFIELDS)
"""
mutable struct DenseSolutionArray{T, N} <: AbstractSolutionArray{T, N}
    u::Array{T, N}
    history::Union{NewtonSolverHistory, Nothing}
end


"""
 $(TYPEDSIGNATURES)
    
`DenseSolutionArray` constructor.
"""
DenseSolutionArray(u::Matrix{T}) where {T} = DenseSolutionArray{T, 2}(u, nothing)

"""
    solutionarray(a::Matrix)
"""
solutionarray(u::Matrix{T}) where {T} = DenseSolutionArray{T, 2}(u, nothing)

"""
 $(TYPEDSIGNATURES)
    
`DenseSolutionArray` constructor.
"""
DenseSolutionArray{T, 2}(nspec::Int, nnodes::Int) where {T} = DenseSolutionArray{T, 2}(Matrix{T}(nspec, nnodes), nothing)

"""
 $(TYPEDSIGNATURES)
    
`DenseSolutionArray` constructor.
"""
DenseSolutionArray{T, 2}(::UndefInitializer, nspec::Int, nnodes::Int) where {T} = DenseSolutionArray{T, 2}(Matrix{T}(undef, nspec, nnodes), nothing)


"""
 $(SIGNATURES)
    
Get degree of freedom number
"""
dof(a::DenseSolutionArray, ispec, K) = (K - 1) * size(a, 1) + ispec

"""
$(SIGNATURES)

Return indices for dense solution array.
"""
unknown_indices(a::DenseSolutionArray) = LinearIndices(a)


Base.copy(a::DenseSolutionArray) = DenseSolutionArray(copy(a.u))

Base.similar(a::DenseSolutionArray) = DenseSolutionArray(similar(a.u))

"""
$(SIGNATURES)

Vector of degrees of freedom in solution array.
"""
dofs(a::DenseSolutionArray) = vec(a.u)
dofs(a::Array) = vec(a)

"""
$(TYPEDSIGNATURES)

Add residual value into global degree of freedom

(Internal method)
"""
_add(U::DenseSolutionArray, idof, val) = U[CartesianIndices(U)[idof]] += val

"""
$(TYPEDSIGNATURES)

Set residual value for global degree of freedom

(Internal method)
"""
_set(U::DenseSolutionArray, idof, val) = U[CartesianIndices(U)[idof]] = val

##################################################################
