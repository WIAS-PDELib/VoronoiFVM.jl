#################################################################
"""
$(TYPEDEF)

Struct holding solution information for SparseSystem. Solution
is stored in a sparse matrix structure.

This class plays well with the abstract array interface.

Fields:

$(TYPEDFIELDS)
"""
mutable struct SparseSolutionArray{T, N, Ti} <: AbstractSolutionArray{T, N}
    """
    Sparse matrix holding actual data.
    """
    u::SparseMatrixCSC{T, Ti}

    history::Union{NewtonSolverHistory, Nothing}
end

"""
    $(TYPEDSIGNATURES)

`SparseSolutionArray` constructor
"""
SparseSolutionArray(a::SparseMatrixCSC{Tv, Ti}) where {Tv, Ti} = SparseSolutionArray{Tv, 2, Ti}(a, nothing)

"""
    solutionarray(a::SparseMatrixCSC)
"""
solutionarray(a::SparseMatrixCSC{Tv, Ti}) where {Tv, Ti} = SparseSolutionArray{Tv, 2, Ti}(a, nothing)

"""
$(TYPEDEF)

Wrapper structure to access sparse solution indices.
"""
struct SparseSolutionIndices
    a::SparseSolutionArray
end

"""
$(SIGNATURES)

Return indices for sparse solution array.
"""
unknown_indices(a::SparseSolutionArray) = SparseSolutionIndices(a)

Base.getindex(idx::SparseSolutionIndices, i, j) = dof(idx.a, i, j)

##################################################################
"""
$(SIGNATURES)

Vector of degrees of freedom in sparse solution array.
"""
dofs(a::SparseSolutionArray) = a.u.nzval
dofs(a::SparseMatrixCSC) = a.nzval

Base.size(a::SparseSolutionArray) = size(a.u)
##################################################################
"""
$(SIGNATURES)
    
Create a copy of sparse solution array
"""
function Base.copy(this::SparseSolutionArray{T, N, Ti}) where {T, N, Ti}
    return SparseSolutionArray{T, N, Ti}(
        SparseMatrixCSC(
            this.u.m,
            this.u.n,
            this.u.colptr,
            this.u.rowval,
            Base.copy(this.u.nzval)
        ), this.history
    )
end


"""
$(SIGNATURES)
    
Create a similar uninitialized sparse solution array
"""
function Base.similar(this::SparseSolutionArray{T, N, Ti}) where {T, N, Ti}
    return SparseSolutionArray{T, N, Ti}(
        SparseMatrixCSC(
            this.u.m,
            this.u.n,
            this.u.colptr,
            this.u.rowval,
            Base.similar(this.u.nzval)
        ), nothing
    )
end
##################################################################
"""
$(SIGNATURES)

Get number of degree of freedom. Return 0 if species is not defined in node.
"""
function dof(a::SparseSolutionArray, i, j)
    A = a.u
    coljfirstk = Int(A.colptr[j])
    coljlastk = Int(A.colptr[j + 1] - 1)
    searchk = searchsortedfirst(A.rowval, i, coljfirstk, coljlastk, Base.Order.Forward)
    if searchk <= coljlastk && A.rowval[searchk] == i
        return searchk
    end
    return 0
end

##################################################################
"""
$(SIGNATURES)

Set value for degree of freedom.
"""
function setdof!(a::SparseSolutionArray, v, i::Integer)
    return a.u.nzval[i] = v
end

##################################################################
"""
$(SIGNATURES)

Return  value for degree of freedom.
"""
getdof(a::SparseSolutionArray, i::Integer) = a.u.nzval[i]

Base.:-(a::SparseSolutionArray, b::SparseSolutionArray) = SparseSolutionArray(a.u - b.u)
Base.:+(a::SparseSolutionArray, b::SparseSolutionArray) = SparseSolutionArray(a.u + b.u)

##################################################################
"""
$(SIGNATURES)

Accessor for sparse solution array.
"""
function Base.setindex!(a::SparseSolutionArray, v, ispec::Int, inode::Int)
    searchk = dof(a, ispec, inode)
    if searchk > 0
        setdof!(a, v, searchk)
        return a
    end
    # TODO: what is the right reacton here ?
    # Ignoring seems to be better, so we can broadcast etc.
    # throw(DomainError("undefined degree of freedom"))
end

##################################################################
"""
$(SIGNATURES)

Accessor for sparse solution array.
"""
function Base.getindex(a::SparseSolutionArray, ispec::Int, inode::Int)
    searchk = dof(a, ispec, inode)
    if searchk > 0
        return getdof(a, searchk)
    end
    #
    # TODO: what is the right reacton here ?
    # Actually, NaN plays well with pyplot...
    return NaN
end

"""
$(TYPEDSIGNATURES)

Add residual value into global degree of freedom

(internal)
"""
_add(U::SparseSolutionArray, idof, val) = U.u.nzval[idof] += val

"""
$(TYPEDSIGNATURES)

Set residual value for global degree of freedom

(internal)
"""
_set(U::SparseSolutionArray, idof, val) = U.u.nzval[idof] = val
