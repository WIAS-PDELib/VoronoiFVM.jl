"""
$(TYPEDEF)

Abstract type for transient solution
"""
abstract type AbstractTransientSolution{T, N, A, B} <: AbstractDiffEqArray{T, N, A} end

"""
$(TYPEDEF)

Transient solution structure

## Fields
$(TYPEDFIELDS)
    

## Interface

Object of this type adhere to the `AbstractDiffEqArray`  interface.
For indexing and interpolation, see [https://diffeq.sciml.ai/stable/basics/solution/](https://diffeq.sciml.ai/stable/basics/solution/).

In particular, a TransientSolution `sol` can be accessed as follows:
- `sol[i]` contains the solution for timestep `i`
- `sol[ispec,:,i]` contains the solution for component `ispec` at timestep `i`
- `sol(t)` returns a (linearly) interpolated solution value for `t`.
- `sol.t[i]` is the corresponding time for timestep `i`
- `sol[ispec,ix,i]` refers to solution of component `ispec` at node `ix` at moment `i`


"""
mutable struct TransientSolution{T, N, A, B} <: AbstractTransientSolution{T, N, A, B}
    """
    Vector of solutions
    """
    u::A

    """
    Vector of times
    """
    t::B

    """
    History
    """
    history::Union{TransientSolverHistory, DiffEqHistory}
end


function TransientSolution(vec::AbstractVector{T}, ts, ::NTuple{N}) where {T, N}
    return TransientSolution{eltype(T), N, typeof(vec), typeof(ts)}(vec, ts, TransientSolverHistory())
end

TransientSolution(vec::AbstractVector, ts::AbstractVector) = TransientSolution(vec, ts, (size(vec[1])..., length(vec)))

@doc """
    append!(tsol::AbstractTransientSolution, t, sol)

Append time step solution `sol` as solution at time `t` to tsol. 
`sol` will be directly references in `tsol`. This method does not copy `sol`.
If during a time steping process it is the same vector, a `copy` should appended.

Defined in VoronoiFVM.jl.
"""
Base.append!

Base.append!(s::AbstractTransientSolution, t::Real, sol::AbstractArray) = push!(s.t, t), push!(s.u, sol)

Base.append!(s::AbstractTransientSolution, t::Real, sol::AbstractSolutionArray) = append!(s, t, sol.u)

(sol::AbstractTransientSolution)(t) = _interpolate(sol, t)

function _interpolate(sol, t)
    if isapprox(t, sol.t[1]; atol = 1.0e-10 * abs(sol.t[2] - sol.t[1]))
        return sol[1]
    end
    idx = searchsortedfirst(sol.t, t)
    if idx == 1 || idx > length(sol)
        return nothing
    end
    if t == sol.t[idx - 1]
        return sol[idx - 1]
    else
        retval = similar(sol.u[idx])
        dt = sol.t[idx] - sol.t[idx - 1]
        a = (sol.t[idx] - t) / dt
        b = (t - sol.t[idx - 1]) / dt
        retval .= a * sol.u[idx - 1] + b * sol.u[idx]
    end
end

mutable struct VectorOfDiskArrays{T} <: AbstractVector{T}
    fname::String
    file::Union{JLD2.JLDFile, Nothing}
    n::Int64
end

function Base.push!(v::VectorOfDiskArrays, obj)
    v.n += 1
    return if isnothing(v.file)
        jldopen(v.fname, "a+") do file
            file[string(v.n)] = obj
        end
    else
        v.file[string(v.n)] = obj
    end
end

Base.size(v::VectorOfDiskArrays) = (v.n,)
Base.length(v::VectorOfDiskArrays) = v.n
Base.eltype(v::VectorOfDiskArrays{T}) where {T} = T
function Base.getindex(v::VectorOfDiskArrays, i)
    return if isnothing(v.file)
        jldopen(v.fname, "r") do file
            file[string(i)]
        end
    else
        v.file[string(i)]
    end
end

_tempname() = Base.VERSION < v"1.4" ? tempname() : tempname(pwd()) * ".jld2"

"""
````
VectorOfDiskArrays(firstobj:AbstractArray;
                   keep_open=true,
                   fname= fname=tempname(pwd())*".jld2")
````
Constructor of vector of arrays stored on disk (via JLD2).

- `keep_open`: if true, disk file is not closed during the existence of the object
- `fname`: file name for the disk file


The disk file is automatically removed if the object is garbage collected.
"""
function VectorOfDiskArrays(obj::AbstractArray{T}; keep_open = true, fname = _tempname()) where {T}
    file = jldopen(fname, "a+")
    file[string(1)] = obj
    if !keep_open
        close(file)
        file = nothing
    end
    v = VectorOfDiskArrays{T}(fname, file, 1)
    finalizer(v -> (isnothing(v.file) ? nothing : close(v.file); rm(v.fname; force = true)), v)
    return v
end

function Base.append!(s::AbstractTransientSolution{T, N, VectorOfDiskArrays{T}, B}, t::Real, sol::AbstractArray) where {T, N, B}
    return push!(s.t, t), push!(s.u, sol)
end

"""
````
TransientSolution(t0,inival;
                  in_memory=true,
                  keep_open=true,
                  fname=tempname(pwd())*".jld2"
````
Constructor of transient solution with initial value and initial time.

- `in_memory`: if true (default), data are kept in main memory, otherwise on disk (via JLD2)
- `keep_open`: if true, disk file is not closed during the existence of the object
- `fname`: file name for the disk file
"""
function TransientSolution(
        t0::Number,
        inival::AbstractArray{T};
        in_memory = true,
        keep_open = true,
        fname = _tempname()
    ) where {T}
    return if !in_memory && !isa(inival, SparseSolutionArray)
        TransientSolution(VectorOfDiskArrays(inival; keep_open = keep_open, fname = fname), [t0])
    else
        TransientSolution([inival], [t0])
    end
end

TransientSolution(t0::Number, inival::AbstractSolutionArray{T, N}; kwargs...) where {T, N} = TransientSolution(t0, inival.u; kwargs...)
