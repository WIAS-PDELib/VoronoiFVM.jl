"""
$(TYPEDEF)

Exception thrown if Newton's method convergence fails.
"""
struct ConvergenceError <: Exception end

"""
$(TYPEDEF)

Exception thrown if error occurred during assembly (e.g. domain error)
"""
struct AssemblyError <: Exception end

"""
$(TYPEDEF)

Exception thrown if error occurred during factorization.
"""
struct LinearSolverError <: Exception end

"""
$(TYPEDEF)

Exception thrown if embedding fails
"""
struct EmbeddingError <: Exception
    msg::String
end

"""
Print error when catching exceptions
"""
function _print_error(err, st)
    io = IOBuffer()
    println(io)
    println(io, err)
    nlines = 5
    for i in 1:min(nlines, length(st))
        line = @sprintf("%s", st[i])
        L = length(line)
        if L < 80
            println(io, line)
        else
            print(io, line[1:35])
            print(io, " ... ")
            println(io, line[(L - 35):L])
        end
    end
    if length(st) > nlines
        println(io, "...")
    end
    println(io)
    return @warn String(take!(io))
end
