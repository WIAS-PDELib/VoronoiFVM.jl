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

log_output::Bool=isdefined(Main,:PlutoRunner)

"""
    print_output!()

Write all subsequent output of the VoronoiFVM package to screen via println().
This is the default if the code is run outside of Pluto notebooks.

If enabled in Pluto notebooks, the output is directed to the terminal  widget
below a pluto cell. Warnings in addition are sent via logging in order to not
miss them.

This behavior can be changed via [`log_output!`](@ref).

For fine-tuning  solver output, see the `verbose` flag in [`SolverControl`](@ref).
"""
print_output!() = global log_output = false


"""
    log_output!()

Write all subsequent output of the VoronoiFVM package to screen via Julia's 
logging methods (using @info or @warn)

This is the default if the code is run in Pluto notebooks.

This behavior can be changed via [`print_output!`](@ref).

For fine-tuning  solver output, see the `verbose` flag in [`SolverControl`](@ref).
"""
log_output!()  = global log_output = true


function _warn(str)
    global log_output
    if log_output
        @warn str
    else
        if isdefined(Main,:PlutoRunner)
            @warn str
        end
        println("WARNING: $(str)")
    end
    return
end

function _info(str)
    global log_output
    if log_output
        @info str
    else
        println(str)
    end
    return
end




"""
Print error when catching exceptions
"""
function _print_error(err, st)
    io = IOBuffer()
    println(io, err)
    nlines = 0
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
    # if length(st) > nlines
    #     println(io, "...")
    # end
    return _warn(String(take!(io)))
end
