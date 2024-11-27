"""
    $(TYPEDEF)

History information for one Newton solve of a nonlinear system.
As an abstract vector it provides the history of the update norms.
See [`summary`](@ref) and [`details`](@ref) for other ways to extract information.

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct NewtonSolverHistory <: AbstractVector{Float64}
    """ number of Jacobi matrix factorizations """
    nlu::Int = 0

    """ number of linear iteration steps / factorization solves"""
    nlin::Int = 0

    """ Elapsed time for solution """
    time::Float64 = 0.0

    """ Elapsed time for assembly """
    tasm::Float64 = 0.0

    """ Elapsed time for linear solve """
    tlinsolve::Float64 = 0.0

    """ History of norms of ``||u_{i+1}-u_i||``"""
    updatenorm = zeros(0)

    """ History of norms of ``|\\;||u_{i+1}||_1 - ||u_{i}||_1\\;|/ ||u_{i}||_1`` """
    l1normdiff = zeros(0)
end

Base.getindex(h::NewtonSolverHistory, i) = h.updatenorm[i]
Base.size(h::NewtonSolverHistory) = size(h.updatenorm)

"""
    summary(h::NewtonSolverHistory)

Return named tuple summarizing history.
"""
function Base.summary(h::NewtonSolverHistory)
    return (
        seconds = round(h.time; sigdigits = 3),
        tasm = round(h.tasm; sigdigits = 3),
        tlinsolve = round(h.tlinsolve; sigdigits = 3),
        iters = length(h.updatenorm),
        absnorm = round(h.updatenorm[end]; sigdigits = 3),
        relnorm = round(h.updatenorm[end] / h.updatenorm[1]; sigdigits = 3),
        roundoff = round(h.l1normdiff[end]; sigdigits = 3),
        factorizations = h.nlu,
        liniters = h.nlin,
    )
end

"""
    details(h::NewtonSolverHistory)

Return array of named tuples  with info on each iteration step
"""
function details(h::NewtonSolverHistory)
    a = []
    for i in 1:length(h)
        push!(
            a,
            (
                update = round(h[i]; sigdigits = 3), contraction = round(i > 1 ? h[i] / h[i - 1] : 1.0; sigdigits = 3),
                round = round(h.l1normdiff[i]; sigdigits = 3),
            )
        )
    end
    return a
end

#################################################################################

"""
    $(TYPEDEF)

History information for transient solution/parameter embedding

As an abstract vector it provides the histories of each implicit Euler/embedding step.
See [`summary`](@ref) and [`details`](@ref) for other ways to extract information.


$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct TransientSolverHistory <: AbstractVector{NewtonSolverHistory}
    """ Histories of each implicit Euler Newton iteration """
    histories = Vector{NewtonSolverHistory}(undef, 0)

    """ Time values """
    times = zeros(0)

    """ Update norms used for step control"""
    updates = zeros(0)
end

Base.getindex(h::TransientSolverHistory, i) = h.histories[i]
Base.size(h::TransientSolverHistory) = size(h.histories)
Base.push!(h::TransientSolverHistory, hx) = push!(h.histories, hx)

"""
    summary(h::TransientSolverHistory)

Return named tuple summarizing history.
"""
function Base.summary(hh::TransientSolverHistory)
    hx = view(hh, 2:length(hh))
    return (
        seconds = round(sum(h -> h.time, hx); sigdigits = 3),
        tasm = round(sum(h -> h.tasm, hx); sigdigits = 3),
        tlinsolve = round(sum(h -> h.tlinsolve, hx); sigdigits = 3),
        steps = length(hh.histories),
        iters = sum(h -> length(h.updatenorm), hx),
        maxabsnorm = round(maximum(h -> h.updatenorm[end], hx); sigdigits = 3),
        maxrelnorm = round(maximum(h -> h.updatenorm[end] / h.updatenorm[1], hx); sigdigits = 3),
        maxroundoff = round(maximum(h -> h.l1normdiff[end], hx); sigdigits = 3),
        iters_per_step = round(mean(h -> length(h.updatenorm), hx); sigdigits = 3),
        facts_per_step = round(mean(h -> h.nlu, hx); sigdigits = 3),
        liniters_per_step = round(mean(h -> h.nlin, hx); sigdigits = 3),
    )
end

"""
    details(h::TransientSolverHistory)

Return array of details of each solver step
"""
function details(hh::TransientSolverHistory)
    a = []
    for i in 1:length(hh)
        push!(
            a,
            (
                t = round(hh.times[i]; sigdigits = 3), Δu = round(hh.updates[i]; sigdigits = 3), summary = summary(hh[i]),
                detail = details(hh[i]),
            )
        )
    end
    return a
end


"""
    $(TYPEDEF)

History information for DiffEqHistory
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct DiffEqHistory
    """ number of combined jacobi/rhs evaluations"""
    nd = 0
    """ number of combined jacobi evaluations"""
    njac = 0
    """ number of rhs evaluations"""
    nf = 0
end

details(h::DiffEqHistory) = (nd = h.nd, njac = h.njac, nf = h.nf)
Base.summary(h::DiffEqHistory) = details(h)

"""
    history(sol)

Return solver history if `log` was set to true.
See  see [`NewtonSolverHistory`](@ref), [`TransientSolverHistory`](@ref).
"""
history(sol) = sol.history

"""
    history_details(sol)

Return details of solver history from last `solve` call, if `log` was set to true.
See [`details`](@ref).
"""
history_details(sol) = details(sol.history)

"""
    history_summary(sol)

Return summary of solver history from last `solve` call, if `log` was set to true.
"""
history_summary(sol) = summary(sol.history)
