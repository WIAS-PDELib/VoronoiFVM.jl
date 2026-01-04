# TODO: these may be not anymore needed but require
# that preconditioners work with AbstractSparseMatrixCSC.
canonical_matrix(A) = A
canonical_matrix(A::AbstractExtendableSparseMatrixCSC) = SparseMatrixCSC(A)

function _solve_linear!(u, state, nlhistory, control, method_linear, A, b, niter)
    if isnothing(state.linear_cache)
        if !isa(method_linear, LinearSolve.SciMLLinearSolveAlgorithm)
            @warn "use of $(typeof(method_linear)) is deprecated, use an algorithm from LinearSolve"
        end
        Pl = nothing
        nlhistory.nlu += 1
        p = LinearProblem(canonical_matrix(A), b)
        state.linear_cache = init(
            p,
            method_linear;
            abstol = control.abstol_linear,
            reltol = control.reltol_linear,
            maxiters = control.maxiters_linear,
            verbose = doprint(control, 'l'),
            Pl,
        )
    else
        reuse_precs = !control.keepcurrent_linear && niter > 1
        reinit!(state.linear_cache; A = canonical_matrix(A), b, reuse_precs)
        if !reuse_precs
            nlhistory.nlu += 1
        end
    end

    try
        sol = LinearSolve.solve!(state.linear_cache)
        u .= sol.u
        nliniter = sol.iters
        nlhistory.nlin = sol.iters
    catch err
        if (control.handle_exceptions)
            _warn(err, stacktrace(catch_backtrace()))
            throw(LinearSolverError())
        else
            rethrow(err)
        end
    end
    return nothing
end
