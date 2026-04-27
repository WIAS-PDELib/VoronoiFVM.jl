# TODO: these may be not anymore needed but require
# that preconditioners work with AbstractSparseMatrixCSC.
canonical_matrix(A) = A
canonical_matrix(A::AbstractExtendableSparseMatrixCSC) = SparseMatrixCSC(A)

function _solve_linear!(u, state, nlhistory, control, method_linear, A, b, reuse_precs)

    if isnothing(state.linear_cache)
        if !isa(method_linear, LinearSolve.SciMLLinearSolveAlgorithm)
            @warn "use of $(typeof(method_linear)) is deprecated, use an algorithm from LinearSolve"
        end
        nlhistory.nlu += 1
        nlhistory.tlinsolve_setup += @elapsed begin
            Pl = nothing
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
        end
        if doprint(control, 'l')
            out = @sprintf("    [l]inear: factorize #%d\n", nlhistory.nlu)
            _info(out)
        end
    else
        nlhistory.tlinsolve_setup += @elapsed begin
            reinit!(state.linear_cache; A = canonical_matrix(A), b, reuse_precs)
        end
        if !reuse_precs
            nlhistory.nlu += 1
            if doprint(control, 'l')
                out = @sprintf("    [l]inear: factorize #%d", nlhistory.nlu)
                _info(out)
            end
        end
    end

    try
        local sol
        nlhistory.tlinsolve_solve += @elapsed begin
            sol = LinearSolve.solve!(state.linear_cache)
        end
        u .= sol.u
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
