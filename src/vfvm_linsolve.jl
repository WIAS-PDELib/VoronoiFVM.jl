################################################################
# These are needed to enable iterative solvers to work with dual numbers
# TODO: Remove this Pirate's nest
Base.Float64(x::ForwardDiff.Dual) = value(x)
function Random.rand(
        rng::AbstractRNG,
        ::Random.SamplerType{ForwardDiff.Dual{T, V, N}}
    ) where {T, V, N}
    return ForwardDiff.Dual{T, V, N}(rand(rng, V))
end

# TODO: these may be not anymore needed
canonical_matrix(A) = A
canonical_matrix(A::AbstractExtendableSparseMatrixCSC) = SparseMatrixCSC(A)

function _solve_linear!(u, state, nlhistory, control, method_linear, A, b)
    if isnothing(state.linear_cache)
        if !isa(method_linear, LinearSolve.SciMLLinearSolveAlgorithm)
            @warn "use of $(typeof(method_linear)) is deprecated, use an algorithm from LinearSolve"
        end
        if hasproperty(method_linear, :precs) && !isnothing(method_linear.precs)
            Pl = nothing
        else
            Pl = control.precon_linear(canonical_matrix(A))
            if !isa(Pl, Identity) && isa(method_linear, LinearSolve.AbstractKrylovSubspaceMethod)
                @warn "Use of control.precon_linear is deprecated. Use the `precs` API of LinearSolve"
            end
        end
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
        if hasproperty(method_linear, :precs) && !isnothing(method_linear.precs)
            reinit!(state.linear_cache; A = canonical_matrix(A), b, reuse_precs = !control.keepcurrent_linear)
            if control.keepcurrent_linear
                nlhistory.nlu += 1
            end
        else
            state.linear_cache.A = canonical_matrix(A)
            state.linear_cache.b = b
            if control.keepcurrent_linear
                nlhistory.nlu += 1
                state.linear_cache.Pl = control.precon_linear(canonical_matrix(A))
            end
        end
    end

    return try
        sol = LinearSolve.solve!(state.linear_cache)
        u .= sol.u
        nliniter = sol.iters
        nlhistory.nlin = sol.iters
    catch err
        if (control.handle_exceptions)
            _print_error(err, stacktrace(catch_backtrace()))
            throw(LinearSolverError())
        else
            rethrow(err)
        end
    end
end
