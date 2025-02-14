#
# See
# https://discourse.julialang.org/t/is-it-possible-to-detect-if-julia-is-ahead-of-time-precompiling/78631
#
is_precompiling() = ccall(:jl_generating_output, Cint, ()) == 1

"""
$(SIGNATURES)

Solve time step problem. This is the core routine
for implicit Euler and stationary solve.
"""
function solve_step!(
        state,
        solution, # old time step solution resp. initial value
        oldsol, # old time step solution resp. initial value
        control,
        time,
        tstep,
        embedparam,
        params
    )
    nlhistory = NewtonSolverHistory()
    tasm = 0.0
    tlinsolve = 0.0
    t = @elapsed begin
        solution .= oldsol
        residual = state.residual
        update = state.update
        _initialize!(solution, state.system, state.data; time, λ = embedparam, params)

        Tv = eltype(solution)
        method_linear = state.system.matrixtype == :sparse ? control.method_linear : nothing
        if isnothing(method_linear) && state.system.matrixtype == :sparse
            if Tv != Float64
                method_linear = SparspakFactorization()
            else
                method_linear = UMFPACKFactorization() # seems to do the best pivoting
            end
        end

        oldnorm = 1.0
        converged = false
        damp = 1.0
        if !state.system.is_linear
            if doprint(control, 'n')
                infoout("  [n]ewton: #it(lin)  |update| cont3tion   |round| #rd")
            end
            damp = control.damp_initial
            rnorm = control.rnorm(solution)
        end

        nlu_reuse = 0
        nround = 0
        tolx = 0.0
        ncalloc = 0
        nballoc = 0
        neval = 0
        niter = 1

        while niter <= control.maxiters
            # Create Jacobi matrix and RHS for Newton iteration
            try
                tasm += @elapsed nca, nba, nev = eval_and_assemble(
                    state.system,
                    solution,
                    oldsol,
                    residual,
                    state.matrix,
                    state.dudp,
                    time,
                    tstep,
                    embedparam,
                    state.data,
                    params;
                    edge_cutoff = control.edge_cutoff,
                )
                ncalloc += nca
                nballoc += nba
                neval += nev
            catch err
                if (control.handle_exceptions)
                    _print_error(err, stacktrace(catch_backtrace()))
                    throw(AssemblyError())
                else
                    rethrow(err)
                end
            end

            tlinsolve += @elapsed _solve_linear!(
                dofs(update),
                state,
                nlhistory,
                control,
                method_linear,
                state.matrix,
                dofs(residual)
            )

            dofs(solution) .-= damp * dofs(update)

            if state.system.is_linear
                converged = true
                break
            end

            damp = min(damp * control.damp_growth, 1.0)
            norm = control.unorm(update)
            if tolx == 0.0
                tolx = norm * control.reltol
            end
            dnorm = 1.0
            rnorm_new = control.rnorm(solution)
            if rnorm > 1.0e-50
                dnorm = abs((rnorm - rnorm_new) / rnorm)
            end

            if dnorm < control.tol_round
                nround = nround + 1
            else
                nround = 0
            end

            if control.log
                push!(nlhistory.l1normdiff, dnorm)
                push!(nlhistory.updatenorm, norm)
            end
            if doprint(control, 'n')
                if control.reltol_linear < 1.0
                    itstring = @sprintf("  [n]ewton: % 3d(% 3d)", niter, nlhistory.nlin)
                else
                    itstring = @sprintf("it=% 3d", niter)
                end
                if control.max_round > 0
                   out=@sprintf(
                        "%s %.3e %.3e %.3e % 2d",
                        itstring,
                        norm,
                        norm / oldnorm,
                        dnorm,
                        nround
                    )
                    infoout(out)
                else
                    out= @sprintf("%s %.3e %.3e", itstring, norm, norm / oldnorm)
                    infoout(out)
                end
            end
            if niter > 1 && norm / oldnorm > 1.0 / control.tol_mono
                converged = false
                break
            end

            if norm < control.abstol || norm < tolx
                converged = true
                break
            end
            oldnorm = norm
            rnorm = rnorm_new

            if nround > control.max_round
                converged = true
                break
            end
            niter = niter + 1
        end
        if !converged
            throw(ConvergenceError())
        end
    end
    if control.log
        nlhistory.time = t
        nlhistory.tlinsolve = tlinsolve
        nlhistory.tasm = tasm
    end

    if neval > 0
        if ncalloc ÷ neval + nballoc ÷ neval > 0 && doprint(control, 'a') && !is_precompiling()
            warnout("[a]llocations in assembly loop: cells: $(ncalloc ÷ neval), bfaces: $(nballoc ÷ neval)")
        end
    end

    if doprint(control, 'n') && !state.system.is_linear
       infoout("  [n]ewton: $(round(t, sigdigits = 3)) seconds asm: $(round(100 * tasm / t, sigdigits = 3))%, linsolve: $(round(100 * tlinsolve / t, sigdigits = 3))%")
    end

    if doprint(control, 'l') && state.system.is_linear
       infoout("  [l]inear($(nameof(typeof(method_linear)))): $(round(t, sigdigits = 3)) seconds")
    end

    solution.history = nlhistory
    return solution
end

function evaluate_residual_and_jacobian!(state, u; params = Float64[], t = 0.0, tstep = Inf, embed = 0.0)
    _initialize_dirichlet!(u, state.system, state.data)
    eval_and_assemble(
        state.system,
        u,
        u,
        state.residual,
        state.matrix,
        state.dudp,
        t,
        tstep,
        embed,
        state.data,
        params
    )
    return flush!(state.matrix)
end


"""
    evaluate_residual_and_jacobian(system,u;
                                   t=0.0, tstep=Inf,embed=0.0)

Evaluate residual and jacobian at solution value u.
Returns a solution vector containing a copy of  residual, and an ExendableSparseMatrix
containing a copy of the linearization at u.

"""
function evaluate_residual_and_jacobian(sys, u; kwargs...)
    state = SystemState(sys)
    eval_and_linearize!(state, u; kwargs...)
    return copy(state.residual), copy(state.matrix)
end


"""
        solve_transient(inival, system, times; kwargs...)
Solve transient or embedding problem.
"""
function solve_transient!(
        state,
        inival,
        lambdas;
        control = SolverControl(),
        transient = true, # choose between transient and stationary (embedding) case
        time = 0.0,
        params = zeros(0),
        kwargs...,
    )

    # rounding in output
    rd(x) = round(x; sigdigits = 5)

    # Set initial value of Δλ etc
    if transient # λ is time
        λstr = "t"
        Δλ = control.Δt
        Δλ_min = control.Δt_min
        Δλ_max = control.Δt_max
        Δλ_grow = control.Δt_grow
        Δλ_decrease = control.Δt_decrease
    else  # λ is embedding parameter
        λstr = "p"
        Δλ = control.Δp
        Δλ_min = control.Δp_min
        Δλ_max = control.Δp_max
        Δλ_grow = control.Δp_grow
        Δλ_decrease = control.Δp_decrease
    end
    Δu_opt = control.Δu_opt
    Δu_max_factor = control.Δu_max_factor

    allhistory = TransientSolverHistory()

    solution = copy(inival)
    oldsolution = copy(inival) # we need a copy as it is later overwritten

    t0 = 0.0
    if transient
        # Initialize transient solution struct
        tsol = TransientSolution(Float64(lambdas[1]), copy(solution); in_memory = control.in_memory)
        # Initialize Dirichlet boundary values
        _initialize_dirichlet!(solution, state.system, state.data; time, λ = Float64(lambdas[1]), params)
    else
        # If not transient, solve for first embedding lambdas[1]
        _initialize_dirichlet!(solution, state.system, state.data; time, λ = Float64(lambdas[1]), params)
        control.pre(solution, Float64(lambdas[1]))
        t0 = @elapsed solution = solve_step!(
            state,
            solution,
            oldsolution,
            control,
            time,
            Inf,
            Float64(lambdas[1]),
            params
        )

        control.post(solution, oldsolution, lambdas[1], 0)

        if control.log
            push!(allhistory, solution.history)
            push!(allhistory.times, lambdas[1])
            Δu = control.delta(state.system, solution, oldsolution, lambdas[1], 0)
            push!(allhistory.updates, Δu)
        end
        oldsolution .= solution

        tsol = TransientSolution(Float64(lambdas[1]), copy(solution); in_memory = control.in_memory)
    end


    if doprint(control, 'e')
        infoout("[e]volution: start in $(extrema(lambdas))")
    end

    λ0 = 0
    istep = 0
    solved = false
    # Outer loop over embedding params/ time values
    t1 = @elapsed for i in 1:(length(lambdas) - 1)
        Δλ = max(Δλ, Δλ_min)
        λstart = lambdas[i]
        λend = lambdas[i + 1]
        λ = Float64(λstart)

        # Inner loop between two embedding params/time values
        while λ < λend
            solved = false
            λ0 = λ
            Δu = 0.0

            # Try to solve, possibly with stepsize decrease
            while !solved
                solved = true
                forced = false
                errored = false
                try # check for non-converging newton
                    λ = λ0 + Δλ
                    control.pre(solution, λ)

                    if transient # λ is time
                        _time = λ
                        _tstep = Δλ
                        _embedparam = 0.0
                    else # λ is embedding parameter
                        _time = time
                        _tstep = Inf
                        _embedparam = λ
                    end
                    solution = solve_step!(
                        state,
                        solution,
                        oldsolution,
                        control,
                        _time,
                        _tstep,
                        _embedparam,
                        params
                    )
                catch err
                    err = "Problem at $(λstr)=$(λ |> rd), Δ$(λstr)=$(Δλ |> rd):\n$(err)"
                    if (control.handle_exceptions)
                        _print_error(err, stacktrace(catch_backtrace()))
                    else
                        rethrow(err)
                    end
                    solved = false
                    errored = true
                end
                if solved
                    Δu = control.delta(state.system, solution, oldsolution, λ, Δλ)
                    if Δu > Δu_max_factor * Δu_opt
                        solved = false
                    end
                end
                if !solved
                    if Δλ ≈ Δλ_min
                        if !(control.force_first_step && istep == 0)
                            err = """
                            At $(λstr)=$(λ |> rd): Δ$(λstr)_min=$(Δλ_min |> rd) reached while Δu/Δu_opt=$(Δu / Δu_opt |> rd).
                            Returning prematurely before $(λstr)[end]=$(lambdas[end] |> rd) 
                            """
                            if control.handle_exceptions
                                warnout(err)
                            else
                                throw(ErrorException(err))
                            end
                            break # give up lowering stepsize, break out if "while !solved" loop
                        elseif !errored
                            if doprint(control, 'e')
                                infoout("[e]volution:  forced first timestep: Δu/Δu_opt=$(Δu / Δu_opt |> rd)")
                            end
                            forced = true
                            solved = true
                        else
                            err = "Convergence problem in first timestep"
                            if control.handle_exceptions
                                warnout(err)
                            else
                                throw(ErrorException(err))
                            end

                            break
                        end
                    else
                        # reduce time step
                        Δλ = max(Δλ_min, Δλ * Δλ_decrease)
                        if doprint(control, 'e')
                            out=@sprintf("[e]volution:  Δu/Δu_opt=%.3e => retry: Δ%s=%.3e\n", Δu / Δu_opt, λstr, Δλ)
                            infoout(out)
                        end
                    end
                end
            end # while !solved

            # Advance step
            if solved
                istep = istep + 1
                if doprint(control, 'e')
                    out= @sprintf(
                        "[e]volution: step=%d %s=%.3e Δ%s=%.3e Δu=%.3e\n",
                        istep,
                        λstr,
                        λ,
                        λstr,
                        Δλ,
                        Δu
                    )
                    infoout(out)
                end
                if control.log
                    push!(allhistory, solution.history)
                    push!(allhistory.updates, Δu)
                    push!(allhistory.times, λ)
                end
                if control.store_all
                    append!(tsol, λ, copy(solution))
                end
                control.post(solution, oldsolution, λ, Δλ)
                oldsolution .= solution

                # Adjust last timesteps, so there will be no accidental
                # very small timestep
                steps_to_go = ceil((λend - λ) / Δλ)
                λ_predict = λend - λ
                if steps_to_go < control.num_final_steps
                    λ_predict = (λend - λ) / steps_to_go
                end

                # see fixed_timesteps!()
                if Δλ_max ≈ Δλ_min
                    λ_predict = Δλ_max
                end

                if λ < λend
                    #  ### account for close last timestep
                    Δλ = min(
                        Δλ_max,
                        Δλ * Δλ_grow,
                        Δλ * Δu_opt / (Δu + 1.0e-14),
                        λ_predict,
                        λend - λ
                    )

                    # adjust for roundoff error, so there will be no accidental
                    # very tiny timestep
                    if isapprox(λ + Δλ, λend; atol = 1.0e-15, rtol = 1.0e-15)
                        Δλ = λend - λ
                    end
                end
            else
                break # break out of inner loop over timestep
            end # if solved
        end # while λ<λ_end

        if !control.store_all # store last solution obtained
            append!(tsol, λ0, copy(solution))
        end

        control.sample(solution, λ0)

        if solved
            if !(λ ≈ lambdas[i + 1]) # check end of interval has been reached in inner loop
                warnout("λ=$(λ), lambdas[i+1]=$(lambdas[i + 1])")
            end
        else
            break # emergency exit
        end
    end # for i = 1:(length(lambdas)-1), end outer loop

    if doprint(control, 'e')
        infoout("[e]volution:  $(round(t0 + t1, sigdigits = 3)) seconds")
    end

    tsol.history = allhistory
    return tsol
end


#####################################################################
"""
    solve!(state; kwargs...)

Built-in solution method for [`VoronoiFVM.System`](@ref).  

Solves finite volume system the satate is belonging to.
Mutates the state and returns the solution.

For the keyword argumentsm see [`VoronoiFVM.solve`](@ref).
"""
function CommonSolve.solve!(
        state::VoronoiFVM.SystemState;
        inival = 0,
        data = nothing,
        params = zeros(0),
        control = VoronoiFVM.SolverControl(),
        time = 0.0,
        tstep = Inf,
        kwargs...,
    )
    fix_deprecations!(control)
    if isa(inival, Number) || isa(inival, Matrix)
        inival = unknowns(state.system; inival = inival)
    elseif !isdensesystem(state.system) && isa(inival, SparseMatrixCSC)
        inival = SparseSolutionArray(inival)
    elseif !VoronoiFVM.isunknownsof(inival, state.system)
        @error "wrong type of inival: $(typeof(inival))"
    end

    for pair in kwargs
        if first(pair) != :times &&
                first(pair) != :embed &&
                hasfield(SolverControl, first(pair))
            setproperty!(control, first(pair), last(pair))
        end
    end

    if !isnothing(data)
        state.data = data
    end
    state.params .= params


    return if haskey(kwargs, :times) && !isnothing(kwargs[:times])
        solve_transient!(
            state,
            inival,
            kwargs[:times];
            control,
            transient = true,
            params,
            time = kwargs[:times][1]
        )
    elseif haskey(kwargs, :embed) && !isnothing(kwargs[:embed])
        solve_transient!(
            state,
            inival,
            kwargs[:embed];
            transient = false,
            control,
            params,
            time,
        )
    else
        solve_step!(
            state,
            unknowns(state.system),
            inival,
            control,
            time,
            tstep,
            0.0,
            params
        )
    end
end

"""
    solve(system; kwargs...)
    
Built-in solution method for [`VoronoiFVM.System`](@ref).  
    
Keyword arguments:
- General for all solvers 
   - `inival` (default: 0) : Array created via [`unknowns`](@ref) or  number giving the initial value.
   - `control` (default: nothing): Pass instance of [`SolverControl`](@ref)
   -  All elements of [`SolverControl`](@ref) can be used as kwargs. Eventually overwrites values given via `control`
   - `params`: Parameters (Parameter handling is experimental and may change)
    
- __Stationary solver__:
  Invoked if neither `times` nor `embed`, nor `tstep` are given as keyword argument.
  - `time` (default: `0.0`): Set time value.
  Returns a [`DenseSolutionArray`](@ref) or [`SparseSolutionArray`](@ref)

- __Embedding (homotopy) solver__: Invoked if `embed` kwarg is given.
  Use homotopy embedding + damped Newton's method  to 
  solve stationary problem or to solve series of parameter dependent problems.
  Parameter step control is performed according to solver control data.  kwargs and default values are:
  - `embed` (default: `nothing` ): vector of parameter values to be reached exactly
  In addition,  all kwargs of the implicit Euler solver (besides `times`) are handled.  
  Returns a transient solution object `sol` containing the stored solution(s),  see [`TransientSolution`](@ref).
  
- __Implicit Euler transient solver__: Invoked if `times` kwarg is given.
  Use implicit Euler method  + damped   Newton's method  to 
  solve time dependent problem. Time step control is performed
  according to solver control data.  kwargs and default values are:
  - `times` (default: `nothing` ): vector of time values to be reached exactly
  - `pre` (default: `(sol,t)->nothing` ):  callback invoked before each time step
  - `post`  (default:  `(sol,oldsol, t, Δt)->nothing` ): callback invoked after each time step
  - `sample` (default:  `(sol,t)->nothing` ): callback invoked after timestep for all times in `times[2:end]`.
  - `delta` (default:  `(system, u,v,t, Δt)->norm(sys,u-v,Inf)` ):  Value  used to control the time step size `Δu`
  If `control.handle_error` is true, if time step solution  throws an error,
  stepsize  is lowered, and  step solution is called again with a smaller time value.
  If `control.Δt<control.Δt_min`, solution is aborted with error.
  Returns a transient solution object `sol` containing the stored solution,  see [`TransientSolution`](@ref).
  
- __Implicit Euler timestep solver__.  Invoked if `tstep` kwarg is given. Solve one time step of the implicit Euler method.
  - `time` (default: `0`): Set time value. 
  - `tstep`: time step
  Returns a [`DenseSolutionArray`](@ref) or [`SparseSolutionArray`](@ref)
"""
function CommonSolve.solve(sys::VoronoiFVM.AbstractSystem; data = sys.physics.data, kwargs...)
    state = SystemState(sys; data)
    return solve!(state; kwargs...)
end
