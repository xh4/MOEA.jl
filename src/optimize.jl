function optimize(
    objective,
    method;
    constraints = nothing,
    population = nothing,
    options::Options = Options(; default_options(method)...),
    state = nothing,
)::OptimizationResult

    if constraints === nothing
        constraints = NoConstraints()
    end

    if isa(objective, Function)
        objective = Objective(objective, first(population))
    end

    if population === nothing
        population = initial_population(method, objective)
    end

    if state === nothing
        state = initial_state(method, objective, population, options)
    end

    iteration = 0
    stopped = false
    converged, counter_tol = false, 0 # tolerance convergence
    is_moo = ismultiobjective(objective)
    start_time = stop_time = time()

    trace = []
    # Ignore first state
    # if options.store_trace
    #     push!(trace, copy(state))
    # end

    # 终止条件
    # 1. 收敛
    # 2. 达到指定的 iteration
    # 3. 达到指定的运行时间
    # 4. terminate(state)?
    # 5. callback

    while !converged && !stopped && iteration < options.iterations
        iteration += 1
        state.iteration = iteration

        state.start_time = time()

        should_break = update_state!(
            method,
            state,
            objective,
            constraints,
            options
        )

        state.stop_time = time()

        converged = assess_convergence!(state)
        counter_tol = converged ? counter_tol + 1 : 0
        converged = converged && (counter_tol > options.successive_f_tol)
        converged = converged || terminate(state)

        if (options.store_trace)
            push!(trace, copy(state))
        end
        if (options.state_callback !== nothing)
            options.state_callback(state)
        end

        should_break && break
    end

    stop_time = time()

    return OptimizationResult(
        method,
        minimizer(state),
        is_moo ? NaN : value(state),
        iteration,
        iteration == options.iterations,
        converged,
        state.metrics,
        f_calls(objective),
        trace,
        options.time_limit,
        stop_time - start_time,
        is_moo,
    )
end
