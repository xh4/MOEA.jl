function optimize(
    problem,
    algorithm;
    options::Options = Options()
)::OptimizationResult
    start_time = stop_time = time()

    state = initial_state(algorithm, problem, options)
    state.start_time = start_time
    state.stop_time = time()

    iteration = 1
    stopped = false

    trace = []
    # Ignore first state
    # if options.store_trace
    #     push!(trace, copy(state))
    # end

    while !stopped
        iteration += 1
        state.iteration = iteration

        state.start_time = time()

        should_break = update_state!(
            algorithm,
            state,
            problem,
            options
        )

        state.stop_time = time()

        if (options.store_trace)
            push!(trace, copy(state))
        end
        if (options.state_callback !== nothing)
            options.state_callback(state)
        end

        if state.fcalls >= problem.maxFE
            stopped = true
        end

        should_break && break
    end

    stop_time = time()

    if options.finish_callback !== nothing
        options.finish_callback(state)
    end

    return OptimizationResult(
        problem,
        algorithm,
        state,
        iterations,
        trace,
        stop_time - start_time
    )
end
