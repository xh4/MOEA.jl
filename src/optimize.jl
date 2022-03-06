"""
    optimize(f[, F], individual, algorithm[, opts])
    optimize(f[, F], constraints, algorithm[, opts])
    optimize(f[, F], constraints, individual, algorithm[, opts])
    optimize(f[, F], constraints, algorithm, population[, opts])

- 对于多目标优化，目标值 `F` 必须提供.
"""
optimize(
    f::TC,
    individual,
    method::M,
    opts::Options = Options(; default_options(method)...),
) where {TC,M<:Optimizer} = optimize(f::TC, NoConstraints(), individual, method, opts)
optimize(
    f,
    F::AbstractVector,
    individual,
    method::M,
    opts::Options = Options(; default_options(method)...),
) where {TC,M<:Optimizer} = optimize(f, F, NoConstraints(), individual, method, opts)
optimize(
    f::TC,
    bounds::ConstraintBounds,
    method::M,
    opts::Options = Options(; default_options(method)...),
) where {TC,M<:Optimizer} = optimize(f, BoxConstraints(bounds), method, opts)
optimize(
    f::TC,
    F::TF,
    bounds::ConstraintBounds,
    method::M,
    opts::Options = Options(; default_options(method)...),
) where {TC,TF,M<:Optimizer} = optimize(f, F, BoxConstraints(bounds), method, opts)
function optimize(
    f::TC,
    constraints::C,
    method::M,
    opts::Options = Options(; default_options(method)...),
) where {TC,M<:Optimizer,C<:AbstractConstraints}
    population = Population(method, bounds(constraints), rng = opts.rng)
    optimize(f, constraints, method, population, opts)
end
function optimize(
    f::TC,
    F::TF,
    constraints::C,
    method::M,
    opts::Options = Options(; default_options(method)...),
) where {TC,TF,M<:Optimizer,C<:AbstractConstraints}
    population = Population(method, bounds(constraints), rng = opts.rng)
    optimize(f, F, constraints, method, population, opts)
end
function optimize(
    f::TC,
    constraints::C,
    individual,
    method::M,
    opts::Options = Options(; default_options(method)...),
) where {TC,M<:Optimizer,C<:AbstractConstraints}
    population = Population(method, individual, rng = opts.rng)
    optimize(f, constraints, method, population, opts)
end
function optimize(
    f::TC,
    F::TF,
    constraints::C,
    individual,
    method::M,
    opts::Options = Options(; default_options(method)...),
) where {TC,TF,M<:Optimizer,C<:AbstractConstraints}
    population = Population(method, individual, rng = opts.rng)
    optimize(f, F, constraints, method, population, opts)
end
function optimize(
    f::TC,
    constraints::C,
    method::M,
    population,
    opts::Options = Options(; default_options(method)...),
) where {TC,M<:Optimizer,C<:AbstractConstraints}
    @assert length(population) > 0 "Population is empty"
    obj = Objective(f, first(population))
    optimize(obj, constraints, method, population, opts)
end
function optimize(
    f::TC,
    F::TF,
    constraints::C,
    method::M,
    population,
    opts::Options = Options(; default_options(method)...),
) where {TC,TF,M<:Optimizer,C<:AbstractConstraints}
    @assert length(population) > 0 "Population is empty"
    obj = Objective(f, first(population), F)
    optimize(obj, constraints, method, population, opts)
end

function optimize(
    objective::D,
    constraints::C,
    method::M,
    population::AbstractArray,
    options::Options = Options(; default_options(method)...),
    state = initial_state(method, objective, population, options),
)::OptimizationResult where {D<:AbstractObjective,C<:AbstractConstraints,M<:Optimizer}

    iteration = 0
    stopped = false
    converged, counter_tol = false, 0 # tolerance convergence
    is_moo = ismultiobjective(objective)
    start_time = stop_time = time()

    trace = []
    if options.store_trace
        push!(trace, copy(state))
    end

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
