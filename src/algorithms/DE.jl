"""
The constructor takes following keyword arguments:
- `populationSize`: The size of the population
- `F`: the differentiation (mutation) scale factor (default: 0.9). It's usually defined in range ``F \\in (0, 1+]``
- `n`: the number of differences used in the perturbation (default: 1)
- `selection`: the selection strategy function (default: [`random`](@ref))
- `recombination`: the recombination functions (default: [`BINX(0.5)`](@ref))
- `K`: the recombination scale factor (default: 0.5*(F+1))
- `metrics` is a collection of convergence metrics.
"""
Base.@kwdef struct DE <: AbstractOptimizer
    populationSize::Integer = 50
    F::Real = 0.9
    n::Integer = 1
    K::Real = 0.5*(F+1)
    selection = random
    recombination = BINX(0.5)
    metrics::ConvergenceMetrics = [AbsDiff(1e-12)]
end

population_size(method::DE) = method.populationSize
default_options(method::DE) = (iterations=1000,)
summary(m::DE) = "DE/$(m.selection)/$(m.n)/$(m.recombination)"
show(io::IO,m::DE) = print(io, summary(m))

mutable struct DEState <: AbstractOptimizerState
    iteration
    start_time
    stop_time
    metrics::ConvergenceMetrics

    N::Int
    population
    fittest
end
value(s::DEState) = objectives(s.fittest)
minimizer(s::DEState) = variables(s.fittest)
copy(s::DEState) = DEState(s.iteration, s.start_time, s.stop_time, copy(s.metrics), s.N, copy(s.population), copy(s.fittest))

function initial_state(method::DE, objfun, population, options)
    T = typeof(value(objfun))
    individual = first(population)
    N = length(variables(individual)) # 决策变量维度
    value!(objfun, population)
    _, i = findmin(objectives(population))
    fittest = population[i]

    return DEState(0, 0, 0, copy(method.metrics), N, copy(population), copy(fittest))
end

function update_state!(
    method::DE,
    state,
    objective,
    constraints,
    options
)
    offspring = copy(state.population)

    targets = method.selection(state.population, method.populationSize)

    for (i, target) in enumerate(targets)

        # mutation
        idxs = randexcl(options.rng, 1:method.populationSize, [i], 3)
        donors = copy(state.population[idxs])
        donor = donors[1]
        differentiation!(donors[1], donors[2:end]; F=method.F)

        # crossover
        trial, _ = method.recombination(donor, target, rng=options.rng)
        offspring[i] = trial
    end

    value!(objective, offspring)

    for i in 1:method.populationSize
        # o = apply!(constraints, offspring[i])
        v = objectives(offspring[i])
        if (v <= objectives(state.population[i]))
            state.population[i] = offspring[i]
            if (v < objectives(state.fittest))
                state.fittest = copy(offspring[i])
            end
        end
    end

    return false
end
