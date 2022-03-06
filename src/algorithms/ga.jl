#=
Implementation of Genetic Algorithm

The constructor takes following keyword arguments:

- `populationSize`: The size of the population
- `crossoverRate`: The fraction of the population at the next generation, not including elite children, that is created by the crossover function.
- `mutationRate`: Probability of chromosome to be mutated
- `ɛ`/`epsilon`: Positive integer specifies how many individuals in the current generation are guaranteed to survive to the next generation. Floating number specifies fraction of population.
- `selection`: [Selection](@ref) function
- `crossover`: [Crossover](@ref) function (default: `identity`)
- `mutation`: [Mutation](@ref) function (default: `identity`)
- `metrics` is a collection of convergence metrics.
=#
struct GA <: Optimizer
    populationSize::Int
    crossoverRate::Float64
    mutationRate::Float64
    ɛ::Real
    selection
    crossover
    mutation
    metrics::ConvergenceMetrics

    GA(;
        populationSize::Int = 50,
        crossoverRate::Float64 = 0.8,
        mutationRate::Float64 = 0.1,
        ɛ::Real = 0,
        epsilon::Real = ɛ,
        selection = ((x, n) -> 1:n),
        crossover = identity,
        mutation = identity,
        metrics = ConvergenceMetric[AbsDiff(1e-12)],
    ) = new(
        populationSize,
        crossoverRate,
        mutationRate,
        epsilon,
        selection,
        crossover,
        mutation,
        metrics,
    )
end

population_size(method::GA) = method.populationSize
default_options(method::GA) = (iterations = 1000,)
summary(m::GA) =
    "GA[P=$(m.populationSize),x=$(m.crossoverRate),μ=$(m.mutationRate),ɛ=$(m.ɛ)]"
show(io::IO, m::GA) = print(io, summary(m))

mutable struct GAState <: AbstractOptimizerState
    iteration
    start_time
    stop_time

    metrics::ConvergenceMetrics # collection of convergence metrics
    eliteSize::Int
    population
    fittest
end

value(s::GAState) = objective(s.fittest)
minimizer(s::GAState) = variable(s.fittest)
copy(s::GAState) = GAState(s.iteration, s.start_time, s.stop_time, copy(s.metrics), s.eliteSize, copy(s.population), copy(s.fittest))

"""Initialization of GA algorithm state"""
function initial_state(method::GA, objective, population, options)
    eliteSize = isa(method.ɛ, Int) ? method.ɛ : round(Int, method.ɛ * method.populationSize)

    fitness = value!(objective, population)
    _, idx = findmin(fitness)

    return GAState(0, 0, 0, copy(method.metrics), eliteSize, population, copy(population[idx]))
end

function update_state!(
    method::GA,
    state,
    obj,
    constraints,
    options
)
    parents = state.population

    populationSize = method.populationSize
    rng = options.rng
    offspring = copy(parents)

    # select offspring
    selected = method.selection(state.population, populationSize, rng = rng)

    # perform mating
    offspringSize = populationSize - state.eliteSize
    recombine!(offspring, parents, selected, method, offspringSize, rng = rng)

    # Elitism (copy population individuals before they pass to the offspring & get mutated)
    idxs = sortperm(vec(objective(state.population)))
    for i = 1:state.eliteSize
        subs = offspringSize + i
        offspring[subs] = copy(parents[idxs[i]])
    end

    # perform mutation
    mutate!(offspring, method, constraints, rng = rng)

    # calculate fitness of the population
    evaluate!(obj, offspring, constraints)

    # select the best individual
    _, idx = findmin(objective(offspring))
    state.fittest = offspring[idx]

    # replace population
    parents .= offspring

    return false
end

function recombine!(
    offspring,
    parents,
    selected,
    method,
    n = length(selected);
    rng::AbstractRNG = Random.default_rng(),
)
    mates = ((i, i == n ? i - 1 : i + 1) for i = 1:2:n)
    for (i, j) in mates
        p1, p2 = parents[selected[i]], parents[selected[j]]
        if rand(rng) < method.crossoverRate
            offspring[i], offspring[j] = method.crossover(p1, p2, rng = rng)
        else
            offspring[i], offspring[j] = p1, p2
        end
    end

end

function mutate!(population, method, constraints; rng::AbstractRNG = Random.default_rng())
    n = length(population)
    for i = 1:n
        if rand(rng) < method.mutationRate
            method.mutation(population[i], rng = rng)
        end

        (constraints, population[i])
    end
end

function evaluate!(objective, population, constraints)
    # calculate fitness of the population
    value!(objective, population)
    # apply penalty to fitness
    # penalty!(fitness, constraints, population)
end
