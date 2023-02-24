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
struct GA <: AbstractOptimizer
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

value(s::GAState) = objectives(s.fittest)
minimizer(s::GAState) = variables(s.fittest)
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
    objective,
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
    idxs = sortperm(vec(objectives(state.population)))
    for i = 1:state.eliteSize
        subs = offspringSize + i
        offspring[subs] = copy(parents[idxs[i]])
    end

    # perform mutation
    mutate!(offspring, method, rng = rng)

    # calculate fitness of the population
    evaluate!(state, objective, offspring)

    # select the best individual
    _, idx = findmin(objectives(offspring))
    state.fittest = offspring[idx]

    # replace population
    parents .= offspring

    return false
end

function recombine!(
    offspring,
    parents,
    selected,
    operator,
    n = length(selected);
    crossover_rate = 1,
    rng::AbstractRNG = Random.default_rng(),
)
    mates = ((i, i == n ? i - 1 : i + 1) for i = 1:2:n)
    for (k, j) in mates
        pk, pj = selected[k], selected[j]
        p1, p2 = parents[pk], parents[pj]
        if rand(rng) < crossover_rate
            offspring[pk], offspring[pj] = operator(p1, p2, rng = rng)
        else
            offspring[pk], offspring[pj] = p1, p2
        end
    end
end

function mutate!(population, operator; mutation_rate = 1, rng::AbstractRNG = Random.default_rng())
    n = length(population)
    for i = 1:n
        if rand(rng) < mutation_rate
            population[i] = operator(population[i], rng = rng)
        end
    end
end

function evaluate!(state, problem, population::Population)
    state.fcalls += evaluate!(problem, population)
end
function evaluate!(problem, population::Population)
    fcalls = 0
    for i in population
        i.objectives = problem.fn(variables(i))
        fcalls += 1
    end
    fcalls
end
function evaluate!(state, problem, individual::AbstractIndividual)
    state.fcalls += evaluate!(problem, individual)
end
function evaluate!(problem, individual::AbstractIndividual)
    individual.objectives = problem.fn(variables(individual))
    1
end