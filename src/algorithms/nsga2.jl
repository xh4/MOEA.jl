#=
Non-dominated Sorting Genetic Algorithm (NSGA-II) for Multi-objective Optimization

The constructor takes following keyword arguments:

- `populationSize`: The size of the population
- `crossoverRate`: The fraction of the population at the next generation, that is created by the crossover function
- `mutationRate`: Probability of chromosome to be mutated
- `selection`: [Selection](@ref) function (default: `tournament`)
- `crossover`: [Crossover](@ref) function (default: `SBX`)
- `mutation`: [Mutation](@ref) function (default: `PLM`)
- `metrics` is a collection of convergence metrics.
=#
struct NSGA2 <: Optimizer
    populationSize::Int
    crossoverRate::Float64
    mutationRate::Float64
    selection
    crossover
    mutation
    metrics::ConvergenceMetrics

    NSGA2(;
        populationSize::Int = 50,
        crossoverRate::Float64 = 0.9,
        mutationRate::Float64 = 0.2,
        selection = tournament(2, select = twowaycomp),
        crossover = SBX(),
        mutation = PLM(),
        metrics = ConvergenceMetric[GD(), GD(true)],
    ) = new(
        populationSize,
        crossoverRate,
        mutationRate,
        selection,
        crossover,
        mutation,
        metrics,
    )
end

population_size(method::NSGA2) = method.populationSize
default_options(method::NSGA2) = (iterations = 1000,)
summary(m::NSGA2) =
    "NSGA-II[P=$(m.populationSize),x=$(m.crossoverRate),μ=$(m.mutationRate)]"
show(io::IO, m::NSGA2) = print(io, summary(m))

mutable struct NSGA2State <: AbstractOptimizerState
    iteration
    start_time
    stop_time

    metrics::ConvergenceMetrics # collection of convergence metrics
    population                  # population
    pfront                      # individuals of the first Pareto front
    ranks::Vector{Int}          # individual ranks
    crowding::Vector{Float64}   # individual crowding distance
end
value(s::NSGA2State) = objective(s.pfront)
minimizer(s::NSGA2State) = variable(s.pfront)
copy(s::NSGA2State) = NSGA2State(s.iteration, s.start_time, s.stop_time,
                                 copy(s.metrics),
                                 copy(s.population),
                                 copy(s.pfront), copy(s.ranks), copy(s.crowding))

"""Initialization of NSGA2 algorithm state"""
function initial_state(method::NSGA2, objective, population, options)
    value!(objective, population)

    # setup initial state
    ranks = vcat(fill(1, method.populationSize), fill(2, method.populationSize))
    crowding =
        vcat(fill(zero(Float64), method.populationSize), fill(typemax(Float64), method.populationSize))
    return NSGA2State(0, 0, 0, copy(method.metrics), population, population, ranks, crowding)
end

function update_state!(
    method::NSGA2,
    state,
    objfun,
    constraints,
    options
)
    parents = state.population

    populationSize = method.populationSize
    rng = options.rng

    # select offspring
    specFit = StackView(state.ranks, state.crowding, dims = 1)
    selected = method.selection(view(specFit, :, 1:populationSize), populationSize; rng = rng)

    offspring = copy(parents)

    # perform mating
    recombine!(offspring, parents, selected, method)

    # perform mutation
    mutate!(offspring, method, constraints, rng = rng)

    # handle constraint
    map(x -> apply!(constraints, variable(x)), offspring)

    evaluate!(objfun, offspring, constraints)

    population = vcat(parents, offspring)
 
    # calculate ranks & crowding for population
    F = nondominatedsort!(state.ranks, objective(population)')
    crowding_distance!(state.crowding, objective(population)', F)

    # select best individuals
    fitidx = Int[] 
    for f in F
        if length(fitidx) + length(f) > populationSize
            idxs = sortperm(view(state.crowding, f))
            append!(fitidx, idxs[1:(populationSize-length(fitidx))])
            break
        else
            append!(fitidx, f)
        end
    end
    # designate the first Pareto front individuals as the fittest
    fidx = length(F[1]) > populationSize ? fitidx : F[1]
    state.pfront = copy(population[fidx])

    # construct new parent population
    parents .= population[fitidx]

    return false
end