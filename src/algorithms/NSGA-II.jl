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
Base.@kwdef struct NSGAII <: AbstractOptimizer
    populationSize::Int = 50
    crossoverRate::Float64 = 0.9
    mutationRate::Float64 = 0.2
    selection = tournament(2, select = twowaycomp)
    crossover = SBX()
    mutation = PLM()
    metrics::ConvergenceMetrics = ConvergenceMetric[GD(), GD(true)]
end

population_size(method::NSGAII) = method.populationSize
default_options(method::NSGAII) = (iterations = 1000,)
summary(m::NSGAII) =
    "NSGA-II[P=$(m.populationSize),x=$(m.crossoverRate),μ=$(m.mutationRate)]"
show(io::IO, m::NSGAII) = print(io, summary(m))

mutable struct NSGAIIIndividual <: AbstractIndividual
    variables::AbstractVector
    objectives::Union{AbstractVector, Number, Nothing}
    constraints::Union{AbstractVector, Nothing}
    rank
    distance

    NSGAIIIndividual(var) = new(var, nothing, nothing)
    NSGAIIIndividual(var, obj) = new(var, obj, nothing)
    NSGAIIIndividual(var, obj, cst) = new(var, obj, cst)
    NSGAIIIndividual(var, obj, cst, rank, distance) = new(var, obj, cst, rank, distance)
end
variables(i::NSGAIIIndividual) = i.variables
objectives(i::NSGAIIIndividual) = i.objectives
constraints(i::NSGAIIIndividual) = i.constraints
rank(i::NSGAIIIndividual) = i.rank
distance(i::NSGAIIIndividual) = i.distance
copy(i::NSGAIIIndividual) = NSGAIIIndividual(copy(i.variables), copy(i.objectives), copy(i.constraints), copy(i.rank), copy(i.distance))
NSGAIIIndividual(i::Individual) = NSGAIIIndividual(copy(i.variables), copy(i.objectives), copy(i.constraints), nothing, nothing)

mutable struct NSGAIIState <: AbstractOptimizerState
    iteration
    start_time
    stop_time
    metrics::ConvergenceMetrics # collection of convergence metrics

    population                  # population
    pfront                      # individuals of the first Pareto front
end
pfront(s::NSGAIIState) = s.pfront
value(s::NSGAIIState) = objectives(pfront(s))
minimizer(s::NSGAIIState) = variables(pfront(s))
copy(s::NSGAIIState) = NSGAIIState(s.iteration, s.start_time, s.stop_time, copy(s.metrics),
                                 copy(s.population),
                                 copy(s.pfront))

function ranks(P::Population)
    map(rank, P)
end
function distances(P::Population)
    map(distance, P)
end
function maxF(P::Population)
    maxF, _ = findmax(map(rank, P))
    maxF
end
function calcF(P::Population)
    F = [Vector{Int64}() for _ in 1:maxF(P)]
    for i = 1:length(P)
        push!(F[rank(P[i])], i)
    end
    F
end

function initial_state(method::NSGAII, objective, population, options)
    population = map(NSGAIIIndividual, population)
    value!(objective, population)
    nondominatedsort!(population)
    crowding_distance!(population)
    pfront = filter(i -> rank(i) == 1, population)
    return NSGAIIState(0, 0, 0, copy(method.metrics), population, pfront)
end

function update_state!(
    method::NSGAII,
    state,
    objective,
    constraints,
    options
)
    parents = state.population

    populationSize = method.populationSize
    rng = options.rng

    F = calcF(parents)
    for f = F
        d = sortperm(distances(parents[f]))
        f = shuffle(f)
    end
    selected = collect(Iterators.flatten(F))

    # select offspring
    # selected = method.selection(hcat(ranks(parents), distances(parents))', populationSize; rng = rng)

    offspring = copy(parents)

    # perform mating
    recombine!(offspring, parents, selected, method)

    # perform mutation
    mutate!(offspring, method, constraints, rng = rng)

    # handle constraint
    map(x -> apply!(constraints, variables(x)), offspring)

    evaluate!(objective, offspring, constraints)

    population = vcat(parents, offspring)

    # calculate ranks & crowding for population
    nondominatedsort!(population)
    crowding_distance!(population)

    ### println(filter(d -> d == -Inf, map(distance, filter(i -> rank(i) == 1, population))))

    # select best individuals
    fitidx = Int[]
    F = calcF(population)
    for f in F
        if length(fitidx) + length(f) > populationSize
            idxs = sortperm(map(rank, population[f]))
            append!(fitidx, idxs[1:(populationSize-length(fitidx))])
            break
        else
            append!(fitidx, f)
        end
    end

    state.pfront = filter(i -> rank(i) == 1, population)

    # construct new parent population
    parents .= population[fitidx]

    return false
end
