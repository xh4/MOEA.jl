Base.@kwdef struct NSGAII <: AbstractAlgorithm
    N::Int = 100
end

mutable struct NSGAIIIndividual <: AbstractIndividual
    variables::AbstractVector
    objectives::Union{AbstractVector, Number, Nothing}
    rank
    distance

    NSGAIIIndividual(var) = new(var, nothing)
    NSGAIIIndividual(var, obj) = new(var, obj)
    NSGAIIIndividual(var, obj, rank, distance) = new(var, obj, rank, distance)
end
variables(i::NSGAIIIndividual) = i.variables
objectives(i::NSGAIIIndividual) = i.objectives
rank(i::NSGAIIIndividual) = i.rank
distance(i::NSGAIIIndividual) = i.distance
copy(i::NSGAIIIndividual) = NSGAIIIndividual(copy(i.variables), copy(i.objectives), copy(i.rank), copy(i.distance))
NSGAIIIndividual(i::Individual) = NSGAIIIndividual(copy(i.variables), copy(i.objectives), nothing, nothing)

mutable struct NSGAIIState <: AbstractOptimizerState
    iteration
    fcalls
    start_time
    stop_time

    population                  # population
    pfront                      # individuals of the first Pareto front
end
pfront(s::NSGAIIState) = s.pfront
value(s::NSGAIIState) = objectives(s.population)
minimizer(s::NSGAIIState) = variables(s.population)
copy(s::NSGAIIState) = NSGAIIState(s.iteration, s.fcalls, s.start_time, s.stop_time,
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

function initial_state(algorithm::NSGAII, problem, options)
    population = [NSGAIIIndividual(rand(problem.D)) for i in 1:algorithm.N]
    fcalls = evaluate!(problem, population)
    nondominatedsort!(population)
    crowding_distance!(population)
    pfront = filter(i -> rank(i) == 1, population)
    return NSGAIIState(0, fcalls, 0, 0, population, pfront)
end

function update_state!(
    algorithm::NSGAII,
    state,
    problem,
    options
)
    parents = state.population
    offspring = copy(parents)

    # select parents for mating
    selected = tournament(2, select = twowaycomp)(hcat(ranks(parents), -distances(parents))', algorithm.N; rng = options.rng)

    # perform mating
    recombine!(offspring, parents, selected, SBX(), rng = options.rng)

    # perform mutation
    mutate!(offspring, PLM(lower=lower(problem), upper=upper(problem)), rng = options.rng)

    # handle constraint
    map(x -> apply!(problem.constraints, variables(x)), offspring)

    evaluate!(state, problem, offspring)

    population = [parents; offspring]

    # calculate ranks & crowding for population
    nondominatedsort!(population)
    crowding_distance!(population)

    ### println(filter(d -> d == -Inf, map(distance, filter(i -> rank(i) == 1, population))))

    # select best individuals
    fitidx = Int[]
    F = calcF(population)
    for f in F
        if length(fitidx) + length(f) > algorithm.N
            idxs = f[reverse(sortperm(map(distance, population[f])))]
            append!(fitidx, idxs[1:(algorithm.N-length(fitidx))])
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
