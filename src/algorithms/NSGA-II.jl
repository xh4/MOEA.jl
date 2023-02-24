Base.@kwdef struct NSGAII <: AbstractAlgorithm
    N::Int = 100
end

name(_::NSGAII) = "NSGA-II"

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

    population
end
pfront(s::NSGAIIState) = get_non_dominated_solutions(s.population)
value(s::NSGAIIState) = objectives(s.population)
minimizer(s::NSGAIIState) = variables(s.population)
copy(s::NSGAIIState) = NSGAIIState(s.iteration, s.fcalls, s.start_time, s.stop_time,
                                 copy(s.population))

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
    population = map(NSGAIIIndividual, initial_population(algorithm, problem))
    fcalls = evaluate!(problem, population)
    nondominatedsort!(population)
    crowding_distance!(population)
    return NSGAIIState(0, fcalls, 0, 0, population)
end

function update_state!(
    algorithm::NSGAII,
    state,
    problem,
    options
)
    parents = state.population
    selected = tournament(2, select = twowaycomp)(hcat(ranks(parents), -distances(parents))', algorithm.N; rng = options.rng)
    offspring = evolute(state, problem, selected, rng=options.rng)
    combine = [parents; offspring]

    # Calculate ranks & crowding distances for individuals
    nondominatedsort!(combine)
    crowding_distance!(combine)

    # Select best individuals
    fitidx = Int[]
    F = calcF(combine)
    for f in F
        if length(fitidx) + length(f) > algorithm.N
            idxs = f[reverse(sortperm(map(distance, combine[f])))]
            append!(fitidx, idxs[1:(algorithm.N-length(fitidx))])
            break
        else
            append!(fitidx, f)
        end
    end

    state.population .= combine[fitidx]
    return false
end
