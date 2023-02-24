Base.@kwdef struct RVEA <: AbstractAlgorithm
    N::Int = 100
    alpha = 2
    fr = 0.1
end

name(_::RVEA) = "RVEA"

mutable struct RVEAState <: AbstractOptimizerState
    iteration::Any
    fcalls::Any
    start_time::Any
    stop_time::Any

    N::Any                           # population size
    population::Population           # population
    V0::Any
    V::Any                           # reference vectors
end
pfront(s::RVEAState) = get_non_dominated_solutions(s.population)
value(s::RVEAState) = objectives(pfront(s))
minimizer(s::RVEAState) = objectives(pfront(s))
copy(s::RVEAState) = RVEAState(
    s.iteration,
    s.fcalls,
    s.start_time,
    s.stop_time,
    s.N,
    copy(s.population),
    copy(s.V0),
    copy(s.V)
)

function initial_state(algorithm::RVEA, problem, options)
    population = initial_population(algorithm, problem)

    # Generate reference vectors
    V, N = TwoLayer(algorithm.N, problem.M)

    fcalls = evaluate!(problem, population)

    return RVEAState(
        0,
        fcalls,
        0,
        0,
        N,
        population,
        V,
        copy(V)
    )
end

function update_state!(algorithm::RVEA, state, problem, options)
    parents = state.population
    selected = shuffle(1:state.N)
    offspring = evolute(state, problem, selected, rng=options.rng)
    combine = [parents; offspring]
    snapshot = copy(combine)

    # Objective value translation
    z = minimum(objectives(combine), dims=1)
    for ind in combine
        ind.objectives = (ind.objectives - z')[:,1]
    end
    
    # Calculate the smallest angle value between each vector and others
    cosine = 1 .- pairwise(CosineDist(), state.V, state.V, dims=1)
    cosine[diagind(cosine)] .= 0
    gamma = minimum(acos.(cosine), dims=2)

    # Associate each solution to a reference vector
    angle = acos.(1 .- pairwise(CosineDist(), objectives(combine), state.V, dims=1))
    _, associate = findmin(angle, dims=2)
    indexes = [index[2] for index in associate]
    
    # Select one solution for each reference vector
    next = fill(0, size(state.V)[1])
    theta = (state.fcalls / problem.maxFE) ^ algorithm.alpha
    for i in unique(indexes)
        current = map(assoc->assoc[1], filter(assoc->i==assoc[2], associate))
        if !isempty(current)
            APD = (1 .+ problem.M*theta*angle[current,i]/gamma[i]) .* sqrt.(sum(objectives(combine[current]).^2, dims=2))
            _, best = findmin(APD)
            next[i] = current[best]
        end
    end

    # Restore objective values
    combine = snapshot

    state.population = combine[filter(i->i!=0,next)]
    state.N = length(state.population)

    if mod(ceil(state.fcalls/state.N), ceil(algorithm.fr*problem.maxFE/state.N)) == 0
        state.V = state.V0 .* repeat(maximum(objectives(state.population), dims=1) - minimum(objectives(state.population), dims=1), size(state.V0)[1], 1)
    end

    return false
end
