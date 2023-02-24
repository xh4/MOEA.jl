Base.@kwdef struct SMSEMOAdp <: AbstractAlgorithm
    N::Int = 100          # number of the subproblems
end

name(_::SMSEMOAdp) = "SMS-EMOA/dp"

Base.@kwdef mutable struct SMSEMOAdpIndividual <: AbstractIndividual
    variables
    objectives
    rank
    contribution

    SMSEMOAdpIndividual(var) = new(var, nothing, nothing, nothing)
    SMSEMOAdpIndividual(var, obj) = new(var, obj, nothing, nothing)
    SMSEMOAdpIndividual(var, obj, rank, contribution) = new(var, obj, rank, contribution)
end
variables(i::SMSEMOAdpIndividual) = i.variables
objectives(i::SMSEMOAdpIndividual) = i.objectives
rank(i::SMSEMOAdpIndividual) = i.rank
contribution(i::SMSEMOAdpIndividual) = i.contribution
copy(i::SMSEMOAdpIndividual) = SMSEMOAdpIndividual(copy(i.variables), copy(i.objectives), i.rank, i.contribution)
SMSEMOAdpIndividual(i::Individual) = SMSEMOAdpIndividual(copy(i.variables), copy(i.objectives), nothing, nothing)

Base.@kwdef mutable struct SMSEMOAdpState <: AbstractOptimizerState
    iteration
    fcalls
    start_time
    stop_time

    population                  # population
    pfront                      # individuals of the first Pareto front
end
pfront(s::SMSEMOAdpState) = s.pfront
value(s::SMSEMOAdpState) = objectives(pfront(s))
minimizer(s::SMSEMOAdpState) = objectives(pfront(s))
copy(s::SMSEMOAdpState) = SMSEMOAdpState(s.iteration, s.fcalls, s.start_time, s.stop_time, copy(s.population), copy(s.pfront))

function initial_state(algorithm::SMSEMOAdp, problem, options)
    population = map(SMSEMOAdpIndividual, initial_population(algorithm, problem))
    fcalls = evaluate!(problem, population)

    return SMSEMOAdpState(0, fcalls, 0, 0, population, pfront)
end

function update_state!(
    algorithm::SMSEMOAdp,
    state,
    problem,
    options
)
    population = state.population

    for i = 1:algorithm.N
        offspring = evolute1(state, problem, rand(options.rng, 1:algorithm.N, 2), rng=options.rng)
        push!(population, offspring)

        nondominatedsort!(population)
        F = calcF(population)
        R = F[end]
        if (length(R)) == 1
            r = R[1]
        else
            if length(F) > 1
                Q = reduce(vcat, F[1:end-1])
                r = R[argmax(map(s -> n_dominates(population[s], population[Q]), R))]
            else
                update_contribution!(population, population[R])
                r = R[argmin(contribution.(population[R]))]
            end
        end
        deleteat!(population, r)

        state.pfront = filter(i -> rank(i) == 1, population)
    end
    return false
end

function n_dominates(s, Q)
    length(filter(q -> dominate(q, s)==1, Q))
end