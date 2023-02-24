Base.@kwdef struct MOEAD <: AbstractAlgorithm
    N::Int = 100          # number of the subproblems
    T::Int = ceil(N/10)   # number of the weight vectors in the neighborhoor of each weight vector
end

name(_::MOEAD) = "MOEA/D"

mutable struct MOEADState <: AbstractOptimizerState
    iteration
    fcalls
    start_time
    stop_time

    N                           # population size
    population::Population      # population
    m                           # number of objectives
    W                           # weight vectors
    B                           # neighbours of each solution
    z                           # the best value found so far
    EP                          # nondominated solutions
end
pfront(s::MOEADState) = get_non_dominated_solutions(s.population)
value(s::MOEADState) = objectives(pfront(s))
minimizer(s::MOEADState) = objectives(pfront(s))
copy(s::MOEADState) = MOEADState(s.iteration, s.fcalls, s.start_time, s.stop_time,
                                 s.N,
                                 copy(s.population),
                                 s.m,
                                 copy(s.W),
                                 copy(s.B),
                                 copy(s.z),
                                 copy(s.EP))

function initial_state(algorithm::MOEAD, problem, options)
    population = initial_population(algorithm, problem)

    fcalls = evaluate!(problem, population)

    m = length(objectives(first(population)))

    # Generate weight vectors
    W, N = TwoLayer(population_size(algorithm), m)

    population = population[1:N]

    # Set the neighbours of each weight vector
    B = pairwise(Euclidean(), W, W, dims=1)
    B = sortperm(B)
    B = B[:,1:algorithm.T]

    # Set reference point
    z = Matrix(minimum(objectives(population), dims=1)')

    EP = Individual[]

    return MOEADState(0, fcalls, 0, 0, N, population, m, W, B, z, EP)
end

function update_state!(
    ::MOEAD,
    state,
    problem,
    options
)
    for i = 1:state.N
        x = state.population[i]

        # Reproduction
        xi = shuffle(options.rng, state.B[i,:])[1:2]
        y = evolute1(state, problem, xi, rng=options.rng)

        # Update reference point
        state.z = minimum(hcat(state.z, objectives(y)), dims=2)

        # Update of neighboring solutions
        for j = state.B[i,:]
            x = state.population[j]
            g_x = maximum(state.W[j,:] .* abs.(objectives(x) - state.z))
            g_y = maximum(state.W[j,:] .* abs.(objectives(y) - state.z))
            if g_y <= g_x
                state.population[j] = y
            end
        end

        # Update of EP
        state.EP = state.population
    end

    return false
end

function Base.sortperm(M::Matrix)
    reduce(vcat, map(sortperm, eachrow(M))')
end

function Base.sortperm!(M::Matrix)
    error("Not implemented")
end
