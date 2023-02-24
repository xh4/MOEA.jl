Base.@kwdef struct MOEADDE <: AbstractAlgorithm
    N::Int = 100            # number of the subproblems
    T::Int = ceil(N/10)     # number of the weight vectors in the neighborhoor of each weight vector
    F::Real = 0.5
    delta::Real = 0.9
    nr::Int = 2
end

name(_::MOEADDE) = "MOEA/D-DE"

mutable struct MOEADDEState <: AbstractOptimizerState
    iteration
    fcalls
    start_time
    stop_time

    N                           # population size
    population                  # population
    W                           # weight vectors
    B                           # neighbours of each solution
    z                           # the best value found so far
    EP                          # nondominated solutions
end
pfront(s::MOEADDEState) = get_non_dominated_solutions(s.population)
value(s::MOEADDEState) = objectives(pfront(s))
minimizer(s::MOEADDEState) = objectives(pfront(s))
copy(s::MOEADDEState) = MOEADDEState(s.iteration, s.fcalls, s.start_time, s.stop_time,
                                 s.N,
                                 copy(s.population),
                                 copy(s.W),
                                 copy(s.B),
                                 copy(s.z),
                                 copy(s.EP))

function initial_state(algorithm::MOEADDE, problem, options)
    population = [Individual(rand(problem.D)) for i in 1:algorithm.N]

    fcalls = evaluate!(problem, population)

    # Generate weight vectors
    W, N = TwoLayer(population_size(algorithm), problem.M)

    population = population[1:N]

    # Set the neighbours of each weight vector
    B = pairwise(Euclidean(), W, W, dims=1)
    B = sortperm(B)
    B = B[:,2:algorithm.T]

    # Set reference point
    z = Matrix(minimum(objectives(population), dims=1)')

    EP = Individual[]
    return MOEADDEState(0, fcalls, 0, 0, N, population, W, B, z, EP)
end

function update_state!(
    algorithm::MOEADDE,
    state,
    problem,
    options
)
    for i = 1:state.N
        # Selection of mating/update range
        if rand() < algorithm.delta
            P = state.B[i,:]
        else
            P = 1:state.N
        end

        # Reproduction
        xi = rand(options.rng, P, 2)
        y = evolute1_DE(state, problem, xi, rng=options.rng)

        # Update reference point
        state.z = minimum(hcat(state.z, objectives(y)), dims=2)

        # Update of neighboring solutions
        c = 0
        while c < algorithm.nr && !isempty(P)
            j = rand(options.rng, P)
            x = state.population[j]
            g_x, _ = findmax(state.W[j,:] .* abs.(objectives(x) - state.z))
            g_y, _ = findmax(state.W[j,:] .* abs.(objectives(y) - state.z))
            if g_y <= g_x
                state.population[j] = y
                c += 1
            end
            P = filter(i -> i!=j, P)
        end

        # Update of EP
        # dominates = map((p) -> dominate(y, p), state.EP)
        # deleteat!(state.EP, findall(isequal(1), dominates))
        # if isempty(findall(isequal(-1), dominates))
        #     pushfirst!(state.EP, y)
        # end
    end

    return false
end
