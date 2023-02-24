Base.@kwdef struct MOEADUR <: AbstractAlgorithm
    N::Int = 100                 # number of the subproblems
    T::Int = ceil(N/10)          # number of the weight vectors in the neighborhoor of each weight vector
    delta = 0.8
    nr = ceil(N/100)
    nEP = 2N
    rate_update_weight = 0.05
    nus = rate_update_weight * N # maximal number of subproblems needed to be adjusted
end

name(_::MOEADUR) = "MOEA/D-UR"

mutable struct MOEADURState <: AbstractOptimizerState
    iteration::Any
    fcalls::Any
    start_time::Any
    stop_time::Any

    N::Any                           # population size
    population::Population           # population
    W::Any                           # weight vectors
    B::Any                           # neighbours of each solution
    z::Any                           # the best value found so far
    EP::Any                          # external population
end
pfront(s::MOEADURState) = get_non_dominated_solutions(s.population)
value(s::MOEADURState) = objectives(pfront(s))
minimizer(s::MOEADURState) = objectives(pfront(s))
copy(s::MOEADURState) = MOEADURState(
    s.iteration,
    s.fcalls,
    s.start_time,
    s.stop_time,
    s.N,
    copy(s.population),
    copy(s.W),
    copy(s.B),
    copy(s.z),
    copy(s.EP),
)

function initial_state(algorithm::MOEADUR, problem, options)
    population = [Individual(rand(problem.D)) for i = 1:algorithm.N]

    fcalls = evaluate!(problem, population)

    # Generate weight vectors
    W, N = uniformly_randomly_points(problem.M, algorithm.N)
    W = transform_weights(W)

    # Set the neighbours of each weight vector
    B = pairwise(Euclidean(), W, W, dims = 1)
    B = sortperm(B)
    B = B[:, 2:algorithm.T]

    # Set reference point
    z = Matrix(minimum(objectives(population), dims = 1)')

    EP = Individual[]

    return MOEADURState(
        0,
        fcalls,
        0,
        0,
        N,
        population,
        W,
        B,
        z,
        EP
    )
end

function update_state!(algorithm::MOEADUR, state, problem, options)
    offspring = copy(state.population)
    O = Individual[]

    for i = 1:state.N
        if rand() < algorithm.delta
            E = state.B[i, :]
        else
            E =
                i == 1 ? (2:state.N) :
                1 == state.N ? (1:state.N-1) : vcat(1:i-1, i+1:state.N)
        end

        # Reproduction
        xi = shuffle(options.rng, E)[1:2]
        y = evolute1(state, problem, xi, rng=options.rng)
        push!(O, y)

        # Update reference point
        state.z = minimum(hcat(state.z, objectives(y)), dims = 2)

        # Update of neighboring solutions
        update = 0
        while update < algorithm.nr && !isempty(E)
            j = rand(options.rng, E)
            x = state.population[j]
            g_x, _ = findmax(state.W[j, :] .* abs.(objectives(x) - state.z))
            g_y, _ = findmax(state.W[j, :] .* abs.(objectives(y) - state.z))
            if g_y <= g_x
                offspring[j] = y
                update += 1
            end
            E = filter(i -> i != j, E)
        end
    end

    state.population = offspring

    GMax = problem.maxFE / algorithm.N

    if state.iteration <= GMax * 0.9
        # println("Update EP on generation $(state.iteration)")
        # for o = O
        #     E = false
        #     Qi = []
        #     for qi = 1:length(state.EP)
        #         q = state.EP[qi]
        #         d = dominate(q, o)
        #         if d == 1
        #             E = true
        #         elseif d == -1
        #             push!(Qi, qi)
        #         end
        #     end
        #     if !E
        #         deleteat!(state.EP, Qi)
        #         push!(state.EP, o)
        #     end
        # end
        state.EP = [state.EP; offspring]
        F = nondominatedsort1(state.EP)
        state.EP = state.EP[F]

        while length(state.EP) > algorithm.nEP
            dis = pairwise(Euclidean(), objectives(state.EP), objectives(state.EP), dims = 1)
            dis[diagind(dis)] .= Inf
            sl = sparsity_level.(eachrow(dis), problem.M)
            worst = sortperm(sl)[end]
            deleteat!(state.EP, worst)
        end
    end

    if mod(state.iteration, GMax * algorithm.rate_update_weight) == 0 && state.iteration <= GMax * 0.9
        # println("Update weight vectors on generation $(state.iteration)")
        adjust = 0

        # Delete overcrowded subproblems
        while adjust < algorithm.nus
            dis = pairwise(Euclidean(), objectives(state.population), objectives(state.population), dims = 1)
            dis[diagind(dis)] .= Inf
            sparsity_level = map(ds -> prod(sort(ds)[1:problem.M]), eachrow(remove_diagonal(dis)))
            worst = sortperm(sparsity_level)[1]
            deleteat!(state.population, worst)
            state.W = state.W[Not(worst), :]
            adjust += 1
        end

        # Add new subproblems
        while adjust > 0
            # Generate a new subproblem using the individual which has the
            # largest sparsity level
            combine = [state.population; state.EP]
            dis = pairwise(Euclidean(), objectives(combine), objectives(combine), dims = 1)
            sl = sparsity_level.(eachrow(dis), problem.M)
            largest_sparsity = sortperm(sl)[end]
            ind = combine[largest_sparsity]
            λ = zeros(problem.M)
            for m = 1:problem.M
                λ[m] = 1 / objectives(ind)[m] - state.z[m] / sum(1 ./ objectives(ind) - state.z)
            end
            push!(state.population, ind)
            state.W = [state.W; λ']
            adjust -= 1
        end

        # Update the neighors of each weight vector
        state.B = pairwise(Euclidean(), state.W, state.W, dims = 1)
        state.B = sortperm(state.B)
        state.B = state.B[:, 2:algorithm.T]
    end

    return false
end
