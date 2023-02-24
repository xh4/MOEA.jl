Base.@kwdef struct MOEADAWA <: AbstractAlgorithm
    N::Int = 100                 # number of the subproblems
    T::Int = ceil(N/10)          # number of the weight vectors in the neighborhoor of each weight vector
    delta = 0.9
    nr = ceil(N/100)
    nEP = ceil(N*1.5)
    wag = 100                    # iteration interval of utilizing the adaptive weight vector adjustment strategy
    rate_update_weight = 0.05
    nus = rate_update_weight * N # maximal number of subproblems needed to be adjusted
    rate_evol = 0.8              # percentage of computing resources assgned to adaptive weight vecgtor adjustment
end

name(_::MOEADAWA) = "MOEA/D-AWA"

mutable struct MOEADAWAState <: AbstractOptimizerState
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
    I::Any
    Pi::Any                          # utility for each subproblem
    lastObjectives::Any              # snapshot of objective values before update π
end
pfront(s::MOEADAWAState) = get_non_dominated_solutions(s.population)
value(s::MOEADAWAState) = objectives(pfront(s))
minimizer(s::MOEADAWAState) = objectives(pfront(s))
copy(s::MOEADAWAState) = MOEADAWAState(
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
    copy(s.I),
    copy(s.Pi),
    copy(s.lastObjectives),
)

function initial_state(algorithm::MOEADAWA, problem, options)
    population = initial_population(algorithm, problem)

    fcalls = evaluate!(problem, population)

    # Generate weight vectors
    W, N = TwoLayer(algorithm.N, problem.M)

    W = transform_weights(W)

    population = population[1:N]

    # Set the neighbours of each weight vector
    B = pairwise(Euclidean(), W, W, dims = 1)
    B = sortperm(B)
    B = B[:, 2:algorithm.T]

    # Set reference point
    z = Matrix(minimum(objectives(population), dims = 1)')

    EP = Individual[]

    Pi = ones(length(population))
    return MOEADAWAState(
        0,
        fcalls,
        0,
        0,
        N,
        population,
        W,
        B,
        z,
        EP,
        nothing,
        Pi,
        copy(objectives(population)),
    )
end

function update_state!(algorithm::MOEADAWA, state, problem, options)
    offspring = copy(state.population)

    I = []
    for i = 1:problem.M
        idx = findfirst(isequal(state.z[i]), objectives(state.population)[:, i])
        if idx !== nothing
            push!(I, idx)
        end
    end
    for i in tournament(10, select = argmax)(
        -state.Pi,
        Int(round(length(state.population) / 5)) - length(I),
    )
        push!(I, i)
    end

    for i in I
        if rand() < algorithm.delta
            P = state.B[i, :]
        else
            P =
                i == 1 ? (2:state.N) :
                1 == state.N ? (1:state.N-1) : vcat(1:i-1, i+1:state.N)
        end

        # Reproduction
        xi = shuffle(options.rng, P)[1:2]
        y = evolute1(state, problem, xi, rng=options.rng)

        # Update reference point
        state.z = minimum(hcat(state.z, objectives(y)), dims = 2)

        # Update of neighboring solutions
        c = 0
        while c < algorithm.nr && !isempty(P)
            j = rand(options.rng, P)
            x = state.population[j]
            g_x, _ = findmax(state.W[j, :] .* abs.(objectives(x) - state.z))
            g_y, _ = findmax(state.W[j, :] .* abs.(objectives(y) - state.z))
            if g_y <= g_x
                offspring[j] = y
                c += 1
            end
            P = filter(i -> i != j, P)
        end
    end

    # Update utility for each subproblem
    if mod(state.iteration, 50) == 0
        delta = (state.lastObjectives .- objectives(offspring)) ./ state.lastObjectives
        for (i, ds) in enumerate(eachrow(delta))
            for d in ds
                if d > 0.001
                    state.Pi[i] = 1
                else
                    state.Pi[i] = (0.95 + 0.05 * d / 0.001) * state.Pi[i]
                end
            end
        end
        state.lastObjectives = objectives(offspring)
    end

    state.population = offspring

    # Adaptive weight adjustment
    if state.fcalls >= algorithm.rate_evol * problem.maxFE && mod(state.iteration, algorithm.wag) == 0
        # println("Update EP on generation $(state.iteration)")
        if isempty(state.EP)
            state.EP = update_ep(state.population, algorithm.nEP, problem.M)
        else
            state.EP = update_ep([state.EP; state.population], algorithm.nEP, problem.M)
        end
        # println("Adjust weights on generation $(state.iteration) fcalls $(state.fcalls)")
        # Delete the overcrowded subproblems
        n = 0
        while n < algorithm.nus
            dis = pairwise(Euclidean(), objectives(state.population), objectives(state.population), dims = 1)
            sparsity_level = map(ds -> prod(sort(ds)[1:problem.M]), eachrow(remove_diagonal(dis)))
            worst = sortperm(sparsity_level)[1]
            deleteat!(state.population, worst)
            state.W = state.W[Not(worst), :]
            n += 1
        end
        
        # Add new subproblems
        # Remove individuals in EP which are dominated by individuals in evolutional population
        for i = 1:length(state.EP)
            for j = 1:length(state.population)
                if dominate(state.population[j], state.EP[i]) == 1
                    deleteat!(state.EP, i)
                    break
                end
            end
        end
        n = 0
        while n < algorithm.nus
            # Generate a new subproblem using the individual which has the
            # largest sparsity level
            combine = [state.population; state.EP]
            dis = pairwise(Euclidean(), objectives(combine), objectives(combine), dims = 1)
            sparsity_level = map(ds -> prod(sort(ds)[1:problem.M]), eachrow(remove_diagonal(dis)))
            largest_sparsity = sortperm(sparsity_level)[end]
            ind = combine[largest_sparsity]
            λ = zeros(problem.M)
            for m = 1:problem.M
                λ[m] = 1 / objectives(ind)[m] - state.z[m] / sum(1 ./ objectives(ind) - state.z)
            end
            push!(state.population, ind)
            state.W = [state.W; λ']
            n += 1
        end
        # Find the T closest weights and build the new B
        state.B = pairwise(Euclidean(), state.W, state.W, dims = 1)
        state.B = sortperm(state.B)
        state.B = state.B[:, 2:algorithm.T]
    end

    return false
end

function transform_weights(W)
    1.0 ./ W ./ repeat(sum(1.0 ./ W, dims = 2), 1, size(W, 2))
end

function update_ep(population, n, m)
    ndx = get_non_dominated_solutions_perm(population)
    while length(ndx) > n
        dis = pairwise(Euclidean(), objectives(population[ndx]), objectives(population[ndx]), dims = 1)
        sparsity_level = map(ds -> prod(sort(ds)[1:m]), eachrow(remove_diagonal(dis)))
        worst = sortperm(sparsity_level)[end]
        deleteat!(ndx, worst)
    end
    population[ndx]
end

function sparsity_level(dis, m)
    prod(sort(dis)[1:m])
end

function remove_diagonal(x)
    sz = LinearIndices(x)
    n = size(x, 1)
    k = [sz[i,i] for i in 1:n]
    b = collect(vec(x'))
    deleteat!(b, k)
    reshape(b, n - 1, n)'
end