Base.@kwdef struct MOEADDRA <: AbstractAlgorithm
    N::Int = 100            # number of the subproblems
    T::Int = ceil(N/10)     # number of the weight vectors in the neighborhoor of each weight vector
    F::Real = 0.5
    delta = 0.9
    nr = ceil(N/100)
end

mutable struct MOEADDRAState <: AbstractOptimizerState
    iteration
    fcalls
    start_time
    stop_time

    N
    population                  # population
    W                           # weight vectors
    B                           # neighbours of each solution
    z                           # the best value found so far
    EP                          # nondominated solutions
    I
    Pi                          # utility for each subproblem
    lastObjectives              # snapshot of objective values before update π
end
pfront(s::MOEADDRAState) = s.population
value(s::MOEADDRAState) = objectives(pfront(s))
minimizer(s::MOEADDRAState) = objectives(pfront(s))
copy(s::MOEADDRAState) = MOEADDRAState(s.iteration, s.fcalls, s.start_time, s.stop_time,
                                 copy(s.N),
                                 copy(s.population),
                                 copy(s.W),
                                 copy(s.B),
                                 copy(s.z),
                                 copy(s.EP),
                                 copy(s.I),
                                 copy(s.Pi),
                                 copy(s.lastObjectives))

function initial_state(algorithm::MOEADDRA, problem, options)
    population = [Individual(rand(problem.D)) for i in 1:algorithm.N]

    fcalls = evaluate!(problem, population)

    # Generate weight vectors
    W, N = NBI(population_size(algorithm), problem.M)

    population = population[1:N]

    # Set the neighbours of each weight vector
    B = pairwise(Euclidean(), W, W, dims=1)
    B = sortperm(B)
    B = B[:,2:algorithm.T]

    # Set reference point
    z = Matrix(minimum(objectives(population), dims=1)')

    EP = Individual[]

    Pi = ones(length(population))

    return MOEADDRAState(0, fcalls, 0, 0, N, population, W, B, z, EP, nothing, Pi, copy(objectives(population)))
end

function update_state!(
    algorithm::MOEADDRA,
    state,
    problem,
    options)

    I = []
    for i in 1:problem.M
        idx = findfirst(isequal(state.z[i]), objectives(state.population)[:,i])
        if idx !== nothing
            push!(I, idx)
        end
    end
    for i in tournament(10,select = argmax)(state.Pi, Int(round(length(state.population)/5))-length(I))
        push!(I, i)
    end

    for i in I
        if rand() < algorithm.delta
            P = state.B[i,:]
        else
            P = i==1 ? (2:state.N) : 1==state.N ? (1:state.N-1) : vcat(1:i-1, i+1:state.N)
        end

        x = state.population[i]

        # Reproduction
        xi = rand(options.rng, P, 2)
        xs = state.population[xi]
        y = copy(x)
        differentiation!(y, xs; F=algorithm.F)

        # Improvement
        y = PLM(lower = lower(problem), upper = upper(problem))(y, rng = options.rng)

        apply!(problem.constraints, variables(y))
        evaluate!(state, problem, y)

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
    end

    # Update utility for each subproblem
    if mod(state.iteration, 50) == 0
        delta = (state.lastObjectives .- objectives(state.population)) ./ state.lastObjectives
        for (i, ds) in enumerate(eachrow(delta))
            for d in ds
                if d > 0.001
                    state.Pi[i] = 1
                else
                    state.Pi[i] = (0.95+0.05*d/0.001)*state.Pi[i]
                end
            end
        end
        state.lastObjectives = objectives(state.population)
    end

    return false
end
