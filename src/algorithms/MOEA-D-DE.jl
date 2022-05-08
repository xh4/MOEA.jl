Base.@kwdef struct MOEADDE <: AbstractOptimizer
    # number of the subproblems
    N::Int = 100
    # number of the weight vectors in the neighborhoor of each weight vector
    T::Int = ceil(N/10)
    crossover = BINX(0.5)
    mutation = PLM()
    metrics::ConvergenceMetrics = ConvergenceMetric[GD(), GD(true)]
end

population_size(method::MOEADDE) = method.N
default_options(method::MOEADDE) = (iterations = 1000,)
summary(m::MOEADDE) =
    "MOEA/D-DE[P=$(population_size(m)),]"
show(io::IO, m::MOEADDE) = print(io, summary(m))

mutable struct MOEADDEState <: AbstractOptimizerState
    iteration
    start_time
    stop_time
    metrics::ConvergenceMetrics # collection of convergence metrics

    population                  # population
    m                           # number of objectives
    W                           # weight vectors
    B                           # neighbours of each solution
    z                           # the best value found so far
    EP                          # nondominated solutions
end
pfront(s::MOEADDEState) = s.EP
value(s::MOEADDEState) = objectives(pfront(s))
minimizer(s::MOEADDEState) = objectives(pfront(s))
copy(s::MOEADDEState) = MOEADDEState(s.iteration, s.start_time, s.stop_time, copy(s.metrics),
                                 copy(s.population),
                                 copy(s.m),
                                 copy(s.W),
                                 copy(s.B),
                                 copy(s.z),
                                 copy(s.EP))

function initial_state(method::MOEADDEState, objective, population, options)
    value!(objective, population)

    m = length(objectives(first(population)))

    # Generate weight vectors
    W = uniform_point(population_size(method), m)

    # Set the neighbours of each weight vector
    B = pairwise(Euclidean(), W, W, dims=1)
    B = sortperm(B)
    B = B[:,2:method.T]

    # Set reference point
    z = Matrix(minimum(objectives(population), dims=1)')

    EP = Individual[]
    return MOEADDEState(0, 0, 0, copy(method.metrics), population, m, W, B, z, EP)
end

function update_state!(
    method::MOEADDE,
    state,
    objective,
    constraints,
    options
)
    for i = 1:method.N
        x = state.population[i]

        # Reproduction
        xi = rand(options.rng, state.B[i,:], 2)
        xs = state.population[xi]
        y,_ = method.crossover(xs[1], xs[2], rng=options.rng)

        # Improvement
        y = method.mutation(y)

        apply!(constraints, variables(y))
        value!(objective, y)

        # Update reference point
        for j = 1:state.m
            if state.z[j] > objectives(y)[j]
                state.z[j] = objectives(y)[j]
            end
        end

        # Update of neighboring solutions
        for j = state.B[i,:]
            x = state.population[j]
            g_x, _ = findmax(state.W[j,:] .* abs.(objectives(x) - state.z))
            g_y, _ = findmax(state.W[j,:] .* abs.(objectives(y) - state.z))
            if g_y <= g_x
                state.population[j] = y
            end
        end

        # Update of EP
        deleteat!(state.EP, findall((p) -> dominate(y, p) == 1, state.EP))
        if isempty(findall((p) -> dominate(p, y) == 1, state.EP))
            pushfirst!(state.EP, y)
        end
    end

    return false
end
