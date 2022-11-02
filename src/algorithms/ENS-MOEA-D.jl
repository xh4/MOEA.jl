Base.@kwdef struct ENSMOEAD <: AbstractOptimizer
    # number of the subproblems
    N::Int = 100
    # number of the weight vectors in the neighborhoor of each weight vector
    T::Int = ceil(N/10)
    F::Real = 1.0
    delta = 0.9
    nr = 1
    crossover = BINX(0.5)
    mutation = PLM()
    NS = 25:25:100
    metrics::ConvergenceMetrics = ConvergenceMetric[]
end

population_size(method::ENSMOEAD) = method.N
default_options(method::ENSMOEAD) = (iterations = 500,)
summary(m::ENSMOEAD) =
    "ENS-MOEA/D[P=$(population_size(m)),]"
show(io::IO, m::ENSMOEAD) = print(io, summary(m))

mutable struct ENSMOEADState <: AbstractOptimizerState
    iteration
    fcalls
    start_time
    stop_time
    metrics::ConvergenceMetrics # collection of convergence metrics

    population                  # population
    m                           # number of objectives
    W                           # weight vectors
    B                           # neighbours of each solution
    z                           # the best value found so far
    EP                          # nondominated solutions
    I
    Pi                          # utility for each subproblem
    lastObjectives              # snapshot of objective values before update π
end
pfront(s::ENSMOEADState) = s.EP
value(s::ENSMOEADState) = objectives(pfront(s))
minimizer(s::ENSMOEADState) = objectives(pfront(s))
copy(s::ENSMOEADState) = ENSMOEADState(s.iteration, s.fcalls, s.start_time, s.stop_time, copy(s.metrics),
                                 copy(s.population),
                                 copy(s.m),
                                 copy(s.W),
                                 copy(s.B),
                                 copy(s.z),
                                 copy(s.EP),
                                 copy(s.I),
                                 copy(s.Pi),
                                 copy(s.lastObjectives))

function initial_state(method::ENSMOEAD, objective, population, options)
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

    Pi = ones(length(population))

    return ENSMOEADState(0, 0, 0, 0, copy(method.metrics), population, m, W, B, z, EP, nothing, Pi, copy(objectives(population)))
end

function update_state!(
    method::ENSMOEAD,
    state,
    objective,
    constraints,
    options)

    I = []
    for i in 1:state.m
        idx = findfirst(isequal(state.z[i]), objectives(state.population)[:,i])
        if idx !== nothing
            push!(I, idx)
        end
    end
    for i in tournament(10,select = argmax)(state.Pi, Int(round(length(state.population)/5))-state.m)
        push!(I, i)
    end

    for i in I
        if rand() < method.delta
            P = state.B[i,:]
        else
            P = i==1 ? (2:method.N) : 1==method.N ? (1:method.N-1) : vcat(1:i-1, i+1:method.N)
        end

        x = state.population[i]

        # Reproduction
        xi = rand(options.rng, P, 2)
        xs = state.population[xi]
        y = copy(x)
        differentiation!(y, xs; F=method.F)

        # Improvement
        y,_ = method.crossover(y, x, rng=options.rng)

        apply!(constraints, variables(y))
        value!(objective, y)

        # Update reference point
        for j = 1:state.m
            if state.z[j] > objectives(y)[j]
                state.z[j] = objectives(y)[j]
            end
        end

        # Update of neighboring solutions
        c = 0
        while c < method.nr && !isempty(P)
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
        deleteat!(state.EP, findall((p) -> dominate(y, p) == 1, state.EP))
        if isempty(findall((p) -> dominate(p, y) == 1, state.EP))
            pushfirst!(state.EP, y)
        end

        # Update utility for each subproblem
        if state.iteration % 50 == 0
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
    end

    return false
end
