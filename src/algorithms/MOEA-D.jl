Base.@kwdef struct MOEAD <: AbstractOptimizer
    # number of the subproblems
    N::Int = 100
    # number of the weight vectors in the neighborhoor of each weight vector
    T::Int = ceil(N/10)
    crossover = BINX(0.5)
    mutation = PLM()
    metrics::ConvergenceMetrics = ConvergenceMetric[GD(), GD(true)]
end

population_size(method::MOEAD) = method.N
default_options(method::MOEAD) = (iterations = 1000,)
summary(m::MOEAD) =
    "MOEA/D[P=$(population_size(m)),]"
show(io::IO, m::MOEAD) = print(io, summary(m))

mutable struct MOEADState <: AbstractOptimizerState
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
pfront(s::MOEADState) = s.EP
value(s::MOEADState) = objectives(pfront(s))
minimizer(s::MOEADState) = objectives(pfront(s))
copy(s::MOEADState) = MOEADState(s.iteration, s.start_time, s.stop_time, copy(s.metrics),
                                 copy(s.population),
                                 copy(s.m),
                                 copy(s.W),
                                 copy(s.B),
                                 copy(s.z),
                                 copy(s.EP))

function initial_state(method::MOEAD, objective, population, options)
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
    return MOEADState(0, 0, 0, copy(method.metrics), population, m, W, B, z, EP)
end

function update_state!(
    method::MOEAD,
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

# 参考 https://github.com/stevengj/Sobol.jl ？
# 这里使用 Das and Dennis's method
function uniform_point(N, M)
    # gen_ref_dirs 会生成 N+1 个...，需修复
    reduce(vcat, gen_ref_dirs(M, N-1)')
    # gen_ref_dirs(M, N-1)
end

function gen_ref_dirs(dimension, n_paritions)
    gen_weights(dimension, n_paritions)
end

function gen_weights(a, b)
    nobj = a;
    H    = b;
    a    = zeros(nobj);
    d    = H;
    w    = [];
    produce_weight!(a, 1, d, H, nobj, w)
    return Array.(w)
end

function  produce_weight!(a, i, d, H, nobj, w)
    for k=0:d
        if i<nobj
            a[i] = k;
            d2   = d - k;
            produce_weight!(a, i+1, d2, H, nobj, w);
        else
            a[i] = d;
            push!(w, a/H)
            break;
        end
    end
end

function Base.sortperm(M::Matrix)
    reduce(vcat, map(sortperm, eachrow(M))')
end

function Base.sortperm!(M::Matrix)
    error("Not implemented")
end
