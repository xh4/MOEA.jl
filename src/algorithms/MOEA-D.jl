Base.@kwdef struct MOEAD <: AbstractAlgorithm
    N::Int = 100          # number of the subproblems
    T::Int = ceil(N/10)   # number of the weight vectors in the neighborhoor of each weight vector
end

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
    population = [Individual(rand(problem.D)) for i in 1:algorithm.N]

    fcalls = evaluate!(problem, population)

    m = length(objectives(first(population)))

    # Generate weight vectors
    W, N = NBI(population_size(algorithm), m)

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
    algorithm::MOEAD,
    state,
    problem,
    options
)
    for i = 1:state.N
        x = state.population[i]

        # Reproduction
        xi = shuffle(options.rng, state.B[i,:])[1:2]
        xs = state.population[xi]
        y,_ = SBX()(xs[1], xs[2], rng=options.rng)

        # Improvement
        y = PLM(lower=lower(problem), upper=upper(problem))(y, rng=options.rng)

        apply!(problem.constraints, variables(y))

        evaluate!(state, problem, y)

        # Update reference point
        state.z = minimum(hcat(state.z, objectives(y)), dims=2)

        # Update of neighboring solutions
        for j = state.B[i,:]
            x = state.population[j]
            g_x = maximum(state.W[j,:] .* abs.(objectives(x) - state.z))
            g_y = maximum(state.W[j,:] .* abs.(objectives(y) - state.z))
            # println(g_y)
            if g_y <= g_x
                state.population[j] = y
            end
        end

        # Update of EP
        state.EP = state.population
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

function produce_weight!(a, i, d, H, nobj, w)
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
