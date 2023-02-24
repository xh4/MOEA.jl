Base.@kwdef struct NSGAIII <: AbstractAlgorithm
    N::Int = 100
end

name(_::NSGAIII) = "NSGA-III"

mutable struct NSGAIIIState <: AbstractOptimizerState
    iteration
    fcalls
    start_time
    stop_time

    N
    population
    Z
    Zmin
end
pfront(s::NSGAIIIState) = get_non_dominated_solutions(s.population)
value(s::NSGAIIIState) = objectives(s.population)
minimizer(s::NSGAIIIState) = variables(s.population)
copy(s::NSGAIIIState) = NSGAIIIState(s.iteration, s.fcalls, s.start_time, s.stop_time,
                                 s.N,
                                 copy(s.population),
                                 copy(s.Z),
                                 copy(s.Zmin))

function initial_state(algorithm::NSGAIII, problem, options)
    population = [Individual(rand(problem.D)) for i in 1:algorithm.N]
    Z, N = TwoLayer(algorithm.N, problem.M)
    fcalls = evaluate!(problem, population)
    Zmin = minimum(objectives(population), dims=1)
    return NSGAIIIState(0, fcalls, 0, 0, N, population, Z, Zmin)
end

function update_state!(
    algorithm::NSGAIII,
    state,
    problem,
    options
)
    parents = state.population
    selected = shuffle(1:state.N)
    offspring = evolute(state, problem, selected, rng=options.rng)
    combine = [parents; offspring]
    snapshot = copy(combine)

    S = []
    Fl = []
    F = nondominatedsort(combine)
    for f in F
        if length(f) + length(S) > algorithm.N
            Fl = f
            break
        else
            S = [S; f]
        end
    end

    if length(S) == algorithm.N
        state.population = snapshot[S]
    else
        state.population = snapshot[S]
        K = algorithm.N - length(S)
        N, M = size(objectives(combine))
        
        # Normalize
        state.Zmin = minimum([state.Zmin; objectives(combine)], dims=1)
        for ind in combine
            ind.objectives = (ind.objectives .- state.Zmin')[:,1]
        end
        extreme = fill(0, M)
        w = zeros(M, M) .+ 1e-6 + 1.0I
        for i = 1:M
            _, ci = findmin(maximum(objectives(combine) ./ repeat(w[i,:]', N), dims=2))
            extreme[i] = ci[1]
        end
        if LinearAlgebra.rank(objectives(combine[extreme])) == problem.M
            hyperplane = objectives(combine[extreme]) \ ones(M)
            a = 1 ./ hyperplane
        else
            a = maximum(objectives(combine), dims=1)
        end
        for ind = combine
            ind.objectives = (ind.objectives ./ a)[:,1]
        end

        # Associate
        cosine = 1 .- pairwise(CosineDist(), objectives(combine), state.Z, dims=1)
        distance = repeat(sqrt.(sum(objectives(combine).^2, dims=2)),1,size(state.Z)[1]) .* sqrt.(1 .- cosine.^2)
        d, pi = findmin(distance', dims=1)

        # Niching
        hist = fit(Histogram, reshape([ci[1] for ci in pi[S]], length(pi[S])), 1:size(state.Z)[1]+1)
        rho = hist.weights
        include = fill(true, length(rho))
        k = 0
        while k < K
            min = minimum(rho[include])
            Jmin = findall(identity, (rho .== min) .& include)
            j = shuffle(Jmin)[1]
            assoc = map(ci -> ci[2], filter(ci -> ci[1] == j, pi[Fl]))
            if length(assoc) > 0
                if rho[j] == 0
                    d = distance[assoc, j]
                    _, i = findmin(d)
                else
                    i = Int(ceil(rand()*length(Fl)))
                end
                push!(state.population, snapshot[Fl[i]])
                rho[j] += 1
                deleteat!(Fl, i)
                k += 1
            else
                include[j] = false
            end
        end
    end
    state.N = length(state.population)

    return false
end
