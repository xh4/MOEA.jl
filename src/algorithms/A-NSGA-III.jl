Base.@kwdef struct ANSGAIII <: AbstractAlgorithm
    N::Int = 100
end

name(_::ANSGAIII) = "A-NSGA-III"

mutable struct ANSGAIIIState <: AbstractOptimizerState
    iteration
    fcalls
    start_time
    stop_time

    N
    population
    Z
    Zmin
    Zinterval
end
pfront(s::ANSGAIIIState) = get_non_dominated_solutions(s.population)
value(s::ANSGAIIIState) = objectives(s.population)
minimizer(s::ANSGAIIIState) = variables(s.population)
copy(s::ANSGAIIIState) = ANSGAIIIState(s.iteration, s.fcalls, s.start_time, s.stop_time,
                                 s.N,
                                 copy(s.population),
                                 copy(s.Z),
                                 copy(s.Zmin),
                                 copy(s.Zinterval))

function initial_state(algorithm::ANSGAIII, problem, options)
    population = initial_population(algorithm, problem)
    Z, N = TwoLayer(algorithm.N, problem.M)
    fcalls = evaluate!(problem, population)
    Zmin = minimum(objectives(population), dims=1)
    Zinterval = Z[1,end] - Z[2,end]
    return ANSGAIIIState(0, fcalls, 0, 0, N, population, Z, Zmin, Zinterval)
end

function update_state!(
    algorithm::ANSGAIII,
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

    function associate()
        cosine = 1 .- pairwise(CosineDist(), objectives(combine), state.Z, dims=1)
        distance = repeat(sqrt.(sum(objectives(combine).^2, dims=2)),1,size(state.Z)[1]) .* sqrt.(1 .- cosine.^2)
        d, pi = findmin(distance', dims=1)
        hist = fit(Histogram, reshape([ci[1] for ci in pi[S]], length(pi[S])), 1:size(state.Z)[1]+1)
        rho = hist.weights
        rho
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
        distance = nothing
        pi = nothing
        rho = associate()

        # Niching
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

    # Addition of reference points
    rho = associate()
    Zold = nothing
    while any(r->r>=2, rho) && state.Z != Zold
        Zold = state.Z
        for i = findall(r->r>=2, rho)
            p = repeat(state.Z[i,:], 1, problem.M)' .- state.Zinterval / problem.M
            p[diagind(p)] .+= state.Zinterval
            state.Z = [state.Z; p]
        end
        indics = map(ci->ci[1], findall(z->any(v->v<0, z), eachrow(state.Z)))
        state.Z = state.Z[Not(indics),:]
        tmp = round.(state.Z .* 1e4) ./ 1e4
        indics = unique(i -> tmp[i,:], 1:size(tmp)[1])
        state.Z = state.Z[indics,:]
        associate()
    end

    # Deletion of reference points
    indics = intersect(state.N+1:size(state.Z)[1], findall(r->r==0, rho))
    state.Z = state.Z[Not(indics),:]

    return false
end
