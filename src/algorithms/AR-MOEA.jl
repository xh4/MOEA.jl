Base.@kwdef struct ARMOEA <: AbstractAlgorithm
    N::Int = 100
end

name(_::ARMOEA) = "AR-MOEA"

mutable struct ARMOEAState <: AbstractOptimizerState
    iteration
    fcalls
    start_time
    stop_time

    N
    population                  # population
    W
    Archive
    RefPoint
    Range
end
pfront(s::ARMOEAState) = get_non_dominated_solutions(s.population)
value(s::ARMOEAState) = objectives(s.population)
minimizer(s::ARMOEAState) = variables(s.population)
copy(s::ARMOEAState) = ARMOEAState(s.iteration, s.fcalls, s.start_time, s.stop_time,
                                 s.N,
                                 copy(s.population),
                                 copy(s.W),
                                 copy(s.Archive),
                                 copy(s.RefPoint),
                                 copy(s.Range))

function initial_state(algorithm::ARMOEA, problem, options)
    population = initial_population(algorithm, problem)
    fcalls = evaluate!(problem, population)
    
    W, N = TwoLayer(algorithm.N, problem.M)

    Archive = nothing
    RefPoint = nothing
    Range = nothing

    Archive, RefPoint, Range = ARMOEA_update_ref_point(population, W, Range)

    return ARMOEAState(0, fcalls, 0, 0, N, population, W, Archive, RefPoint, Range)
end

function update_state!(
    algorithm::ARMOEA,
    state,
    problem,
    options
)
    parents = state.population
    
    selected = ARMOEA_mating_selection(parents, state.RefPoint, state.Range)
    offspring = evolute(state, problem, selected, rng=options.rng)

    state.Archive, state.RefPoint, state.Range = ARMOEA_update_ref_point([state.Archive; offspring], state.W, state.Range)

    combine = [parents; offspring]
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
        state.population = combine[S]
    else
        K = algorithm.N - length(S)
        # Select K individuals from last Pareto front
        N = length(Fl)
        NR = size(state.RefPoint)[1]
        distance = ARMOEA_cal_distance(objectives(combine[Fl]) - repeat(state.Range[1,:]', N, 1), state.RefPoint)
        convergence = minimum(distance, dims=2)
        rank = reduce(hcat, map(sortperm, eachcol(distance)))
        dis = reduce(hcat, map(t -> distance[t[2],t[1]], enumerate(eachcol(rank))))
        remain = fill(true, N)
        while sum(remain) > K
            # println("Fl: $(length(Fl)), K: $(K), remain: $(sum(remain)), rank: $(size(rank)[1])")
            # Calculate the fitness of noncontributing solutions
            Non_contributing = copy(remain)
            Non_contributing[rank[1,:]] .= false
            METRIC = sum(dis[1,:]) + sum(convergence[Non_contributing])
            Metric = fill(Inf, N)
            Metric[Non_contributing] .= METRIC .- convergence[Non_contributing]
            # Calculate the fitness of contributing solutions
            for p = findall(remain .& .!Non_contributing)
                temp = rank[1,:] .== p
                non_contributing = fill(false, N)
                non_contributing[rank[2,temp]] .= true
                non_contributing = non_contributing .& Non_contributing
                Metric[p] = METRIC - sum(dis[1,temp]) + sum(dis[2,temp]) - sum(convergence[non_contributing])
            end
            # Delete the worst solution and update the variables
            _, del = findmin(Metric)
            temp = rank .!= del
            dis = reshape(dis[temp], sum(remain)-1, NR)
            rank = reshape(rank[temp], sum(remain)-1, NR)
            remain[del] = false
        end
        state.population = combine[[S; Fl[remain]]]
    end

    state.Range[2,:] .= maximum(objectives(state.population), dims=1)'
    state.Range[2, (state.Range[2,:] - state.Range[1,:]) .< 1e-6] .= 1

    return false
end

function ARMOEA_mating_selection(parents, ref_point, range; rng=Random.default_rng())
    N = length(parents)

    # Calculate the fitness of each feasible solution based on IGD-NS
    distance = ARMOEA_cal_distance(objectives(parents) - repeat(range[1,:]', N, 1), ref_point)
    convergence = minimum(distance, dims=2)
    rank = reduce(hcat, map(sortperm, eachcol(distance)))
    dis = reduce(hcat, map(t -> distance[t[2],t[1]], enumerate(eachcol(rank))))
    Non_contributing = fill(true, N)
    Non_contributing[rank[1,:]] .= false
    METRIC = sum(dis[1,:]) + sum(convergence[Non_contributing])
    fitness = fill(Inf, N)
    fitness[Non_contributing] = METRIC .- convergence[Non_contributing]
    for p = findall(.!Non_contributing)
        temp = rank[1,:] .== p
        non_contributing = fill(false, N)
        non_contributing[rank[2,temp]] .= true
        non_contributing = non_contributing .& Non_contributing
        fitness[p] = METRIC - sum(dis[1,temp]) + sum(dis[2,temp]) - sum(convergence[non_contributing])
    end
    Fitness = fitness

    tournament(2, select = twowaycomp)(hcat(fill(0, N), -Fitness)', N; rng=rng)
end

function ARMOEA_update_ref_point(archive, W, range)
    F1 = nondominatedsort1(archive)
    F1 = unique(i -> objectives(archive[i]), F1)
    archive = archive[F1]
    NA = length(F1)
    NW = size(W)[1]

    if !isnothing(range)
        range[1,:] .= minimum([range[1,:]'; objectives(archive)], dims=1)'
    elseif !isnothing(archive)
        range = [minimum(objectives(archive), dims=1);
                 maximum(objectives(archive), dims=1)]
    end

    if length(archive) <= 1
        ref_point = W
    else
        # Find contributing solutions and valid weight vectors
        t_archive = objectives(archive) - repeat(range[1,:]', NA)
        W = W .* repeat((range[2,:] - range[1,:])', NW)
        distance = ARMOEA_cal_distance(t_archive, W)
        nearest_p = map(ci -> ci[1], findmin(distance, dims=1)[2])[1,:]
        contributing_s = unique(nearest_p)
        nearest_w = map(ci -> ci[2], findmin(distance, dims=2)[2])[:,1]
        valid_w = unique(nearest_w[contributing_s])

        # Update archive
        choose = in.(collect(1:NA), [contributing_s])
        cosine = 1 .- pairwise(CosineDist(), t_archive, t_archive, dims=1)
        cosine[diagind(cosine)] .= 0
        while sum(choose) < min(3*NW, size(t_archive)[1])
            un_selected = findall(.!choose)
            x = map(ci -> ci[1], findmin(maximum(cosine[.!choose, choose], dims=2), dims=1)[2])[1,:]
            choose[un_selected[x]] .= 1
        end
        archive = archive[choose]
        t_archive = t_archive[choose,:]

        # Update reference points
        ref_point = [W[valid_w,:]; t_archive]
        choose = [fill(true, length(valid_w)); fill(false, size(t_archive)[1])]
        cosine = 1 .- pairwise(CosineDist(), ref_point, ref_point, dims=1)
        cosine[diagind(cosine)] .= 0
        while sum(choose) < min(NW, size(ref_point)[1])
            selected = findall(.!choose)
            x = map(ci -> ci[1], findmin(maximum(cosine[.!choose, choose], dims=2), dims=1)[2])[1,:]
            choose[selected[x]] .= true
        end
        ref_point = ref_point[choose,:] 
    end
    (archive, ref_point, range)
end

function ARMOEA_cal_distance(pop_obj, ref_point)
    N = size(pop_obj)[1]
    NR = size(ref_point)[1]
    pop_obj = max.(pop_obj, 1e-6)
    ref_point = max.(ref_point, 1e-6)

    cosine = 1 .- pairwise(CosineDist(), pop_obj, ref_point, dims=1)
    norm_r = sqrt.(sum(ref_point .^ 2, dims=2))
    norm_p = sqrt.(sum(pop_obj .^ 2, dims=2))
    d1 = repeat(norm_p, 1, NR) .* cosine
    d2 = repeat(norm_p, 1, NR) .* sqrt.(1 .- cosine .^ 2)
    nearest = map(ci -> ci[1], findmin(d2, dims=1)[2])[1,:]

    ref_point = ref_point .* repeat(d1[N .* collect(0:NR-1) + nearest] ./ norm_r, 1, size(ref_point)[2])

    distance = pairwise(Euclidean(), pop_obj, ref_point, dims=1)
    distance
end