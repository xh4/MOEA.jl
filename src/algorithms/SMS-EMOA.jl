Base.@kwdef struct SMSEMOA <: AbstractAlgorithm
    N::Int = 100          # number of the subproblems
end

name(_::SMSEMOA) = "SMS-EMOA"

Base.@kwdef mutable struct SMSEMOAIndividual <: AbstractIndividual
    variables
    objectives
    rank
    contribution

    SMSEMOAIndividual(var) = new(var, nothing, nothing, nothing)
    SMSEMOAIndividual(var, obj) = new(var, obj, nothing, nothing)
    SMSEMOAIndividual(var, obj, rank, contribution) = new(var, obj, rank, contribution)
end
variables(i::SMSEMOAIndividual) = i.variables
objectives(i::SMSEMOAIndividual) = i.objectives
rank(i::SMSEMOAIndividual) = i.rank
contribution(i::SMSEMOAIndividual) = i.contribution
copy(i::SMSEMOAIndividual) = SMSEMOAIndividual(copy(i.variables), copy(i.objectives), i.rank, i.contribution)
SMSEMOAIndividual(i::Individual) = SMSEMOAIndividual(copy(i.variables), copy(i.objectives), nothing, nothing)

Base.@kwdef mutable struct SMSEMOAState <: AbstractOptimizerState
    iteration
    fcalls
    start_time
    stop_time

    population                  # population
    pfront                      # individuals of the first Pareto front
end
pfront(s::SMSEMOAState) = s.pfront
value(s::SMSEMOAState) = objectives(pfront(s))
minimizer(s::SMSEMOAState) = objectives(pfront(s))
copy(s::SMSEMOAState) = SMSEMOAState(s.iteration, s.fcalls, s.start_time, s.stop_time, copy(s.population), copy(s.pfront))

function initial_state(algorithm::SMSEMOA, problem, options)
    population = map(SMSEMOAIndividual, initial_population(algorithm, problem))
    fcalls = evaluate!(problem, population)

    return SMSEMOAState(0, fcalls, 0, 0, population, pfront)
end

function update_state!(
    algorithm::SMSEMOA,
    state,
    problem,
    options
)
    population = state.population

    for i = 1:algorithm.N
        offspring = evolute1(state, problem, rand(options.rng, 1:algorithm.N, 2), rng=options.rng)
        push!(population, offspring)

        nondominatedsort!(population)
        F = calcF(population)
        R = F[end]
        if (length(R)) == 1
            r = R[1]
        else
            update_contribution!(population, population[R])
            r = R[argmin(contribution.(population[R]))]
        end
        deleteat!(population, r)

        state.pfront = filter(i -> rank(i) == 1, population)
    end
    return false
end

function update_contribution!(population, last_front, n_samples = 10_000)
    # reset contribution
    for sol in population
        sol.contribution = Inf
    end

    ΔS = compute_contribution(last_front, n_samples)
    for i in eachindex(ΔS)
        last_front[i].contribution = ΔS[i]
    end
end

function compute_contribution(population, n_samples = 10_000)
    if isempty(population)
        # nothing to do
        return 
    end

    M = length(objectives(population[1]))

    # objetive function values
    Fs = objectives(population)
    N = size(Fs, 1)
    ΔS = fill(Inf, N)

    if M == 2
        rank = sortperm(view(Fs, :, 1))

        for i = 2 : N-1
            ΔS[rank[i]] = (Fs[rank[i+1],1] .- Fs[rank[i],1]) .* (Fs[rank[i-1],2] .- Fs[rank[i],2])
        end
    elseif N > 1
        ΔS = calculate_hv(Fs, nadir(Matrix(objectives(population)))*1.1, 1, n_samples)
    end

    return ΔS
end

function calculate_hv(points, bounds, k, n_sample)

    # This function is modified from the code in
    # http://www.tik.ee.ethz.ch/sop/download/supplementary/hype/

    N, M = size(points)
    if M > 2
        # Use the estimated method for three or more objectives
        alpha = zeros(1,N); 
        for i = 1:k 
            J = 1:i-1
            alpha[i] = prod((k .- J)./(N .- J ))./i
        end

        f_min = ideal(Matrix(points))
        S = f_min' .+ (bounds - f_min)' .* rand(n_sample, M)

        PdS  = zeros(Bool, N,n_sample)
        dS   = zeros(Int, n_sample)
        for i = 1:N
            x = sum(points[i,:]' .- S .<= 0, dims = 2) .== M
            mask = view(x, :,1)
            PdS[i, mask] .= true
            dS[mask] = dS[mask] .+ 1
        end

        F = zeros(N)
        for i = 1:N
            mask = view(dS, view(PdS, i,:))
            F[i] = sum(alpha[mask])
        end

        # ΔS
        F .* prod(bounds - f_min) / n_sample

    else
        # Use the accurate method for two objectives
        pvec  = 1:size(points,1)
        alpha = zeros(1,k)
        for i = 1 : k 
            j = 1 : i-1
            alpha[i] = prod((k.-j)./(N.-j))./i
        end
        hypesub(N,points,M,bounds,pvec,alpha,k)
    end
end

function hypesub(l,A,M,bounds,pvec,alpha,k)
    # The recursive function for the accurate method
    h = zeros(1,l)
    #[S,mask] = sortrows(A,M)
    mask = sortperm(view(A, :, M))
    S = A[mask,:]
    pvec  = pvec[mask]
    for i = 1 : size(S,1) 
        if i < size(S,1) 
            extrusion = S[i+1,M] - S[i,M]
        else
            extrusion = bounds[M] - S[i,M]
        end
        if M == 1
            if i > k
                break
            end
            if all(isone(alpha .>= 0))
                h[pvec[1:i]] = h[pvec[1:i]] .+ extrusion*alpha[i]
            end
        elseif extrusion > 0
            h += extrusion*hypesub(l,S[1:i,:],M-1,bounds,pvec[1:i],alpha,k)
        end
    end
    return h
end