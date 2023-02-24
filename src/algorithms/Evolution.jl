function evolute(state, problem, selected; rng::AbstractRNG = Random.default_rng())
    parents = state.population
    offspring = copy(parents)

    if encoding(problem) == :real
        recombine!(offspring, parents, selected, SBX(), rng=rng)
        mutate!(offspring, PLM(lower=lower(problem), upper=upper(problem)), rng=rng)
        map(x -> apply!(problem.constraints, variables(x)), offspring)
    elseif encoding(problem) == :binary
        recombine!(offspring, parents, selected, UX(), rng=rng)
        mutate!(offspring, Flip(), rng=rng)
    elseif encoding(problem) == :permutation
        recombine!(offspring, parents, selected, OX1(), rng=rng)
        mutate!(offspring, Swap2(), rng=rng)
    else
        error("unknown problem encoding $(encoding(problem))")
    end

    evaluate!(state, problem, offspring)
    offspring
end

function evolute1(state, problem, xi; rng::AbstractRNG = Random.default_rng())
    xs = state.population[xi]
    
    if encoding(problem) == :real
        y,_ = SBX()(xs[1], xs[2], rng=rng)
        y = PLM(lower=lower(problem), upper=upper(problem))(y, rng=rng)
        apply!(problem.constraints, variables(y))
    elseif encoding(problem) == :binary
        y,_ = UX()(xs[1], xs[2], rng=rng)
        y = Flip()(y, rng=rng)
    elseif encoding(problem) == :permutation
        y,_ = OX1()(xs[1], xs[2], rng=rng)
        y = Swap2()(y, rng=rng)
    else
        error("unknown problem encoding $(encoding(problem))")
    end
    
    evaluate!(state, problem, y)
    y
end

function evolute1_DE(state, problem, xi; F::Real = 0.5, rng::AbstractRNG = Random.default_rng())
    xs = state.population[xi]
    y = copy(xs[1])
    differentiation!(y, xs; F=F)
    y = PLM(lower=lower(problem), upper=upper(problem))(y, rng=rng)
    apply!(problem.constraints, variables(y))
    evaluate!(state, problem, y)
    y
end