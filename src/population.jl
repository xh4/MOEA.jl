const Population = Vector{Individual}

variable(pop::Population) = mapreduce(variable, hcat, pop)'

objective(pop::Population) = mapreduce(objective, hcat, pop)'

constraint(pop::Population) = map(constraint, pop)

copy(pop::Population) = map(copy, pop)

Population(method::M, individual::I; kwargs...) where {M<:Optimizer,I<:Individual} =
    [copy(individual) for i = 1:population_size(method)]

Population(
    method::M,
    variable::I;
    kwargs...,
) where {M<:Optimizer,I<:AbstractVector} =
    [Individual(variable) for i = 1:population_size(method)]

function Population(
    method::M,
    individuals::AbstractVector{I};
    kwargs...,
) where {M<:Optimizer,I<:AbstractVector}
    n = population_size(method)
    @assert length(individuals) == n "Size of initial population must be $n"
    individuals
end

Population(method::M, individualFunc::Function; kwargs...) where {M<:Optimizer} =
    [individualFunc() for i = 1:population_size(method)]

function Population(
    method::M,
    individual::I;
    kwargs...,
) where {M<:Optimizer,I<:AbstractMatrix}
    [copy(individual) for i = 1:population_size(method)]
end

function Population(
    method::M,
    bounds::ConstraintBounds;
    rng::AbstractRNG = Random.GLOBAL_RNG,
) where {M<:Optimizer}
    n = population_size(method)
    cn = nconstraints_x(bounds)
    indv = rand(rng, cn, n)
    if length(bounds.eqx) > 0
        indv[bounds.eqx, :] .= bounds.valx
    end
    T = eltype(bounds)
    rngs = Dict(i => zeros(T, 2) for i in bounds.ineqx)
    for (j, v, s) in zip(bounds.ineqx, bounds.bx, bounds.σx)
        rngs[j][1] -= s * v
        if s > 0
            rngs[j][2] += v
        end
    end
    for i = 1:n
        for (j, (r, l)) in rngs
            indv[j, i] *= r
            indv[j, i] += l
        end
    end
    Population(method, [collect(i) for i in eachcol(indv)])
end