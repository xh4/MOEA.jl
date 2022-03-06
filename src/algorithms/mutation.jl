# Mutation operators
# ==================


# Genetic mutations
# =================

# Binary mutations
# ---------------------

"""
    flip(recombinant)

Returns an in-place mutated binary `recombinant` with a bit flips at random positions.
"""
function flip(
    recombinant::T;
    rng::AbstractRNG = Random.default_rng(),
) where {T<:AbstractVector{Bool}}
    s = length(recombinant)
    pos = rand(rng, 1:s)
    recombinant[pos] = !recombinant[pos]
    return recombinant
end
flip(recombinant::Individual; rng::AbstractRNG = Random.default_rng(),) = 
    Individual(flip(variable(recombinant), rng=rng))

"""
    bitinversion(recombinant)

Returns an in-place mutated binary `recombinant` with its bits inverted.
"""
bitinversion(recombinant::T) where {T<:AbstractVector{Bool}} =
    map!(!, recombinant, recombinant)
bitinversion(recombinant::Individual) = 
    Individual(bitinversion(variable(recombinant)))

# Real-valued mutations
# ---------------------

"""
    uniform(r = 1.0)

Returns an in-place real valued mutation function that performs the uniform distributed mutation [^1].

The mutation operator randomly chooses a number ``z`` in from the uniform distribution on the interval ``[-r,r]``, the mutation range.
The mutated individual is given by

- ``x_i^\\prime = x_i + z_i``

"""
function uniform(r::Real = 1.0)
    function mutation(
        recombinant::T;
        rng::AbstractRNG = Random.default_rng(),
    ) where {T<:AbstractVector}
        d = length(recombinant)
        recombinant .+= 2r .* rand(rng, d) .- r
        return recombinant
    end
    function mutation(
        recombinant::T;
        rng::AbstractRNG = Random.default_rng(),
    ) where {T<:Individual}
        Individual(mutation(variable(recombinant), rng=rng))
    end
    return mutation
end

"""
    gaussian(σ = 1.0)

Returns an in-place real valued mutation function that performs the normal distributed mutation [^1].

The mutation operator randomly chooses a number ``z`` in from the normal distribution ``\\mathcal{N}(0,\\sigma)`` with standard deviation ``\\sigma``.
The mutated individual is given by

- ``x_i^\\prime = x_i + z_i``

"""
function gaussian(σ::Real = 1.0)
    function mutation(
        recombinant::T;
        rng::AbstractRNG = Random.default_rng(),
    ) where {T<:AbstractVector}
        d = length(recombinant)
        recombinant .+= σ .* randn(rng, d)
        return recombinant
    end
    function mutation(
        recombinant::T;
        rng::AbstractRNG = Random.default_rng(),
    ) where {T<:Individual}
        Individual(mutation(variable(recombinant), rng=rng))
    end
    return mutation
end

"""
    BGA(valrange, m = 20)

Returns an in-place real valued mutation function that performs the BGA mutation scheme with the mutation range `valrange` and the mutation probability `1/m` [^1].
"""
function BGA(valrange::Vector, m::Int = 20)
    prob = 1.0 / m
    function mutation(
        recombinant::T;
        rng::AbstractRNG = Random.default_rng(),
    ) where {T<:AbstractVector}
        d = length(recombinant)
        @assert length(valrange) == d "Range matrix must have $(d) columns"
        δ = zeros(m)
        for i = 1:length(recombinant)
            for j = 1:m
                δ[j] = (rand(rng) < prob) ? δ[j] = 2.0^(-j) : 0.0
            end
            if rand(rng, Bool)
                recombinant[i] += sum(δ) * valrange[i]
            else
                recombinant[i] -= sum(δ) * valrange[i]
            end
        end
        return recombinant
    end
    function mutation(
        recombinant::T;
        rng::AbstractRNG = Random.default_rng(),
    ) where {T<:Individual}
        Individual(mutation(variable(recombinant), rng=rng))
    end
    return mutation
end

"""
    PM(lower, upper, p = 2)

Returns an in-place real valued mutation function that performs the Power Mutation (PM)
scheme within `lower` and `upper` bound, and an index of
mutation `p`[^3].

*Note:* The implementation is a degenerate case of Mixed Integer Power Mutation ([`MIPM`](@ref))
"""
function PM(lower::Vector, upper::Vector, p::Float64 = 5.0) # index of distribution p
    return mipmmutation(lower, upper, p)
end

"""
    MIPM(lower, upper, p_real = 10, p_int = 4)

Returns an in-place real valued mutation function that performs the Mixed Integer Power Mutation (MI-PM) scheme within `lower` and `upper` bound, and an index of mutation `p_real` for real value and `p_int` for integer values[^4].
"""
function MIPM(
    lowerBounds::Vector,
    upperBounds::Vector,
    p_real::Float64 = 10.0,
    p_int::Float64 = 4.0,
) # index of distribution p
    return mipmmutation(lowerBounds, upperBounds, p_real, p_int)
end

function mipmmutation(
    lowerBounds::Vector,
    upperBounds::Vector,
    p_real::Float64,
    p_int::Union{Nothing,Float64} = nothing,
)
    function pm_mutation(rng, x, l, u, s, d)
        x̄ = d < rand(rng) ? x - s * (x - l) : x + s * (u - x)
        if isa(x, Integer)
            if isinteger(x̄)
                Int(x̄)
            else
                floor(Int, x̄) + (rand(rng) > 0.5)
            end
        else
            x̄
        end
    end
    function mutation(
        recombinant::T;
        rng::AbstractRNG = Random.default_rng(),
    ) where {T<:Vector}
        d = length(recombinant)
        @assert length(lowerBounds) == d "Bounds vector must have $(d) columns"
        @assert length(upperBounds) == d "Bounds vector must have $(d) columns"
        @assert !(p_int === nothing && any(isa(x, Integer) for x in recombinant)) "Need to set p_int for integer variables"
        u = rand(rng)
        P = (isa(x, Integer) ? p_int : p_real for x in recombinant)
        S = u .^ (1 ./ P) # random var following power distribution
        D = (recombinant - lowerBounds) ./ (upperBounds - lowerBounds)
        broadcast!(
            (x, l, u, s, d) -> pm_mutation(rng, x, l, u, s, d),
            recombinant,
            recombinant,
            lowerBounds,
            upperBounds,
            S,
            D,
        )
        return recombinant
    end
    function mutation(
        recombinant::T;
        rng::AbstractRNG = Random.default_rng(),
    ) where {T<:Individual}
        Individual(mutation(variable(recombinant), rng=rng))
    end
    return mutation
end

"""
    PLM(lower, upper, η = 2)

Returns an in-place real valued mutation function that performs the Polynomial Mutation (PLM) scheme
within `lower` and `upper` bounds, and a mutation distribution index `η`[^9].
"""
function PLM(Δ::Union{Real,Vector} = 1.0; η = 2, pm::Real = NaN) # index of distribution p
    function mutation(
        recombinant::T;
        rng::AbstractRNG = Random.default_rng(),
    ) where {T<:AbstractVector}
        d = length(recombinant)
        pm = isnan(pm) ? 1 / d : pm
        mask = rand(rng, d) .< pm
        u = rand(rng, d)
        mask_p = u .<= 0.5
        mask_n = u .> 0.5
        δpw = 1 / (η + 1)
        δ = similar(recombinant)
        δ[mask_p] .= (2 * u[mask_p]) .^ δpw .- 1
        δ[mask_n] .= 1 .- (2(1 .- u[mask_n])) .^ δpw
        recombinant[mask] .+= (Δ.*δ)[mask]
        return recombinant
    end
    function mutation(
        recombinant::T;
        rng::AbstractRNG = Random.default_rng(),
    ) where {T<:Individual}
        Individual(mutation(variable(recombinant), rng=rng))
    end
    return mutation
end
PLM(lower::Vector, upper::Vector; η::Real = 2, pm::Real = NaN) =
    PLM(upper - lower; η = η, pm = pm)


# Combinatorial mutations (applicable to binary vectors)
# ------------------------------------------------------

"""
    inversion(recombinant)

Returns an in-place mutated individual with a random arbitrary length segment of the genome in the reverse order.
"""
function inversion(
    recombinant::T;
    rng::AbstractRNG = Random.default_rng(),
) where {T<:AbstractVector}
    l = length(recombinant)
    from, to = randseg(rng, l)
    l = round(Int, (to - from) / 2)
    if from + 1 == to
        swap!(recombinant, from, to)
    else
        for i = 0:(l-1)
            swap!(recombinant, from + i, to - i)
        end
    end
    return recombinant
end
inversion(recombinant::T; rng::AbstractRNG = Random.default_rng(),) where {T<:Individual} = 
    Individual(inversion(variable(recombinant), rng = rng))

"""
    insertion(recombinant)

Returns an in-place mutated individual with an arbitrary element of the genome moved in a random position.
"""
function insertion(
    recombinant::T;
    rng::AbstractRNG = Random.default_rng(),
) where {T<:AbstractVector}
    l = length(recombinant)
    from, to = randseg(rng, l)
    val = recombinant[from]
    deleteat!(recombinant, from)
    return insert!(recombinant, to, val)
end
insertion(recombinant::T; rng::AbstractRNG = Random.default_rng(),) where {T<:Individual} = 
    Individual(insertion(variable(recombinant), rng = rng))

"""
    swap2(recombinant)

Returns an in-place mutated individual with a two random elements of the genome are swapped.
"""
function swap2(
    recombinant::T;
    rng::AbstractRNG = Random.default_rng(),
) where {T<:AbstractVector}
    l = length(recombinant)
    from, to = randseg(rng, l)
    swap!(recombinant, from, to)
    return recombinant
end
swap2(recombinant::T; rng::AbstractRNG = Random.default_rng(),) where {T<:Individual} = 
    Individual(swap2(variable(recombinant), rng = rng))

"""
    scramble(recombinant)

Returns an in-place mutated individual with elements, on a random arbitrary length segment of the genome, been scrambled.
"""
function scramble(
    recombinant::T;
    rng::AbstractRNG = Random.default_rng(),
) where {T<:AbstractVector}
    l = length(recombinant)
    from, to = randseg(rng, l)
    diff = to - from + 1
    if diff > 1
        patch = recombinant[from:to]
        idx = randperm(rng, diff)
        for i = 1:diff
            recombinant[from+i-1] = patch[idx[i]]
        end
    end
    return recombinant
end
scramble(recombinant::T; rng::AbstractRNG = Random.default_rng(),) where {T<:Individual} = 
    Individual(scramble(variable(recombinant), rng = rng))

"""
    shifting(recombinant)

Returns an in-place mutated individual with a random arbitrary length segment of the genome been shifted to an arbitrary position.
"""
function shifting(
    recombinant::T;
    rng::AbstractRNG = Random.default_rng(),
) where {T<:AbstractVector}
    l = length(recombinant)
    from, to, where = sort(rand(rng, 1:l, 3))
    patch = recombinant[from:to]
    diff = where - to
    if diff > 0
        # move values after tail of patch to the patch head position
        for i = 1:diff
            recombinant[from+i-1] = recombinant[to+i]
        end
        # place patch values in order
        start = from + diff
        for i = 1:length(patch)
            recombinant[start+i-1] = patch[i]
        end
    end
    return recombinant
end
shifting(recombinant::T; rng::AbstractRNG = Random.default_rng(),) where {T<:Individual} = 
    Individual(shifting(variable(recombinant), rng = rng))

"""
    replace(pool,[minchange=1])(recombinant)

Replacement mutation operator changes an arbitrary number, no smaller then `minchange`,
of elements in the individual by replacing them with elements from the predefined `pool` that are not in the individual.
"""
function replace(pool::Vector{P}; minchange = 1) where {P}
    function rplc(
        recombinant::T;
        rng::AbstractRNG = Random.default_rng(),
    ) where {T<:AbstractVector}
        l = length(recombinant)
        p = length(pool)
        # how many values to change
        nchg = max(minchange, rand(rng, 1:min(l, p - l)))
        # select new elements
        idxs = randperm(rng, p)
        new_vals = P[]
        for i in idxs
            if pool[i] ∉ recombinant
                push!(new_vals, pool[i])
            end
            length(new_vals) == nchg && break
        end
        # update arbitrary positions with new values
        new_idxs = randperm(rng, l)[1:nchg]
        recombinant[new_idxs] .= new_vals
        return recombinant
    end
    function rplc(
        recombinant::T;
        rng::AbstractRNG = Random.default_rng(),
    ) where {T<:Individual}
        Individual(rplc(variable(recombinant), rng = rng))
    end
    return rplc
end


# Differential Evolution
# ======================

"""
    differentiation(recombinant, mutators; F = 1.0)

Returns an in-place differently mutated individual ``x^\\prime`` from `recombinant` ``x``  by `mutators` ``\\{\\xi_1, \\ldots, \\xi_n \\}`` as follows

- ``x^\\prime = x + \\sum_{i=1}^{n/2} F (\\xi_{2i-1} - \\xi_{2i})``

"""
function differentiation(
    recombinant::T,
    mutators::AbstractVector{T};
    F::Real = 1.0,
) where {T<:AbstractVector}
    m = length(mutators)
    @assert m % 2 == 0 "Must be even number of target mutators"
    for i = 1:2:m
        recombinant .+= F .* (mutators[i] .- mutators[i+1])
    end
    return recombinant
end
differentiation(recombinant::T, mutators::AbstractVector{T}; F::Real = 1.0,) where {T<:Individual} = 
    Individual(differentiation(variable(recombinant), map(variable, mutators), F=F))
