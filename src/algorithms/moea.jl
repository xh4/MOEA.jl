"""
    dominate(p, q)

Returns `1` if `p` dominates `q`, `-1` if otherwise, and `0` if dominance cannot be determined.
"""
function dominate(p::T, q::T) where {T<:AbstractArray}
    ret = 0
    for (i, j) in zip(p, q)
        if i < j
            ret == -1 && return 0
            ret = 1
        elseif j < i
            ret == 1 && return 0
            ret = -1
        end
    end
    return ret
end
dominate(p::T, q::T) where {T<:Individual} = dominate(objectives(p), objectives(q))

"""
dominations(P::AbstractVector)

Returns a domination matrix of all elements in the input collection `P`.
"""
function dominations(P::AbstractVector{T}) where {T<:AbstractArray}
    l = length(P)
    D = zeros(Int8, l, l)
    for i = 1:l
        for j = (i+1):l
            D[i, j] = dominate(P[i], P[j])
            D[j, i] = -D[i, j]
        end
    end
    D
end
dominations(P::Population) = dominations(objectives(P))

# todo: 把 P 转置一下
"""
    nondominatedsort!(R, F)

Calculate fronts for fitness values `F`, and store ranks of the individuals into `R`.
"""
function nondominatedsort!(R, P)
    n = size(P, 2)
    @assert length(R) == n "Ranks must be defined for the whole population"

    Sₚ = Dict(i => Set() for i = 1:n)
    C = zeros(Int, n)

    # construct first front
    F = [Int[]]
    for i = 1:n
        for j = i+1:n
            r = dominate(view(P, :, i), view(P, :, j)) #M[i,j]
            if r == 1
                push!(Sₚ[i], j)
                C[j] += 1
            elseif r == -1
                push!(Sₚ[j], i)
                C[i] += 1
            end
        end
        if C[i] == 0
            R[i] = 1
            push!(F[1], i)
        end
    end

    # construct rest of the fronts
    while !isempty(last(F))
        Q = Int[]
        for i in last(F)
            for j in Sₚ[i]
                C[j] -= 1
                if C[j] == 0
                    push!(Q, j)
                    R[j] = length(F) + 1
                end
            end
        end
        push!(F, Q)
    end
    isempty(last(F)) && pop!(F)

    F #, R #, Sₚ
end
function nondominatedsort!(P::Population)
    R = fill(0, length(P))
    nondominatedsort!(R, objectives(P)')
    for i = 1:length(P)
        P[i].rank = R[i]
    end
    P
end

"""
    crowding_distance!((D, P)

Calculate crowding distance for individuals and save the results into `D`
given the fitness values `P` and collection of `F`.
"""
function crowding_distance!(D::AbstractVector, I::AbstractMatrix) where {T}
    N, _ = size(I)

    for m = eachcol(I)
        s = sortperm(m)
        D[s[1]] = Inf
        D[s[end]] = Inf

        max, _ = findmax(m)
        min, _ = findmin(m)

        for i = 2:N-1
            D[i] += (m[s[i+1]] - m[s[i-1]]) / (max - min)
        end
    end
end
function crowding_distance!(P::Population)
    maxF, _ = findmax(map(rank, P))
    F = [Vector{Int64}() for _ in 1:maxF]
    for i = 1:length(P)
        push!(F[rank(P[i])], i)
    end
    for f = F
        C = fill(0.0, length(f))
        crowding_distance!(C, objectives(P[f]))
        for i = 1:length(f)
            P[f][i].distance = C[i]
        end
    end
    P
end

"""
    compare(a, b)
    compares whether two vectors are dominated or not.
    Output:
    `1` if argument 1 (a) dominates argument 2 (b).
    `2` if argument 2 (b) dominates argument 1 (a).
    `3` if both arguments 1 (a) and 2 (b) are incomparable.
    `0` if both arguments 1 (a) and 2 (b) are equal.
"""
function compare(a::Vector, b::Vector)
    k = length(a)
    @assert k == length(b)

    i = 1
    while i <= k && a[i] == b[i]
        i += 1;
    end

    if i > k
        return 0 # equals
    end

    if a[i] < b[i]

        for j = i+1:k# (j = i+1; j <= k; ++j)
            if b[j] < a[j]
                return 3 #a and b are incomparable
            end
        end

        return 1; #  a dominates b
    end

    for j = i+1:k #(j = i+1; j < k; ++j) 
        if (a[j] < b[j])
            return 3 #; // a and b are incomparable
        end
    end

    return 2 # b dominates a
    
end


function compare(a::AbstractIndividual, b::AbstractIndividual)
    compare(objectives(a), objectives(b))
end

"""
    get_non_dominated_solutions_perm(population)
Return a vector of integers `v` such that `population[v]` are the non dominated
solutions contained in `population`.
"""
function get_non_dominated_solutions_perm(population)
    ids = Int[1]
    n = length(population)

    for i in 2:n
        j = 1
        while j <= length(ids)
            jj = ids[j]
            relation = compare(population[i], population[jj])
 
            if relation == 2 # j dominates i
                break
            elseif relation == 1 # i dominates j
                deleteat!(ids, j)
                continue
            end

            j += 1
        end

        if j > length(ids)
            push!(ids, i)
        end
        
    end

    return ids
end


"""
    get_non_dominated_solutions(population)
Return the non dominated solutions contained in `population`.
"""
function get_non_dominated_solutions(population)
    return population[get_non_dominated_solutions_perm(population)]
end