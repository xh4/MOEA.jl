"""
	hypervolume(front, reference_point)
Computes the hypervolume indicator, i.e., volume between points in `front` and `reference_point`.
Note that each point in `front` must (weakly) dominates to `reference_point`. Also, `front`
is a non-dominated set.
If `front::State` and `reference_point::Vector`, then computes `hypervolume(front.population, reference_point)` after
ignoring solutions in `front` that do not dominate `reference_point`.
"""
function hypervolume(front::Array{Vector{T}}, reference_point::Vector) where {T<:Real}
    weaklyDominates(point, other) = begin
        for i in 1:length(point)
            if point[i] > other[i]
                return false
            end
        end
        return true
    end

    relevantPoints = Vector[]
    for point in front
        # only consider points that dominate the reference point
        if weaklyDominates(point, reference_point)
            push!(relevantPoints, point)
        end
    end

    if length(relevantPoints) != length(front)
        ign = length(front) - length(relevantPoints)
        rel = length(relevantPoints)
        # @warn "Ignoring $ign points dominated by the reference point ($rel points are used)."
    end

    return FPL(relevantPoints, reference_point)
end

function hypervolume(front::AbstractMatrix, reference_point::Vector)
    front_ = [ front[i,:] for i in 1:size(front,1) ]
    hypervolume(front_, reference_point)
end

function HV(pfront::T, truepf::T) where {T<:AbstractMatrix}
    N, M = size(pfront)
    fmin = minimum(vcat(minimum(pfront, dims=1), zeros(1,M)), dims=1)
    fmax = maximum(truepf, dims=1)
    pfront = (pfront - repeat(fmin, N, 1)) ./ repeat((fmax-fmin)*1.1, N, 1)
    ref_point = fill(1.0, M)
    hypervolume(pfront, ref_point)
end

function HV(pfront::Population, truepf::Population)
    HV(objectives(pfront), objectives(truepf))
end