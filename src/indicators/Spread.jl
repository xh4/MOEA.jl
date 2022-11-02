"""
    spread(S,R)

Returns a diversity metric of a population of set `S` to the reference set `R`.
"""
function spread(S::AbstractMatrix, R::AbstractMatrix)
    n = size(S, 2)
    m = size(R, 2)
    Δₖ = [
        minimum(
            norm(view(S, :, i) - view(R, :, j)) for
            i in 1:n if view(S, :, i) != view(R, :, j)
        ) for j = 1:m
    ]
    Δ = mean(Δₖ)
    sum(abs.(Δₖ .- Δ)) / m * Δ
end

function spread(S::AbstractMatrix)
    n = size(S, 2)
    n == 1 && return NaN
    Δₖ = [
        minimum(
            norm(view(S, :, i) - view(S, :, j)) for
            j in 1:n if view(S, :, i) != view(S, :, j)
        ) for i = 1:n
    ]
    Δ = mean(Δₖ)
    sum(abs.(Δₖ .- Δ)) / n * Δ
end