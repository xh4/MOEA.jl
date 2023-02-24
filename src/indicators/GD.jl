"""
    GD(A,R)

Calculate a generational distance between set `A` and the reference set `R`.
This metric measures the convergence, i.e. closeness of the non-dominated solutions
to the Pareto front, of a population.
"""
function GD(A::AbstractMatrix, R::AbstractMatrix)
    na, da = size(A)
    nr, dr = size(R)
    (na == 0 || nr == 0) && return Inf
    sum = 0
    for a in eachrow(A)
        sum += minimum(norm(a - r) for r in eachrow(R))
    end
    sum / na
end
GD(A::Population, R::Population) = GD(objectives(A), objectives(R))