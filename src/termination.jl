abstract type ConvergenceMetric end

description(m::ConvergenceMetric) =
    error("`description` is not implemented for $(typeof(m)).")

"""
    converged(metric)

判断是否已收敛
"""
converged(m::ConvergenceMetric) = diff(m) <= tolerance(m)

diff(m::ConvergenceMetric) = m.Δ #error("`diff` is not implemented for $(typeof(m)).")

tolerance(m::ConvergenceMetric) = m.tol #error("`tolerance` is not implemented for $(typeof(m)).")

"""
    assess!(metric, state)

根据状态计算并判断是否已收敛
"""
assess!(m::ConvergenceMetric, s::AbstractOptimizerState) =
    error("`assess!` is not implemented for $(typeof(m)).")


const ConvergenceMetrics = Vector{ConvergenceMetric}

copy(cm::ConvergenceMetrics) = map(copy, cm)


# Utilities

abschange(curr, prev) = Float64(abs(curr - prev))
relchange(curr, prev) = Float64(abs(curr - prev) / abs(curr))

maxdiff(x::AbstractArray, y::AbstractArray) = mapreduce((a, b) -> abs(a - b), max, x, y)
abschange(curr::T, prev) where {T<:AbstractArray} = maxdiff(curr, prev)
relchange(curr::T, prev) where {T<:AbstractArray} = maxdiff(curr, prev) / maximum(abs, curr)

# Single value

"""
Absolute difference convergence for single objective optimization.

This convergence metric allows to estimate an absolute difference between consecutive
states of the optimization algorithm, and triggers convergence when,

- `|f(x) - f(x')| < ε`

where `ε` is a tolerance value, `x` and `x'` are previous and current minimizers
found by the optimization algorithm.
"""
mutable struct AbsDiff{T} <: ConvergenceMetric
    tol::T
    Δ::Float64
    value::T
end
AbsDiff(tol::T) where {T<:AbstractFloat} = AbsDiff(tol, Inf, zero(T))
AbsDiff() = AbsDiff(1e-12)
description(m::AbsDiff) = "|f(x) - f(x')|"
function assess!(m::AbsDiff, state::AbstractOptimizerState)
    val = value(state)
    m.Δ = abschange(val, m.value)
    m.value = val
    converged(m)
end
copy(m::AbsDiff) = AbsDiff(copy(m.tol), copy(m.Δ), copy(m.value))

"""
Relative difference convergence metric for single objective optimization.

This convergence metric allows to estimate a relative difference between consecutive
states of the optimization algorithm, and triggers convergence when,

- `|f(x) - f(x')|/|f(x')| < ε`

where `ε` is a tolerance value, `x` and `x'` are previous and current minimizers
found by the optimization algorithm.
"""
mutable struct RelDiff{T} <: ConvergenceMetric
    tol::T
    Δ::Float64
    value::T
end
RelDiff(tol::T) where {T<:AbstractFloat} = RelDiff(tol, Inf, zero(T))
RelDiff() = RelDiff(1e-12)
description(m::RelDiff) = "|f(x) - f(x')|/|f(x')|"
function assess!(m::RelDiff, state::AbstractOptimizerState)
    val = value(state)
    m.Δ = relchange(val, m.value)
    m.value = val
    converged(m)
end
copy(m::RelDiff) = RelDiff(copy(m.tol), copy(m.Δ), copy(m.value))



##########################
# Convergence Assessment #
##########################

function assess_convergence!(state)
    for metric in state.metrics
        converged = assess!(metric, state)
        if converged
            return true
        end
    end
    return false
end
