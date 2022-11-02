
struct OptimizationTraceRecord
    iteration::Int
    population
    value
    minimizer
    metadata::Dict{String,Any}
end

value(tr::OptimizationTraceRecord) = tr.value

function show(io::IO, t::OptimizationTraceRecord)
    print(io, lpad("$(t.iteration)", 6))
    print(io, "   ")
    print(io, lpad("$(t.value)", 14))
    for (key, value) in t.metadata
        print(io, "\n * $key: $value")
    end
    return
end

const OptimizationTrace = Vector{OptimizationTraceRecord}

function show(io::IO, tr::OptimizationTrace)
    print(io, "Iter     Function value\n")
    print(io, "------   --------------\n")
    for rec in tr
        show(io, rec)
        print(io, "\n")
    end
    return
end


abstract type AbstractOptimizationResult end


summary(or::AbstractOptimizationResult) = summary(or.method)

minimizer(r::AbstractOptimizationResult) = r.minimizer

minimum(r::AbstractOptimizationResult) = r.minimum

iterations(r::AbstractOptimizationResult) = r.iterations

iteration_limit_reached(r::AbstractOptimizationResult) = r.iteration_converged

trace(r::AbstractOptimizationResult) =
    length(r.trace) > 0 ? r.trace :
    error(
        "No trace in optimization results. To get a trace, run optimize() with store_trace = true.",
    )

f_calls(r::AbstractOptimizationResult) = r.f_calls

abstol(r::AbstractOptimizationResult) =
    error("`abstol` is not implemented for $(summary(r)).")

reltol(r::AbstractOptimizationResult) =
    error("`reltol` is not implemented for $(summary(r)).")

abschange(r::AbstractOptimizationResult) =
    error("`abschange` is not implemented for $(summary(r)).")
relchange(r::AbstractOptimizationResult) =
    error("`relchange` is not implemented for $(summary(r)).")


mutable struct OptimizationResult <: AbstractOptimizationResult
    problem
    algorithm
    state
    iterations
    trace
    time::Float64
end

function show(io::IO, r::OptimizationResult)
    print(io, " * Work counters\n")
    tr = round(r.time; digits = 4)
    print(io, "    Seconds run:   $tr\n")
    print(io, "    Iterations:    $(r.state.iteration)\n")
    print(io, "    f(x) calls:    $(r.state.fcalls)\n")
    return
end
