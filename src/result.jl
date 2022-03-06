
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
    method
    minimizer
    minimum
    iterations::Int
    iteration_converged::Bool
    converged::Bool
    metrics::ConvergenceMetrics
    f_calls::Int
    trace
    time_limit::Float64
    time_run::Float64
    is_moo::Bool
end

converged(r::OptimizationResult) = r.converged
time_limit(r::OptimizationResult) = r.time_limit
time_run(r::OptimizationResult) = r.time_run
is_moo(r::OptimizationResult) = r.is_moo

function show(io::IO, r::OptimizationResult)
    failure_string = "failure"
    if iteration_limit_reached(r)
        failure_string *= " (reached maximum number of iterations)"
    end
    if time_run(r) > time_limit(r)
        failure_string *= " (exceeded time limit of $(time_limit(r)))"
    end
    print(io, "\n")
    print(io, " * Status: ", converged(r) ? "success" : failure_string, "\n\n")
    print(io, " * Candidate solution\n")
    mzr = minimizer(r)
    if is_moo(r)
        pfsize = size(mzr)[1]
        pl = pfsize > 1 ? "s" : ""
        print(io, "    Pareto front: $(pfsize) element$pl\n")
    elseif mzr isa AbstractVector
        nx = length(mzr)
        str_x_elements = ["$_x" for _x in Iterators.take(mzr, min(nx, 3))]
        if nx >= 4
            push!(str_x_elements, " ...")
        end
        print(io, "    Minimizer:  [", join(str_x_elements, ", "), "]\n")
    else
        print(io, "    Minimizer:  $mzr\n")
    end
    !is_moo(r) && print(io, "    Minimum:    $(minimum(r))\n")
    print(io, "    Iterations:", rpad(" ", is_moo(r) ? 3 : 0), "$(iterations(r))\n")
    print(io, "\n")
    print(io, " * Found with\n")
    print(io, "    Algorithm: $(summary(r))\n")
    print(io, "\n")
    if length(r.metrics) > 0
        print(io, " * Convergence measures\n")
        maxdsclen = maximum(length(description(cm)) for cm in r.metrics)
        rpd = " "^4
        for cm in r.metrics
            sgn = converged(cm) ? "≤" : "≰"
            dsc = description(cm)
            lpd = " "^(maxdsclen + 1 - length(dsc))
            print(io, "$rpd$dsc$lpd= $(diff(cm)) $sgn $(tolerance(cm))\n")
        end
        print(io, "\n")
    end
    print(io, " * Work counters\n")
    tr = round(time_run(r); digits = 4)
    tl = isnan(time_limit(r)) ? Inf : round(time_limit(r); digits = 4)
    print(io, "    Seconds run:   $tr (vs limit $tl)\n")
    print(io, "    Iterations:    $(iterations(r))\n")
    print(io, "    f(x) calls:    $(f_calls(r))\n")
    return
end
