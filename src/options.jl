"""
There are following options available:
- `abstol::Float64`: the absolute tolerance used in the convergence test (*deprecated, use `metrics` parameter of a particular optimization algorithm*)
- `reltol::Float64`: the relative tolerance used in the convergence test (*deprecated, use `metrics` parameter of a particular optimization algorithm*)
- `successive_f_tol::Integer`: the additional number of the iterations of the optimization algorithm after the convergence test is satisfied (*default: 10*)
- `iterations::Integer`: the total number of the iterations of the optimization algorithm (*default: 1000*)
- `show_trace::Bool`: enable the trace information display during the optimization (*default: false*).
- `store_trace::Bool`: enable the trace information capturing during the optimization (*default: false*). The trace can be accessed by using [`trace`](@ref) function after optimization is finished.
- `show_every::Integer`: show every `n`s successive trace message (*default: 1*)
- `time_limit::Float64`: the time limit for the optimization run in seconds. If the value set to `NaN` then the limit is not set. (*default: NaN*)
- `callback`: the callback function that is called after each iteration of the optimization algorithm. The function accepts as parameter a trace dictionary, and **must** return a `Bool` value which if `true` terminates the optimization. (*default: nothing*)
- `parallelization::Symbol`: allows parallelization of the population fitness evaluation if set to `:thread` using multiple threads (*default: `:serial`*)
- `rng::AbstractRNG`: a random number generator object that is used to control generation of random data during the evolutionary optimization (*default: `Random.default_rng()`*)
"""
@kwdef struct Options{TCallback<:Union{Nothing,Function},TRNG<:AbstractRNG}
    abstol::Float64 = Inf
    reltol::Float64 = Inf
    successive_f_tol::Int = 10
    iterations::Int = 1000
    store_trace::Bool = false
    show_trace::Bool = false
    show_every::Int = 1
    callback::TCallback = nothing
    time_limit::Float64 = NaN
    parallelization::Symbol = :serial
    rng::TRNG = Random.default_rng()
end

function show(io::IO, o::Options)
    for k in fieldnames(typeof(o))
        v = getfield(o, k)
        if v === nothing
            print(io, lpad("$(k)", 24) * " = nothing\n")
        else
            print(io, lpad("$(k)", 24) * " = $v\n")
        end
    end
end
