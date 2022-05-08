#########
# TRACE #
#########

function update!(
    tr::OptimizationTrace,
    state::S,
    iteration::Integer,
    population,
    val,
    minimizer,
    dt::Dict,
    store_trace::Bool,
    show_trace::Bool,
    show_every::Int = 1,
    callback = nothing,
) where {S<:AbstractOptimizerState,T,O}
    os = OptimizationTraceRecord(
        iteration,
        copy(population),
        val,
        minimizer,
        dt,
    )
    store_trace && push!(tr, os)
    if show_trace
        if iteration % show_every == 0
            show(os)
            print("\n")
            flush(stdout)
        end
    end
    if callback !== nothing && (iteration % show_every == 0)
        if store_trace
            stopped = callback(tr)
        else
            stopped = callback(os)
        end
    else
        stopped = false
    end
    stopped
end

function trace!(
    tr,
    iteration,
    objfun,
    state,
    method,
    options,
    curr_time = time(),
)
    dt = Dict{String,Any}()
    dt["time"] = curr_time
    # set additional trace value
    trace!(dt, objfun, state, method, options)
    update!(
        tr,
        state,
        iteration,
        state.fittest,
        value(state),
        minimizer(state),
        dt,
        options.store_trace,
        options.show_trace,
        options.show_every,
        options.callback,
    )
end

"""
    trace!(record::Dict{String,Any}, objfun, state, population, method, options)

Update the trace `record`. This function allows to supplement an additional information into the optimization algorithm trace by modifying a trace `record`. It can be overridden by specifying particular parameter types.
"""
trace!(record::Dict{String,Any}, objfun, state, method, options) = ()

########
# MISC #
########

default_values(x::AbstractArray{T}) where {T} = fill!(similar(x), zero(T))
default_values(x::AbstractArray{T}) where {T<:AbstractFloat} = fill!(similar(x), T(NaN))

function funargnum(f)
    fobj = try
        first(methods(f))
    catch
        @error "Cannot find function $f"
    end
    fobj.nargs
end

function vswap!(v1::T, v2::T, idx::Int) where {T<:AbstractVector}
    val = v1[idx]
    v1[idx] = v2[idx]
    v2[idx] = val
end

function dump_result(result::OptimizationResult, directory)
    mkpath(directory)
    result_file = joinpath(directory, "result.toml")
    save(joinpath(directory, "trace.jld"), "data", result.trace)
    if length(result.metrics) > 0
        save(joinpath(directory, "convergence.jld"), "data", result.metrics)
    end
    open(result_file, "w") do io
        data = Dict(
            "algorithm" => string(result.method),
            "iterations" => result.iterations,
            "converged" => result.converged,
            "trace" => "trace.jld",
            "time_run" => time_run(result),
            "time_limit" => time_limit(result),
            "f_calls" => f_calls(result),
            "minimum" => minimum(result),
            "minimizer" => minimizer(result),
        )
        if length(result.metrics) > 0
            data["convergence"] = "convergence.jld"
        end
        TOML.print(io, data)
    end
end

function draw_convergence(result::OptimizationResult, directory)
    images_directory = joinpath(directory, "convergence")
    mkpath(images_directory)

end

function draw_generations(result::OptimizationResult, directory=pwd())
    images_directory = joinpath(directory, "tmp", "generations")
    mkpath(images_directory)
    gen = 0
    anim = @animate for state in result.trace
        gen += 1
        println("Plotting generation $(gen) / $(length(result.trace))")
        for i in state.population
            if length(variables(i)) > 3 || length(variables(i)) < 2
                error("Can't draw individual with dimensions of $(length(variables(i)))")
            end
        end
        fig = Plots.plot(
            [variables(i)[1] for i in state.population],
            [variables(i)[2] for i in state.population],
            seriestype = :scatter,
            title = "Generation $(gen)",
        )
        xlims!((-3, 3))
        ylims!((-3, 3))
        # contour!(-3:0.1:3, -3:0.1:3, (x, y) -> x^2+y^2, legend = false, colorbar = false)
        savefig(fig, joinpath(images_directory, "$(gen).png"))
    end
    gif(anim, joinpath(directory, "generations.gif"), fps = 15, loop = 1)
end

function draw_pareto_fronts(result::OptimizationResult, directory=pwd())
    images_directory = joinpath(directory, "tmp", "objectives")
    mkpath(images_directory)
    gen = 0
    anim = @animate for state in result.trace[2:end]
        gen += 1
        # println("Plotting generation $(gen) / $(length(result.trace))")
        mF = maxF(Population(state.population))
        F = calcF(Population(state.population))
        cols = distinguishable_colors(mF, [RGB(1,1,1), RGB(0,0,0)], dropseed=false)
        fig = nothing
        for i = 1:length(F)
            f = F[i]
            if length(f) == 0
                continue
            end
            if i == 1
                fig = Plots.plot(
                    objectives(Population(state.population[f]))[:, 1],
                    objectives(Population(state.population[f]))[:, 2],
                    seriestype = :scatter,
                    color = cols[i],
                    title = "Generation $(gen)",
                    legend = false,
                )
            else
                Plots.plot!(
                    objectives(Population(state.population[f]))[:, 1],
                    objectives(Population(state.population[f]))[:, 2],
                    seriestype = :scatter,
                    color = cols[i]
                )
            end
        end
        # xlims!((0, 5))
        # ylims!((0, 5))
        savefig(fig, joinpath(images_directory, "$(gen).png"))
    end
    gif(anim, joinpath(images_directory, "objectives.gif"), fps = 20, loop = 1)
end

function draw_objectives(result::OptimizationResult, directory=pwd())
    images_directory = joinpath(directory, "tmp", "objectives")
    mkpath(images_directory)
    gen = 0
    anim = @animate for state in result.trace[2:end]
        gen += 1
        # println("Plotting generation $(gen) / $(length(result.trace))")
        fig = Plots.plot(
            objectives(Population(pfront(state)))[:, 1],
            objectives(Population(pfront(state)))[:, 2],
            seriestype = :scatter,
            title = "Generation $(gen)",
            legend = false,
        )
        savefig(fig, joinpath(images_directory, "$(gen).png"))
    end
    gif(anim, joinpath(images_directory, "objectives.gif"), fps = 20, loop = 1)
end

function draw_metrics(result::OptimizationResult)
    plot(
        1:length(result.trace),
        [state.metrics[2].value for state in result.trace],
        legend = false,
    )
end

function draw_igd(result, truepf)
    Plots.plot(
        1:length(result.trace),
        [state_igd(state, truepf) for state in result.trace],
        legend = false,
    )
end

function state_igd(state, truepf)
    igd(pfront(state), truepf)
end
