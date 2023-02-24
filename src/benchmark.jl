Base.@kwdef struct BenchmarkResult
    problem::Any
    algorithm::Any
    runs::Any
    igds::Any
    hvs::Any
end

function benchmark(problem, algorithm; runs = 30)
    igds = []
    hvs = []
    progress = Progress(
        runs * problem.maxFE,
        barglyphs = BarGlyphs("[=> ]"),
        barlen = 50,
        color = :yellow,
    )
    queue = ConcurrentQueue()
    tasks = []

    for i in 1:runs
        push!(queue, i)
    end

    l = Threads.ReentrantLock()
    for i = 1:Threads.nthreads()
        task = Threads.@spawn begin
            while (run = maybepopfirst!(queue)) !== nothing
                fcalls = 0
                state_callback = function (state)
                    next!(progress, step=state.fcalls - fcalls)
                    fcalls = state.fcalls
                end
                result =
                    optimize(problem, algorithm, options = Options(state_callback = state_callback))
                state = result.state
                igd = IGD(pfront(state), problem.truepf)
                push!(igds, igd)
                hv = HV(pfront(state), problem.truepf)
                push!(hvs, hv)
                lock(l)
                try
                    save_run(problem, algorithm, igd, hv)
                finally
                    unlock(l)
                end
            end
        end
        push!(tasks, task)
    end

    map(fetch, tasks)
    finish!(progress)

    @assert length(igds) == runs
    @assert length(hvs) == runs

    igd_min = minimum(igds)
    igd_max = maximum(igds)
    igd_mean = mean(igds)
    igd_std = std(igds)
    hv_min = minimum(hvs)
    hv_max = maximum(hvs)
    hv_mean = mean(hvs)
    hv_std = std(hvs)

    fmt7 = n -> @sprintf("%.7f", n)
    fmt4 = n -> @sprintf("%.4f", n)

    table_data = [[1:runs; "Means (Std.)"] [map(fmt7, igds); "$(fmt7(igd_mean)) ($(fmt7(igd_std)))"] [
        map(fmt7, hvs)
        "$(fmt7(hv_mean)) ($(fmt7(hv_std)))"
    ]]
    table_header = [typeof(problem), "IGD", "HV"]
    h1 = Highlighter((data, i, j) -> (i == size(data)[1]), bold = true)
    h2 = Highlighter(
        (data, i, j) -> (j == 2 && data[i, j] == igd_max),
        bold = true,
        foreground = :red,
    )
    h3 = Highlighter(
        (data, i, j) -> (j == 2 && data[i, j] == igd_min),
        bold = true,
        foreground = :green,
    )
    h4 = Highlighter(
        (data, i, j) -> (j == 3 && data[i, j] == hv_max),
        bold = true,
        foreground = :green,
    )
    h5 = Highlighter(
        (data, i, j) -> (j == 3 && data[i, j] == hv_min),
        bold = true,
        foreground = :red,
    )
    pretty_table(
        table_data,
        header = table_header,
        highlighters = (h1, h2, h3, h4, h5),
    )
    # BenchmarkResult(problem, algorithm, runs, igds, hvs)
end
