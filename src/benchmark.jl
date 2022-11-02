Base.@kwdef struct BenchmarkResult
    problem
    algorithm
    runs
    igds
    hvs
end

function benchmark(problem, algorithm; runs = 30)
    igds = []
    hvs = []
    times = []
    progress = Progress(runs, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow)
    Threads.@threads for i in 1:runs
        time_start = time()
        result = optimize(problem, algorithm)
        time_stop = time()
        duration = time_stop - time_start
        push!(times, duration)
        state = result.trace[end]
        igd = IGD(pfront(state), problem.truepf)
        push!(igds, igd)
        hv = HV(pfront(state), problem.truepf)
        push!(hvs, hv)
        save_run(problem, algorithm, igd, hv)
        next!(progress)
    end
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
    time_min = minimum(times)
    time_max = maximum(times)
    time_mean = mean(times)

    table_data = [[1:runs; "Means (Std.)"] [igds; "$(igd_mean) ($(igd_std))"] [hvs; "$(hv_mean) ($(hv_std))"] [times; time_mean]]
    table_header = [typeof(problem), "IGD", "HV", "Time"]
    h1 = Highlighter((data,i,j) -> (i == size(data)[1]),
                        bold = true)
    h2 = Highlighter((data,i,j) -> (j == 2 && data[i,j] == igd_max),
                        bold = true,
                        foreground = :red)
    h3 = Highlighter((data,i,j) -> (j == 2 && data[i,j] == igd_min),
                        bold = true,
                        foreground = :green)
    h4 = Highlighter((data,i,j) -> (j == 3 && data[i,j] == hv_max),
                        bold = true,
                        foreground = :green)
    h5 = Highlighter((data,i,j) -> (j == 3 && data[i,j] == hv_min),
                        bold = true,
                        foreground = :red)
    h6 = Highlighter((data,i,j) -> (j == 4 && data[i,j] == time_min),
                        bold = true,
                        foreground = :green)
    h7 = Highlighter((data,i,j) -> (j == 4 && data[i,j] == time_max),
                        bold = true,
                        foreground = :red)
    pretty_table(table_data, header = table_header, highlighters = (h1, h2, h3, h4, h5, h6, h7))
    # BenchmarkResult(problem, algorithm, runs, igds, hvs)
end
