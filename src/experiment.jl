abstract type Experiment end

function experiment(problems, algorithms; runs = 30)
    indicators = [:IGD, :HV]

    queue = ConcurrentQueue()
    tasks = []
    result = fill(NaN, length(indicators), length(problems), length(algorithms), runs)

    progress = Progress(
        reduce((n, p) -> n + runs * p.maxFE * length(algorithms), problems, init=0),
        barglyphs = BarGlyphs("[=> ]"),
        barlen = 50,
        color = :yellow,
    )

    # 加载数据库中保存的结果
    for (problem_index, algorithm_index) in
        Iterators.product(1:length(problems), 1:length(algorithms))
        problem = problems[problem_index]
        algorithm = algorithms[algorithm_index]
        records = fetch_runs(problem, algorithm, runs)
        for run = 1:size(records)[1]
            igd = records[run, 1]
            hv = records[run, 2]
            result[1, problem_index, algorithm_index, run] = igd
            result[2, problem_index, algorithm_index, run] = hv
        end
    end

    # 生成任务，放入队列
    for v in shuffle([
        v for v in Iterators.product(1:length(problems), 1:length(algorithms), 1:runs)
    ])
        push!(queue, v)
    end

    # 运行实验
    l = Threads.ReentrantLock()
    for i = 1:Threads.nthreads()
        task = Threads.@spawn begin
            while (v = maybepopfirst!(queue)) !== nothing
                problem_index, algorithm_index, run = something(v)
                problem = problems[problem_index]
                algorithm = algorithms[algorithm_index]

                if !isnan(result[1, problem_index, algorithm_index, run])
                    continue
                end

                fcalls = 0
                state_callback = function (state)
                    try
                        next!(progress, step=state.fcalls-fcalls)
                        fcalls = state.fcalls
                    catch e
                        @error "Error in state callback!" exception = e
                        Base.show_backtrace(stderr, catch_backtrace())
                    end
                end

                finish_callback = function (state)
                    igd = IGD(pfront(state), problem.truepf)
                    hv = HV(pfront(state), problem.truepf)

                    result[1, problem_index, algorithm_index, run] = igd
                    result[2, problem_index, algorithm_index, run] = hv

                    lock(l)
                    try
                        save_run(problem, algorithm, igd, hv)
                    finally
                        unlock(l)
                    end
                end

                population_size = algorithm.N
                try
                    # println("Optimize $(problem) with $(algorithm) run $(run)")
                    optimize(
                        problem,
                        algorithm,
                        options = Options(
                            state_callback = state_callback,
                            finish_callback = finish_callback,
                        ),
                    )
                catch e
                    @error "Error in state callback!" exception = e
                    Base.show_backtrace(stderr, catch_backtrace())
                end
            end
        end
        push!(tasks, task)
    end
    map(fetch, tasks)
    finish!(progress)

    # 处理结果
    for (indicator_index, indicator) in enumerate(indicators)
        table_header = [indicator; map(typeof, algorithms)]
        table_highlighters = []
        main_algorithm_index = length(algorithms)
        table_data = fill("", length(problems), length(algorithms))
        ranksum_data = fill(0, 3, length(algorithms))
        for (problem_index, problem) in enumerate(problems)
            for (algorithm_index, algorithm) in enumerate(algorithms)
                is_main_algorithm = algorithm_index == main_algorithm_index
                this_values = [
                    getindex(result, indicator_index, problem_index, algorithm_index, run) for run = 1:runs
                ]
                this_mean = mean(this_values)
                if indicator == :IGD
                    comparator = <
                elseif indicator == :HV
                    comparator = >
                end
                # 秩和检验
                ranksum_symbol = ""
                if !is_main_algorithm
                    main_values = [
                        getindex(
                            result,
                            indicator_index,
                            problem_index,
                            main_algorithm_index,
                            run,
                        ) for run = 1:runs
                    ]
                    main_mean = mean(main_values)
                    ranksum_symbol = "?"
                    p = pvalue(MannWhitneyUTest(this_values, main_values))
                    # 有明显差异
                    if p <= 0.05
                        if comparator(this_mean, main_mean)
                            ranksum_symbol = "+"
                            ranksum_data[1, algorithm_index] += 1
                        else
                            ranksum_symbol = "-"
                            ranksum_data[2, algorithm_index] += 1
                        end
                    else # 无明显差异
                        ranksum_symbol = "≈"
                        ranksum_data[3, algorithm_index] += 1
                    end
                end
                table_data[
                    problem_index,
                    algorithm_index,
                ] = "$(this_mean) $(ranksum_symbol)"
                # Compare with others
                other_algorithm_indexes =
                    filter(i -> i != algorithm_index, 1:length(algorithms))
                other_algorithm_means = []
                for other_algorithm_index in other_algorithm_indexes
                    other_values = [
                        getindex(
                            result,
                            indicator_index,
                            problem_index,
                            other_algorithm_index,
                            run,
                        ) for run = 1:runs
                    ]
                    other_mean = mean(other_values)
                    push!(other_algorithm_means, other_mean)
                end
                is_best = all(
                    map(
                        other_mean -> comparator(this_mean, other_mean),
                        other_algorithm_means,
                    ),
                )
                if is_best
                    h = Highlighter(
                        f = (data, i, j) ->
                            (i == problem_index && j - 1 == algorithm_index),
                        crayon = crayon"green bold",
                    )
                    push!(table_highlighters, h)
                end
            end
        end
        ranksums = map(
            algorithm_index ->
                "$(getindex(ranksum_data, 1, algorithm_index))/$(getindex(ranksum_data, 2, algorithm_index))/$(getindex(ranksum_data, 3, algorithm_index))",
            1:length(algorithms),
        )
        ranksums[main_algorithm_index] = ""
        footer = permutedims(["+/-/≈"; ranksums])
        table_data = [[map(typeof, problems) table_data]; footer]
        pretty_table(
            table_data,
            header = table_header,
            highlighters = Tuple(table_highlighters),
        )
    end
end

platemo_data_path = "C:/Users/XH/PlatEMO 3.4/Data/"

function find_platemo_mat_files(problem, algorithm)
    "$(typeof(algorithm).name.name)/$(typeof(algorithm).name.name)_$(typeof(problem).name.name).mat"
end
