abstract type Experiment end

@tags table tr th td thead tbody tfoot head meta body h1 h2

function experiment(problems, algorithms; runs = 10, indicators = [:IGD, :HV])
    queue = ConcurrentQueue()
    tasks = []
    result = fill(NaN, length(indicators), length(problems), length(algorithms), runs)
    progress_counter = reduce((n, p) -> n + runs * p.maxFE * length(algorithms), problems, init=0)

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
            progress_counter -= problem.maxFE
        end
    end

    progress = Progress(
        progress_counter,
        barglyphs = BarGlyphs("[=> ]"),
        barlen = 50,
        color = :yellow,
    )

    # 生成任务，放入队列
    for v in shuffle([
        v for v in Iterators.product(1:length(problems), 1:length(algorithms), 1:runs)
    ])
        problem_index = v[1]
        algorithm_index = v[2]
        run = v[3]
        if isnan(result[1, problem_index, algorithm_index, run])
            push!(queue, v)
        end
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
                    igd = IGD(pfront(state), problem)
                    hv = HV(pfront(state), problem)

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
                            show_progress = false
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
        fmt7 = n -> @sprintf("%.7f", n)
        for (problem_index, problem) in enumerate(problems)
            for (algorithm_index, algorithm) in enumerate(algorithms)
                is_main_algorithm = algorithm_index == main_algorithm_index
                this_values = [
                    getindex(result, indicator_index, problem_index, algorithm_index, run) for run = 1:runs
                ]
                this_mean = mean(this_values)
                this_std = std(this_values)
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
                            ranksum_data[2, algorithm_index] += 1
                        else
                            ranksum_symbol = "-"
                            ranksum_data[1, algorithm_index] += 1
                        end
                    else # 无明显差异
                        ranksum_symbol = "≈"
                        ranksum_data[3, algorithm_index] += 1
                    end
                end
                table_data[
                    problem_index,
                    algorithm_index,
                ] = "$(fmt7(this_mean)) ($(fmt7(this_std))) $(ranksum_symbol)"
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

        # HTML report
        # tables = []
        # for (indicator_index, indicator) in enumerate(indicators)
        #     h = thead(tr([th(""); map(alg -> th(identifier(alg)), algorithms)]))
        #     b = tbody(tr([td("1"), td("2")]))
        #     f = tfoot(tr([th("+/-/~"), th("-")]))
        #     t = table(h, b, f)
        #     tables = [tables; [h2(indicator), t]]
        # end
        # doc = [
        #     head(
        #     meta(charset="UTF-8"),
        #     ),
        #     body(
        #         [
        #         h1("Experiment"); tables
        #         ] )
        # ]
        # savehtml("c:/users/xh/report.html", doc)
        # DefaultApplication.open("c:/users/xh/report.html")

        # LaTeX report
        file = open("result.tex", "w")
        for (indicator_index, indicator) in enumerate(indicators)
            if indicator_index > 1
                println(file, "")
                println(file, "")
            end
            println(file, "\\begin{table}[!h]")
            println(file, "\\centering")
            print(file, "\\begin{tabular}{l")
            for i = 1:length(algorithms)
                print(file, " l")
            end
            println(file, "}")
            println(file, "\\hline")
            print(file, "Problem")
            for (algorithm_index, algorithm) in enumerate(algorithms)
                print(file, " & $(name(algorithm))")
            end
            println(file, "\\\\")
            println(file, "\\hline")
            ranksum_data = fill(0, 3, length(algorithms))
            fmt6 = n -> @sprintf("%.6f", n)
            for (problem_index, problem) in enumerate(problems)
                print(file, "$(typeof(problem))")
                for (algorithm_index, algorithm) in enumerate(algorithms)
                    is_main_algorithm = algorithm_index == main_algorithm_index
                    this_values = [
                        getindex(result, indicator_index, problem_index, algorithm_index, run) for run = 1:runs
                    ]
                    this_mean = mean(this_values)
                    this_std = std(this_values)
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
                                ranksum_data[2, algorithm_index] += 1
                            else
                                ranksum_symbol = "-"
                                ranksum_data[1, algorithm_index] += 1
                            end
                        else # 无明显差异
                            ranksum_symbol = "\\approx"
                            ranksum_data[3, algorithm_index] += 1
                        end
                    end
                    
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
                    print(file, " & ")
                    if is_best
                        print(file, "\\textbf{")
                    end
                    print(file, "$(fmt6(this_mean)) ($(fmt6(this_std)))\$^{$(ranksum_symbol)}\$")
                    if is_best
                        print(file, "}")
                    end
                end
                println(file, "\\\\")
            end
            print(file, "\$+\$/\$-\$/\$\\approx\$")
            for (algorithm_index, algorithm) in enumerate(algorithms)
                print(file, " & ")
                if !(algorithm_index == length(algorithms))
                    print(file, "$(getindex(ranksum_data, 1, algorithm_index))/$(getindex(ranksum_data, 2, algorithm_index))/$(getindex(ranksum_data, 3, algorithm_index))")
                end
            end
            println(file, "\\\\")
            println(file, "\\hline")
            println(file, "\\end{tabular}")
            println(file, "\\caption{\\label{tab:$(typeof(indicator))}$(typeof(indicator))}")
            println(file, "\\end{table}")
        end
        close(file)

        # LaTeX report 2
        file = open("result-2.tex", "w")
        for (indicator_index, indicator) in enumerate(indicators)
            if indicator_index > 1
                println(file, "")
                println(file, "")
            end
            println(file, "\\begin{table}[!h]")
            println(file, "\\centering")
            print(file, "\\begin{tabular}{c")
            for i = 1:length(algorithms)
                print(file, " c")
            end
            println(file, "}")
            println(file, "\\hline")
            print(file, "Problem")
            for (algorithm_index, algorithm) in enumerate(algorithms)
                print(file, " & $(name(algorithm))")
            end
            println(file, "\\\\")
            println(file, "\\hline")
            ranksum_data = fill(0, 3, length(algorithms))
            fmt6 = n -> @sprintf("%.6f", n)
            for (problem_index, problem) in enumerate(problems)
                print(file, "$(typeof(problem))")
                for (algorithm_index, algorithm) in enumerate(algorithms)
                    is_main_algorithm = algorithm_index == main_algorithm_index
                    this_values = [
                        getindex(result, indicator_index, problem_index, algorithm_index, run) for run = 1:runs
                    ]
                    this_mean = mean(this_values)
                    this_std = std(this_values)
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
                        p = pvalue(MannWhitneyUTest(this_values, main_values))
                        # 有明显差异
                        if p <= 0.05
                            if comparator(this_mean, main_mean)
                                ranksum_data[2, algorithm_index] += 1
                            else
                                ranksum_data[1, algorithm_index] += 1
                            end
                        else # 无明显差异
                            ranksum_data[3, algorithm_index] += 1
                        end
                    end
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
                    print(file, " & ")
                    if is_best
                        print(file, "\\textbf{")
                    end
                    print(file, "$(@sprintf("%.3e", this_mean)) ($(@sprintf("%.3e", this_std)))")
                    if is_best
                        print(file, "}")
                    end
                end
                println(file, "\\\\")
                # p
                print(file, "\$p\$")
                for (algorithm_index, algorithm) in enumerate(algorithms)
                    is_main_algorithm = algorithm_index == main_algorithm_index
                    this_values = [
                        getindex(result, indicator_index, problem_index, algorithm_index, run) for run = 1:runs
                    ]
                    this_mean = mean(this_values)
                    this_std = std(this_values)
                    if indicator == :IGD
                        comparator = <
                    elseif indicator == :HV
                        comparator = >
                    end
                    # 秩和检验
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
                        p = pvalue(MannWhitneyUTest(this_values, main_values))
                        # 有明显差异
                        print(file, " & $(@sprintf("%.3e", p))")
                    else
                        print(file, " & -")
                    end
                end
                println(file, "\\\\")
                # h
                print(file, "\$h\$")
                for (algorithm_index, algorithm) in enumerate(algorithms)
                    is_main_algorithm = algorithm_index == main_algorithm_index
                    this_values = [
                        getindex(result, indicator_index, problem_index, algorithm_index, run) for run = 1:runs
                    ]
                    this_mean = mean(this_values)
                    this_std = std(this_values)
                    if indicator == :IGD
                        comparator = <
                    elseif indicator == :HV
                        comparator = >
                    end
                    # 秩和检验
                    h = -1
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
                        p = pvalue(MannWhitneyUTest(this_values, main_values))
                        # 有明显差异
                        if p <= 0.05
                            h = 1
                        else # 无明显差异
                            h = 0
                        end
                        print(file, " & $(h)")
                    else
                        print(file, " & -")
                    end
                end
                println(file, "\\\\")
            end
            print(file, "\$+\$/\$-\$/\$\\approx\$")
            for (algorithm_index, algorithm) in enumerate(algorithms)
                print(file, " & ")
                if !(algorithm_index == length(algorithms))
                    print(file, "$(getindex(ranksum_data, 1, algorithm_index))/$(getindex(ranksum_data, 2, algorithm_index))/$(getindex(ranksum_data, 3, algorithm_index))")
                else
                    print(file, "-")
                end
            end
            println(file, "\\\\")
            println(file, "\\hline")
            println(file, "\\end{tabular}")
            println(file, "\\caption{\\label{tab:$(typeof(indicator))}$(typeof(indicator))}")
            println(file, "\\end{table}")
        end
        close(file)

        pretty_table(
            table_data,
            header = table_header,
            highlighters = Tuple(table_highlighters),
        )
    end
end

function table_report()

end

function html_report()

end
