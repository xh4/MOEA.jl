Base.@kwdef mutable struct ExperimentDesc
    algorithms = []
    problems = []
    indicators = []
end

function get_experiment()
    ExperimentDesc(
        algorithms = [get_algorithm_desc("NSGA-II"),
                      get_algorithm_desc("MOEA/D"),
                      # get_algorithm_desc("MOEA/D-DE"),
                      get_algorithm_desc("MOEA/D-MY"),
                      # get_algorithm_desc("MOEA/D-DRA"),
                      ],
        problems = [# get_problem_desc("ZDT1"),
                    # get_problem_desc("ZDT2"),
                    # get_problem_desc("ZDT3"),
                    # get_problem_desc("ZDT4"),
                    # get_problem_desc("ZDT6"),
                    get_problem_desc("DTLZ1"),
                    # get_problem_desc("DTLZ2"),
                    # get_problem_desc("DTLZ3"),
                    # get_problem_desc("DTLZ4"),
                    # get_problem_desc("DTLZ5"),
                    # get_problem_desc("DTLZ6")
                    ],
        indicators = [igd_indicator_desc(), hv_indicator_desc()]
    )
end

function experiment_window()

    experiment = get_experiment()

    algorithms = experiment.algorithms
    problems = experiment.problems
    indicators = experiment.indicators

    runs = 30
    queue = ConcurrentQueue()
    tasks = []
    result = fill(NaN, length(indicators), length(problems), length(algorithms), runs)

    function ()
        CImGui.Begin("Experiment")

        CImGui.Text("=============== Algorithms ===============")

        for algorithm in algorithms
            CImGui.Text(algorithm.name)

            for param in algorithm.params
                CImGui.Text("")
                CImGui.SameLine(100)
                CImGui.Text(param.name)
                CImGui.SameLine(200)
                CImGui.Text(string(param.value))
            end
        end

        CImGui.Text("")
        CImGui.Text("=============== Problems ===============")
        CImGui.Text("")

        for problem in problems
            CImGui.Text(problem.name)

            for param in problem.params
                CImGui.Text("")
                CImGui.SameLine(100)
                CImGui.Text(param.name)
                CImGui.SameLine(200)
                CImGui.Text(string(param.value))
            end
        end

        CImGui.Text("")
        CImGui.Text("=============== Results ===============")

        for (indicator_index, indicator) in enumerate(indicators)
            CImGui.Text("")
            if CImGui.BeginTable("algorithms", length(algorithms)+1)
                flags = CImGui.ImGuiTableColumnFlags_WidthFixed
                CImGui.TableSetupColumn(indicator.name, flags, Float32(80.0))
                for (algorithm_index, algorithm) in enumerate(algorithms)
                    CImGui.TableSetupColumn(algorithm.name)
                end
                CImGui.TableHeadersRow()
                CImGui.TableNextRow()
                CImGui.TableSetColumnIndex(0)
                if CImGui.BeginTable("problems", 1)
                    CImGui.TableNextRow()
                    CImGui.TableSetColumnIndex(0)
                    CImGui.Text("")
                    for (problem_index, problem) in enumerate(problems)
                        CImGui.TableNextRow()
                        CImGui.TableSetColumnIndex(0)
                        CImGui.Text(problem.name)
                    end
                    CImGui.EndTable()
                end
                for (algorithm_index, algorithm) in enumerate(algorithms)
                    CImGui.TableSetColumnIndex(algorithm_index)
                    flags = CImGui.ImGuiTableFlags_RowBg
                    if CImGui.BeginTable("values", 4, flags)
                        CImGui.TableSetupColumn("Mean")
                        CImGui.TableSetupColumn("Min")
                        CImGui.TableSetupColumn("Max")
                        CImGui.TableSetupColumn("Std.")
                        CImGui.TableHeadersRow()

                        for (problem_index, problem) in enumerate(problems)
                            values = [getindex(result, indicator_index, problem_index, algorithm_index, run) for run in 1:runs]
                            notrun = all(map(isnan, values))
                            filter!(v->!isnan(v), values)
                            CImGui.TableNextRow()
                            # Mean
                            CImGui.TableSetColumnIndex(0)
                            if notrun
                                CImGui.Text("-")
                            else
                                v = mean(values)

                                # Compare with others
                                other_algorithm_indexes = filter(i -> i!=algorithm_index, 1:length(algorithms))
                                other_algorithm_means = []
                                for other_algorithm_index in other_algorithm_indexes
                                    other_values = [getindex(result, indicator_index, problem_index, other_algorithm_index, run)
                                                    for run in 1:runs]
                                    filter!(v->!isnan(v), other_values)
                                    other_v = mean(other_values)
                                    push!(other_algorithm_means, other_v)
                                end
                                all_better = all(map(other_v -> indicator.better(v, other_v), other_algorithm_means))
                                if all_better
                                    bg_color = CImGui.GetColorU32(0, 255, 0, 0.65)
                                    CImGui.TableSetBgColor(CImGui.ImGuiTableBgTarget_CellBg, bg_color)
                                end

                                CImGui.Text(@sprintf("%.6f", v))
                            end

                            # Min
                            CImGui.TableSetColumnIndex(1)
                            if notrun
                                CImGui.Text("-")
                            else
                                v,_ = findmin(values)
                                CImGui.Text(@sprintf("%.6f", v))
                            end

                            # Max
                            CImGui.TableSetColumnIndex(2)
                            if notrun
                                CImGui.Text("-")
                            else
                                v, _ = findmax(values)
                                CImGui.Text(@sprintf("%.6f", v))
                            end

                            # Std.
                            CImGui.TableSetColumnIndex(3)
                            if notrun
                                CImGui.Text("-")
                            else
                                v = std(values)
                                CImGui.Text(@sprintf("%.6f", v))
                            end
                        end
                        CImGui.EndTable()
                    end
                end
                CImGui.EndTable()
            end
        end

        CImGui.Text("")
        CImGui.Text("=============== Runs ===============")

        if CImGui.BeginTable("algorithms", length(algorithms)+1)
            flags = CImGui.ImGuiTableColumnFlags_WidthFixed
            CImGui.TableSetupColumn("", flags, Float32(80.0))
            for (algorithm_index, algorithm) in enumerate(algorithms)
                CImGui.TableSetupColumn(algorithm.name)
            end
            CImGui.TableHeadersRow()
            for (problem_index, problem) in enumerate(problems)
                CImGui.TableNextRow()
                CImGui.TableSetColumnIndex(0)
                CImGui.Text(problem.name)

                for (algorithm_index, algorithm) in enumerate(algorithms)
                    CImGui.TableSetColumnIndex(algorithm_index)
                    values = [getindex(result, 1, problem_index, algorithm_index, run) for run in 1:runs]
                    filter!(v->!isnan(v), values)
                    CImGui.Text(@sprintf("%d / %d", length(values), runs))
                end
            end
            CImGui.EndTable()
        end

        if CImGui.Button("Run")
            queue = ConcurrentQueue()
            tasks = []
            result = fill(NaN, length(indicators), length(problems), length(algorithms), runs)

            for v in shuffle([v for v in Iterators.product(1:length(problems), 1:length(algorithms), 1:runs)])
                push!(queue, v)
            end

            for i = 1:Threads.nthreads()
                task = Threads.@spawn begin
                    while (v = maybepopfirst!(queue)) !== nothing
                        problem_index, algorithm_index, run = something(v)
                        problem = problems[problem_index]
                        algorithm = algorithms[algorithm_index]
                        constraints = BoxConstraints(problem.bounds)

                        state_callback = function (state)
                            try
                                generation = state.iteration
                                
                            catch e
                                @error "Error in state callback!" exception = e
                                Base.show_backtrace(stderr, catch_backtrace())
                            end
                        end

                        finish_callback = function(state)

                            igd = IGD(pfront(state), problem.truepf)
                            hv = HV(pfront(state), problem.truepf)

                            @info "" run igd hv

                            result[1, problem_index, algorithm_index, run] = igd
                            result[2, problem_index, algorithm_index, run] = hv
                        end

                        # TODO: initial population
                        population_size = get_algorithm_param(algorithm, "N")
                        try
                            _, m = size(variables(problem.truepf))
                            optimize(problem.fn,
                                    algorithm.method,
                                    constraints = constraints,
                                    population = [Individual(rand(m)) for i in 1:population_size],
                                    options = Options(maxFE = 10000, state_callback = state_callback, 
                                                    finish_callback = finish_callback))
                        catch e
                            @error "Error in state callback!" exception = e
                            Base.show_backtrace(stderr, catch_backtrace())
                        end
                    end
                end
                push!(tasks, task)
            end
        end

        CImGui.End()
    end
end
