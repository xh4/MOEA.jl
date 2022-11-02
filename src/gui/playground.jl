function playground_window()
    problem_desc = get_setting("gui.playground.problem", get_problem_desc("ZDT1")) 
    problem_renderer = problem_gui(problem_desc)
    
    algorithm_desc = get_setting("gui.playground.algorithm", get_algorithm_desc("NSGA-II"))
    algorithm_renderer = algorithm_gui(algorithm_desc)

    runtime_image_id = nothing
    runtime_image = nothing
    generation = 0
    fcalls = 0
    runtime = 0.0
    igd = 0.0
    igds::Vector{Float64} = []
    igd_image_id = nothing
    igd_image = nothing
    hv = 0.0
    hvs::Vector{Float64} = []
    hv_image_id = nothing
    hv_image = nothing

    function ()
        CImGui.Begin("Playground")

        problem_renderer()

        algorithm_renderer()

        if CImGui.Button("Evaluate")
            set_setting("gui.playground.problem", problem_desc)
            set_setting("gui.playground.algorithm", algorithm_desc)

            problem = problem_desc.problem
            algorithm = algorithm_desc.algorithm

            runtime = 0
            igds = []
            hvs = []
            fig_size = get_figure_size()
            calc_metrics = function(state)
                igd = IGD(pfront(state), problem.truepf)
                push!(igds, igd)
                fig = Plots.plot(
                    1:length(igds),
                    igds,
                    legend = false,
                )
                igd_image = fit_image(plot_image(fig), fig_size)

                hv = HV(pfront(state), problem.truepf)
                push!(hvs, hv)
                fig = Plots.plot(
                    1:length(hvs),
                    hvs,
                    legend = false,
                )
                hv_image = fit_image(plot_image(fig), fig_size)
            end
            state_callback = function(state)
                try
                    generation = state.iteration
                    fcalls = state.fcalls
                    runtime += (state.stop_time - state.start_time)

                    if problem.M == 2
                        fig = Plots.plot(
                            objectives(problem.truepf[1:100:end])[:, 1],
                            objectives(problem.truepf[1:100:end])[:, 2],
                            seriestype = :scatter,
                            color = :blue,
                            legend = false
                        )
                        Plots.plot!(
                            objectives(pfront(state))[:, 1],
                            objectives(pfront(state))[:, 2],
                            seriestype = :scatter
                        )
                    elseif problem.M == 3
                        fig = Plots.plot(
                            objectives(problem.truepf[1:100:end])[:, 1],
                            objectives(problem.truepf[1:100:end])[:, 2],
                            objectives(problem.truepf[1:100:end])[:, 3],
                            seriestype = :scatter,
                            color = :blue,
                            legend = false,
                            camera = (30,60)
                        )
                        Plots.plot!(
                            objectives(pfront(state))[:, 1],
                            objectives(pfront(state))[:, 2],
                            objectives(pfront(state))[:, 3],
                            seriestype = :scatter
                        )
                    end
                    runtime_image = fit_image(plot_image(fig), fig_size)

                    if state.fcalls % 500 == 0
                        calc_metrics(state)
                    end
                catch e
                    @error "Error in state callback!" exception=e
                    Base.show_backtrace(stderr, catch_backtrace())
                end
            end
            finish_callback = function(state)
                calc_metrics(state)
            end
            Threads.@spawn try
                optimize(problem,
                         algorithm,
                         options = Options(state_callback = state_callback,
                                           finish_callback = finish_callback))
            catch e
                @error "Error when optimize" exception=e
                Base.show_backtrace(stderr, catch_backtrace())
            end
        end

        @cstatic v = Cint(0) begin
            v = Int32(generation)
            @c CImGui.InputInt("Generation", &v)
        end
        @cstatic v = Cint(0) begin
            v = Int32(fcalls)
            @c CImGui.InputInt("Evaluations", &v)
        end
        @cstatic v = Cfloat(0.0) begin
            v = Float32(runtime)
            @c CImGui.InputFloat("Runtime", &v, 0, 0, "%.3f")
        end
        @cstatic v = Cfloat(0.0) begin
            v = Float32(igd)
            @c CImGui.InputFloat("IGD", &v, 0, 0, "%.6f")
        end
        @cstatic v = Cfloat(0.0) begin
            v = Float32(hv)
            @c CImGui.InputFloat("HV", &v, 0, 0, "%.6f")
        end

        if runtime_image !== nothing
            if runtime_image_id !== nothing
                DestroyImageTexture(runtime_image_id)
            end
            img_width, img_height = size(runtime_image)
            runtime_image_id = CreateImageTexture(img_width, img_height)
            UpdateImageTexture(runtime_image_id, image_uint8(runtime_image), img_width, img_height)
            CImGui.Image(Ptr{Cvoid}(runtime_image_id), (img_width, img_height))
        end

        if igd_image !== nothing && length(igds) > 0
            if igd_image_id !== nothing
                DestroyImageTexture(igd_image_id)
            end
            img_width, img_height = size(igd_image)
            igd_image_id = CreateImageTexture(img_width, img_height)
            UpdateImageTexture(igd_image_id, image_uint8(igd_image), img_width, img_height)
            CImGui.Image(Ptr{Cvoid}(igd_image_id), (img_width, img_height))
        end

        if hv_image !== nothing && length(hvs) > 0
            if hv_image_id !== nothing
                DestroyImageTexture(hv_image_id)
            end
            img_width, img_height = size(hv_image)
            hv_image_id = CreateImageTexture(img_width, img_height)
            UpdateImageTexture(hv_image_id, image_uint8(hv_image), img_width, img_height)
            CImGui.Image(Ptr{Cvoid}(hv_image_id), (img_width, img_height))
        end

        CImGui.End()
    end
end
