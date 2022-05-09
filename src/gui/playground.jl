function playground_window()
    problem = get_problem_desc("ZDT1")
    problem_renderer = problem_gui(problem)

    algorithm = get_algorithm_desc("MOEA/D")
    algorithm_renderer = algorithm_gui(algorithm)

    runtime_image_id = nothing
    runtime_image = nothing
    generation = 0
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
            igds = []
            hvs = []
            constraints = BoxConstraints(problem.bounds)
            fig_size = get_figure_size()
            state_callback = function (state)
                try
                    generation = state.iteration

                    fig = Plots.plot(
                        objectives(pfront(state))[:, 1],
                        objectives(pfront(state))[:, 2],
                        seriestype = :scatter,
                        legend = false,
                    )
                    runtime_image = fit_image(plot_image(fig), fig_size)

                    igd = MOEA.igd(pfront(state), problem.pfront)
                    push!(igds, igd)
                    fig = Plots.plot(
                        1:length(igds),
                        igds,
                        legend = false,
                    )
                    igd_image = fit_image(plot_image(fig), fig_size)

                    hv = hypervolume(objectives(pfront(state)), ones(size(objectives(pfront(state)))[2]))
                    push!(hvs, hv)
                    fig = Plots.plot(
                        1:length(hvs),
                        hvs,
                        legend = false,
                    )
                    hv_image = fit_image(plot_image(fig), fig_size)
                catch e
                    @error "Error in state callback!" exception=e
                    Base.show_backtrace(stderr, catch_backtrace())
                end
            end
            # TODO: initial population
            population_size = get_algorithm_param(algorithm, "N")
            Threads.@spawn optimize(problem.fn,
                                    algorithm.method,
                                    constraints = constraints,
                                    population = [Individual(rand(30)) for i in 1:population_size],
                                    options = Options(state_callback = state_callback))
        end

        @cstatic v = Cint(0) begin
            v = Int32(generation)
            @c CImGui.InputInt("Generation", &v)
        end
        @cstatic v = Cfloat(0.0) begin
            v = Float32(igd)
            @c CImGui.InputFloat("IGD", &v)
        end
        @cstatic v = Cfloat(0.0) begin
            v = Float32(hv)
            @c CImGui.InputFloat("HV", &v)
        end

        if runtime_image !== nothing
            if runtime_image_id !== nothing
                ImGui_ImplOpenGL3_DestroyImageTexture(runtime_image_id)
            end
            img_width, img_height = size(runtime_image)
            runtime_image_id = ImGui_ImplOpenGL3_CreateImageTexture(img_width, img_height)
            ImGui_ImplOpenGL3_UpdateImageTexture(runtime_image_id, image_uint8(runtime_image), img_width, img_height)
            CImGui.Image(Ptr{Cvoid}(runtime_image_id), (img_width, img_height))
        end

        if igd_image !== nothing && length(igds) > 0
            if igd_image_id !== nothing
                ImGui_ImplOpenGL3_DestroyImageTexture(igd_image_id)
            end
            img_width, img_height = size(igd_image)
            igd_image_id = ImGui_ImplOpenGL3_CreateImageTexture(img_width, img_height)
            ImGui_ImplOpenGL3_UpdateImageTexture(igd_image_id, image_uint8(igd_image), img_width, img_height)
            CImGui.Image(Ptr{Cvoid}(igd_image_id), (img_width, img_height))
        end

        if hv_image !== nothing && length(hvs) > 0
            if hv_image_id !== nothing
                ImGui_ImplOpenGL3_DestroyImageTexture(hv_image_id)
            end
            img_width, img_height = size(hv_image)
            hv_image_id = ImGui_ImplOpenGL3_CreateImageTexture(img_width, img_height)
            ImGui_ImplOpenGL3_UpdateImageTexture(hv_image_id, image_uint8(hv_image), img_width, img_height)
            CImGui.Image(Ptr{Cvoid}(hv_image_id), (img_width, img_height))
        end

        CImGui.End()
    end
end
