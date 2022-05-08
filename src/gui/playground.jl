Base.@kwdef mutable struct ParamDesc
    name::String
    type::String
    default_value
    value
    desc::String = ""
    min = nothing
    max = nothing
end

Base.@kwdef mutable struct ProblemDesc
    name::String
    params = []
    fn = nothing
    bounds = nothing
    pfront = nothing
    pfront_image = nothing
end

Base.@kwdef mutable struct AlgorithmDesc
    name::String
    params = []
    method = nothing
end

function get_problem_desc(problem_name)
    ProblemDesc(name = problem_name,
                params = [ParamDesc(name="D", type="Int",
                                    value=30, default_value=30,
                                    desc="Dimension",
                                    min=2, max=100),
                          ParamDesc(name="n", type="Int",
                                    value=100, default_value=100,
                                    desc="Number of solutions",
                                    min=2, max=1000)])
end

function get_algorithm_desc(algorithm_name)
    if algorithm_name == "MOEA/D"
        AlgorithmDesc(name = algorithm_name,
                      params = [ParamDesc(name="N", type="Int",
                                          value=100, default_value=100,
                                          desc="Population size",
                                          min=2, max=10000),
                                ParamDesc(name="T", type="Int",
                                          value=10, default_value=10,
                                          desc="Number of solutions",
                                          min=2, max=100)])
    elseif algorithm_name == "NSGA-II"
        AlgorithmDesc(name = algorithm_name,
                      params = [ParamDesc(name = "N", type = "Int",
                                          value = 100, default_value = 100,
                                          desc = "Population size",
                                          min = 2, max = 100)])
    end
end

function update_algorithm_desc(algorithm_desc, algorithm_name, param_name=nothing, param_value=nothing)
    should_update = algorithm_desc.method === nothing ? true : false
    if algorithm_name !== nothing && algorithm_desc.name != algorithm_name
        new_algorithm_desc = get_algorithm_desc(algorithm_name)
        algorithm_desc.name = new_algorithm_desc.name
        algorithm_desc.params = new_algorithm_desc.params
        should_update = true
    else
        for param in algorithm_desc.params
            if param.name == param_name
                if param.max !== nothing && param_value > param.max
                        param_value = param.max
                    end
                    if param.min !== nothing && param_value < param.min
                        param_value = param.min
                    end
                if param.value != param_value
                    param.value = param_value
                    should_update = true
                end
            end
        end
    end
    if should_update
        if algorithm_desc.name == "MOEA/D"
            N = get_algorithm_param(algorithm_desc, "N")
            T = get_algorithm_param(algorithm_desc, "T")
            algorithm_desc.method = MOEAD(N=N, T=T)
        elseif algorithm_desc.name == "NSGA-II"
            N = get_algorithm_param(algorithm_desc, "N")
            algorithm_desc.method = NSGAII(populationSize=N)
        end
    end
end

function get_algorithm_param(algorithm_desc, param_name)
    for param in algorithm_desc.params
        if param.name == param_name
            return param.value === nothing ? param.default_value : param.value
        end
    end
end

function get_problem_param(problem_desc, param_name)
    for param in problem_desc.params
        if param.name == param_name
            return param.value === nothing ? param.default_value : param.value
        end
    end
end

function update_problem_desc(problem_desc, problem_name, param_name=nothing, param_value=nothing)
    # println("update_problem_desc ", problem_name, " ", param_name, " ", param_value)
    should_update_fn = problem_desc.fn === nothing ? true : false
    if problem_name !== nothing && problem_desc.name != problem_name
        new_problem_desc = get_problem_desc(problem_name)
        problem_desc.name = new_problem_desc.name
        problem_desc.params = new_problem_desc.params
        should_update_fn = true
    else
        for param in problem_desc.params
            if param.name == param_name
                if param.value != param_value
                    should_update_fn = true
                    if param.max !== nothing && param_value > param.max
                        param_value = param.max
                    end
                    if param.min !== nothing && param_value < param.min
                        param_value = param.min
                    end
                    param.value = param_value
                end
            end
        end
    end
    if should_update_fn
        if problem_desc.name in ["ZDT1", "ZDT2", "ZDT3", "ZDT4", "ZDT6"]
            fn, bounds, pfront = eval(Expr(:call, Symbol(problem_desc.name),
                                           get_problem_param(problem_desc, "D"),
                                           get_problem_param(problem_desc, "n")))
            problem_desc.fn = fn
            problem_desc.bounds = bounds
            problem_desc.pfront = pfront
            size = CImGui.GetWindowWidth() / 5 * 4
            pfront_fig = plot(
                objectives(pfront)[:, 1],
                objectives(pfront)[:, 2],
                seriestype = :scatter,
                legend = false,
                width = size,
                height = size
            )
            pfront_image = plot_image(pfront_fig)
            pfront_image = fit_image(pfront_image, size)
            problem_desc.pfront_image = image_uint8(pfront_image)
        end
    end
end

function param_gui(desc, param_desc)
    quote
        if $param_desc.type == "Int"
            @cstatic v = Cint($param_desc.value) begin
                @c CImGui.InputInt($param_desc.name, &v)
                # if $param_desc.value === nothing
                #     v = Int32($param_desc.default_value)
                # end
                if isa($desc, ProblemDesc)
                    update_problem_desc($desc, nothing, $param_desc.name, v)
                elseif isa($desc, AlgorithmDesc)
                    update_algorithm_desc($desc, nothing, $param_desc.name, v)
                end
                CImGui.SameLine()
                HelpMarker($param_desc.desc)
            end
        end
    end
end

function test_image(img_width, img_height)
    rand(GLubyte, 4, img_width, img_height)
end

function load_image()
    img = load("/mnt/c/Users/XH/Desktop/AwesomeFace.png")'
    img = coloralpha.(img)
    view = channelview(img)
    reinterpret(UInt8, view)
end

function plot_image(fig)
    io = IOBuffer()
    show(io, MIME("image/png"), fig)
    load(io)'
end

function image_uint8(img)
    img = coloralpha.(img)
    view = channelview(img)
    reinterpret(UInt8, view)
end

function fit_image(img, size)
    size = Int(round(size))
    width, height = Base.size(img)
    if width > size
        new_width = size
        new_height = Int(round(size/width*height))
        img = imresize(img, (new_width, new_height))
    end
    img
end

function playground_window()
    problem = get_problem_desc("ZDT1")
    problem_pfront_image_id = nothing

    algorithm = get_algorithm_desc("MOEA/D")

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
        problems = ["ZDT1", "ZDT2", "ZDT3", "ZDT4", "ZDT6",
                    "DTLZ1", "DTLZ2", "DTLZ3", "DTLZ4", "DTLZ5", "DTLZ6"]
        @cstatic item_current = "ZDT1" begin
            if CImGui.BeginCombo("Problem", item_current)
                for n = 0:length(problems)-1
                    is_selected = item_current == problems[n+1]
                    if CImGui.Selectable(problems[n+1], is_selected)
                        item_current = problems[n+1]
                        update_problem_desc(problem, problems[n+1])
                    end
                    # set the initial focus when opening the combo (scrolling + for keyboard navigation support in the upcoming navigation branch)
                    is_selected && CImGui.SetItemDefaultFocus()
                end
                CImGui.EndCombo()
            end
        end
        for param_desc in problem.params
            eval(param_gui(problem, param_desc))
        end
        if problem_pfront_image_id !== nothing
            ImGui_ImplOpenGL3_DestroyImageTexture(problem_pfront_image_id)
        end
        if problem.pfront_image !== nothing
            img = problem.pfront_image
            _, img_width, img_height = size(img)
            problem_pfront_image_id = ImGui_ImplOpenGL3_CreateImageTexture(img_width, img_height)
            ImGui_ImplOpenGL3_UpdateImageTexture(problem_pfront_image_id, img, img_width, img_height)
            CImGui.Image(Ptr{Cvoid}(problem_pfront_image_id), (img_width, img_height))
        end

        algorithms = ["NSGA-II", "MOEA/D"]
        @cstatic item_current = "MOEA/D" begin
            if CImGui.BeginCombo("Algorithm", item_current)
                for n = 0:length(algorithms)-1
                    is_selected = item_current == algorithms[n+1]
                    if CImGui.Selectable(algorithms[n+1], is_selected)
                        item_current = algorithms[n+1]
                        update_algorithm_desc(algorithm, algorithms[n+1])
                    end
                    # set the initial focus when opening the combo (scrolling + for keyboard navigation support in the upcoming navigation branch)
                    is_selected && CImGui.SetItemDefaultFocus()
                end
                CImGui.EndCombo()
            end
        end
        for param_desc in algorithm.params
            eval(param_gui(algorithm, param_desc))
        end

        if CImGui.Button("Evaluate")
            igds = []
            hvs = []
            constraints = BoxConstraints(problem.bounds)
            state_callback = function (state)
                try
                    generation = state.iteration

                    fig = Plots.plot(
                        objectives(pfront(state))[:, 1],
                        objectives(pfront(state))[:, 2],
                        seriestype = :scatter,
                        legend = false,
                    )
                    runtime_image = plot_image(fig)

                    igd = MOEA.igd(pfront(state), problem.pfront)
                    push!(igds, igd)
                    fig = Plots.plot(
                        1:length(igds),
                        igds,
                        legend = false,
                    )
                    igd_image = plot_image(fig)

                    hv = hypervolume(objectives(pfront(state)), ones(size(objectives(pfront(state)))[2]))
                    push!(hvs, hv)
                    fig = Plots.plot(
                        1:length(hvs),
                        hvs,
                        legend = false,
                    )
                    hv_image = plot_image(fig)
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
