Base.@kwdef mutable struct ProblemDesc
    name::String
    params = []
    fn = nothing
    bounds = nothing
    pfront = nothing
    pfront_image = nothing
end

function problem_gui(problem_desc)
    problem_pfront_image_id = nothing

    function ()
        problems = ["ZDT1", "ZDT2", "ZDT3", "ZDT4", "ZDT6",
            "DTLZ1", "DTLZ2", "DTLZ3", "DTLZ4", "DTLZ5", "DTLZ6"]
        @cstatic item_current = "ZDT1" begin
            if CImGui.BeginCombo("Problem", item_current)
                for n = 0:length(problems)-1
                    is_selected = item_current == problems[n+1]
                    if CImGui.Selectable(problems[n+1], is_selected)
                        item_current = problems[n+1]
                        update_problem_desc(problem_desc, problems[n+1])
                    end
                    # set the initial focus when opening the combo (scrolling + for keyboard navigation support in the upcoming navigation branch)
                    is_selected && CImGui.SetItemDefaultFocus()
                end
                CImGui.EndCombo()
            end
        end
        for param_desc in problem_desc.params
            eval(param_gui(problem_desc, param_desc))
        end
        if problem_pfront_image_id !== nothing
            ImGui_ImplOpenGL3_DestroyImageTexture(problem_pfront_image_id)
        end
        if problem_desc.pfront_image !== nothing
            img = problem_desc.pfront_image
            _, img_width, img_height = size(img)
            problem_pfront_image_id = ImGui_ImplOpenGL3_CreateImageTexture(img_width, img_height)
            ImGui_ImplOpenGL3_UpdateImageTexture(problem_pfront_image_id, img, img_width, img_height)
            CImGui.Image(Ptr{Cvoid}(problem_pfront_image_id), (img_width, img_height))
        end
    end
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
            size = get_figure_size()
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
