Base.@kwdef mutable struct ProblemDesc
    name::String
    params = []
    problem = nothing
    truepf_image = nothing
end

function problem_gui(problem_desc)
    problem_truepf_image_id = nothing

    function ()
        problems = ["ZDT1", "ZDT2", "ZDT3", "ZDT4", "ZDT6",
                    "DTLZ1", "DTLZ2", "DTLZ3", "DTLZ4", "DTLZ5", "DTLZ6",
                    "WFG1", "WFG2", "WFG3", "WFG4", "WFG5", "WFG6", "WFG7", "WFG8", "WFG9"]
        @cstatic item_current = "" begin
            item_current = problem_desc.name
            if CImGui.BeginCombo("Problem", item_current)
                for n = 0:length(problems)-1
                    is_selected = item_current == problems[n+1]
                    if CImGui.Selectable(problems[n+1], is_selected)
                        item_current = problems[n+1]
                        update_problem_param(problem_desc, problems[n+1])
                    end
                    # set the initial focus when opening the combo (scrolling + for keyboard navigation support in the upcoming navigation branch)
                    is_selected && CImGui.SetItemDefaultFocus()
                end
                CImGui.EndCombo()
            end
        end
        for param_desc in problem_desc.params
            param_gui(problem_desc, param_desc)
        end
        if problem_truepf_image_id !== nothing
            DestroyImageTexture(problem_truepf_image_id)
        end
        if problem_desc.truepf_image !== nothing
            img = problem_desc.truepf_image
            _, img_width, img_height = size(img)
            problem_truepf_image_id = CreateImageTexture(img_width, img_height)
            UpdateImageTexture(problem_truepf_image_id, img, img_width, img_height)
            CImGui.Image(Ptr{Cvoid}(problem_truepf_image_id), (img_width, img_height))
        end
    end
end

function get_problem_desc(problem_name)
    common_params = [ParamDesc(name = "maxFE", type = "Int",
                               value = 10_000, default_value = 10_000,
                               desc = "Max Function Evaluations",
                               min = 100, max = 10_000_000_000)]
    if problem_name in ["ZDT1", "ZDT2", "ZDT3", "ZDT4", "ZDT6"]
        problem_desc = ProblemDesc(name = problem_name,
                                   params = [[ParamDesc(name = "D", type = "Int",
                                                        value = 30, default_value = 30,
                                                        desc = "Dimension",
                                                        min = 2, max = 100)]; 
                                             common_params])
    elseif problem_name in ["DTLZ1", "DTLZ2", "DTLZ3", "DTLZ4", "DTLZ5", "DTLZ6"]
        problem_desc = ProblemDesc(name = problem_name,
                                   params = [[ParamDesc(name = "M", type = "Int",
                                                        value = 3, default_value = 3,
                                                        desc = "Dimension",
                                                        min = 3, max = 100)];
                                             common_params])
    elseif problem_name in ["WFG1", "WFG2", "WFG3", "WFG4", "WFG5", "WFG6", "WFG7", "WFG8", "WFG9"]
        problem_desc = ProblemDesc(name = problem_name,
                                   params = [[ParamDesc(name = "M", type = "Int",
                                                        value = 3, default_value = 3,
                                                        desc = "Dimension",
                                                        min = 3, max = 100)];
                                             common_params])
    end
    update_problem_desc(problem_desc)
end

function get_problem_param(problem_desc, param_name)
    for param in problem_desc.params
        if param.name == param_name
            return param.value === nothing ? param.default_value : param.value
        end
    end
end

function update_problem_param(problem_desc, problem_name, param_name=nothing, param_value=nothing)
    # println("update_problem_desc ", problem_name, " ", param_name, " ", param_value)
    should_update = problem_desc.problem === nothing ? true : false
    if problem_name !== nothing && problem_desc.name != problem_name
        new_problem_desc = get_problem_desc(problem_name)
        problem_desc.name = new_problem_desc.name
        problem_desc.params = new_problem_desc.params
        should_update = true
    else
        for param in problem_desc.params
            if param.name == param_name
                if param.value != param_value
                    should_update = true
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
    if should_update
        update_problem_desc(problem_desc)
    end
    problem_desc
end

function update_problem_desc(problem_desc)
    maxFE = get_problem_param(problem_desc, "maxFE")
    if problem_desc.name in ["ZDT1", "ZDT2", "ZDT3", "ZDT4", "ZDT6"]
        D = get_problem_param(problem_desc, "D")
        problem = eval(Expr(:call, Symbol(problem_desc.name),
                            D, Expr(:kw, :maxFE, maxFE)))
    elseif problem_desc.name in ["DTLZ1", "DTLZ2", "DTLZ3", "DTLZ4", "DTLZ5", "DTLZ6"]
        M = get_problem_param(problem_desc, "M")
        problem = eval(Expr(:call, Symbol(problem_desc.name),
                            M, Expr(:kw, :maxFE, maxFE)))
    elseif problem_desc.name in ["WFG1", "WFG2", "WFG3", "WFG4", "WFG5", "WFG6", "WFG7", "WFG8", "WFG9"]
        M = get_problem_param(problem_desc, "M")
        problem = eval(Expr(:call, Symbol(problem_desc.name),
                            M, Expr(:kw, :maxFE, maxFE)))
    end
    problem_desc.problem = problem

    fig_size = get_figure_size()
    if 2 <= problem.M <= 3
        if problem.M == 2
            truepf_fig = plot(
                objectives(problem.truepf[1:100:end])[:, 1],
                objectives(problem.truepf[1:100:end])[:, 2],
                seriestype = :scatter,
                legend = false,
                width = fig_size,
                height = fig_size,
                color = :blue
            )
        elseif problem.M == 3
            truepf_fig = plot(
                objectives(problem.truepf[1:100:end])[:, 1],
                objectives(problem.truepf[1:100:end])[:, 2],
                objectives(problem.truepf[1:100:end])[:, 3],
                seriestype = :scatter,
                legend = false,
                width = fig_size,
                height = fig_size,
                camera = (30,45),
                color = :blue
            )
        end
        truepf_image = plot_image(truepf_fig)
        truepf_image = fit_image(truepf_image, fig_size)
        problem_desc.truepf_image = image_uint8(truepf_image)
    end
    
    problem_desc
end
