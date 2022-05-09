Base.@kwdef mutable struct AlgorithmDesc
    name::String
    params = []
    method = nothing
end

function algorithm_gui(algorithm_desc)
    function ()
        algorithms = ["NSGA-II", "MOEA/D"]
        @cstatic item_current = "MOEA/D" begin
            if CImGui.BeginCombo("Algorithm", item_current)
                for n = 0:length(algorithms)-1
                    is_selected = item_current == algorithms[n+1]
                    if CImGui.Selectable(algorithms[n+1], is_selected)
                        item_current = algorithms[n+1]
                        update_algorithm_desc(algorithm_desc, algorithms[n+1])
                    end
                    # set the initial focus when opening the combo (scrolling + for keyboard navigation support in the upcoming navigation branch)
                    is_selected && CImGui.SetItemDefaultFocus()
                end
                CImGui.EndCombo()
            end
        end
        for param_desc in algorithm_desc.params
            eval(param_gui(algorithm_desc, param_desc))
        end
    end
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
