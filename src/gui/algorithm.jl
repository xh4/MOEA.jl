Base.@kwdef mutable struct AlgorithmDesc
    name::String
    params = []
    algorithm = nothing
end

function algorithm_gui(algorithm_desc)
    algorithms = ["NSGA-II", "MOEA/D", "MOEA/D-EP", "MOEA/D-EPHV", "MOEA/D-DE", "MOEA/D-DRA", "SMS-EMOA", "MOEA/D-MY"]

    function ()
        @cstatic item_current = "MOEA/D" begin
            if CImGui.BeginCombo("Algorithm", item_current)
                for n = 0:length(algorithms)-1
                    is_selected = item_current == algorithms[n+1]
                    if CImGui.Selectable(algorithms[n+1], is_selected)
                        item_current = algorithms[n+1]
                        update_algorithm_param(algorithm_desc, algorithms[n+1])
                    end
                    # set the initial focus when opening the combo (scrolling + for keyboard navigation support in the upcoming navigation branch)
                    is_selected && CImGui.SetItemDefaultFocus()
                end
                CImGui.EndCombo()
            end
        end
        for param_desc in algorithm_desc.params
            param_gui(algorithm_desc, param_desc)
        end
    end
end

function get_algorithm_desc(algorithm_name)
    if algorithm_name == "NSGA-II"
        algorithm_desc = AlgorithmDesc(name = algorithm_name,
                                       params = [ParamDesc(name = "N", type = "Int",
                                                           value = 100, default_value = 100,
                                                           desc = "Population size",
                                                           min = 2, max = 10000)])
    elseif algorithm_name == "MOEA/D"
        algorithm_desc = AlgorithmDesc(name = algorithm_name,
                                       params = [ParamDesc(name = "N", type = "Int",
                                                           value = 100, default_value = 100,
                                                           desc = "Population size",
                                                           min = 2, max = 10000),
                                                 ParamDesc(name = "T", type = "Int",
                                                           value = 10, default_value = 10,
                                                           desc = "Number of solutions",
                                                           min = 2, max = 100)])
    elseif algorithm_name == "MOEA/D"
        algorithm_desc = AlgorithmDesc(name = algorithm_name,
                                       params = [ParamDesc(name = "N", type = "Int",
                                                           value = 100, default_value = 100,
                                                           desc = "Population size",
                                                           min = 2, max = 10000),
                                                 ParamDesc(name = "T", type = "Int",
                                                           value = 10, default_value = 10,
                                                           desc = "Number of solutions",
                                                           min = 2, max = 100)])
    elseif algorithm_name == "MOEA/D-EP"
        algorithm_desc = AlgorithmDesc(name = algorithm_name,
                                        params = [ParamDesc(name = "N", type = "Int",
                                                            value = 100, default_value = 100,
                                                            desc = "Population size",
                                                            min = 2, max = 10000),
                                                    ParamDesc(name = "T", type = "Int",
                                                            value = 10, default_value = 10,
                                                            desc = "Number of solutions",
                                                            min = 2, max = 100)])
    elseif algorithm_name == "MOEA/D-EPHV"
        algorithm_desc = AlgorithmDesc(name = algorithm_name,
                                        params = [ParamDesc(name = "N", type = "Int",
                                                            value = 100, default_value = 100,
                                                            desc = "Population size",
                                                            min = 2, max = 10000),
                                                    ParamDesc(name = "T", type = "Int",
                                                            value = 10, default_value = 10,
                                                            desc = "Number of solutions",
                                                            min = 2, max = 100)])
    elseif algorithm_name == "MOEA/D-DE"
        algorithm_desc = AlgorithmDesc(name = algorithm_name,
                                       params = [ParamDesc(name = "N", type = "Int",
                                                           value = 100, default_value = 100,
                                                           desc = "Population size",
                                                           min = 2, max = 10000),
                                                 ParamDesc(name = "T", type = "Int",
                                                           value = 10, default_value = 10,
                                                           desc = "Number of solutions",
                                                           min = 2, max = 100)])
    elseif algorithm_name == "MOEA/D-DRA"
        algorithm_desc = AlgorithmDesc(name = algorithm_name,
                                       params = [ParamDesc(name = "N", type = "Int",
                                                           value = 100, default_value = 100,
                                                           desc = "Population size",
                                                           min = 2, max = 10000),
                                                 ParamDesc(name = "T", type = "Int",
                                                           value = 10, default_value = 10,
                                                           desc = "Number of solutions",
                                                           min = 2, max = 100)])
    elseif algorithm_name == "SMS-EMOA"
        algorithm_desc = AlgorithmDesc(name = algorithm_name,
                                        params = [ParamDesc(name = "N", type = "Int",
                                                            value = 100, default_value = 100,
                                                            desc = "Population size",
                                                            min = 2, max = 10000)])
    elseif algorithm_name == "MOEA/D-MY"
        algorithm_desc = AlgorithmDesc(name = algorithm_name,
                                        params = [ParamDesc(name = "N", type = "Int",
                                                            value = 100, default_value = 100,
                                                            desc = "Population size",
                                                            min = 2, max = 10000),
                                                    ParamDesc(name = "T", type = "Int",
                                                            value = 10, default_value = 10,
                                                            desc = "Number of solutions",
                                                            min = 2, max = 100)])
    end
    update_algorithm_desc(algorithm_desc)
end

function update_algorithm_param(algorithm_desc, algorithm_name, param_name=nothing, param_value=nothing)
    should_update = algorithm_desc.algorithm === nothing ? true : false
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
        update_algorithm_desc(algorithm_desc)
    end
    algorithm_desc
end

function get_algorithm_param(algorithm_desc, param_name)
    for param in algorithm_desc.params
        if param.name == param_name
            return param.value === nothing ? param.default_value : param.value
        end
    end
end

function update_algorithm_desc(algorithm_desc)
    if algorithm_desc.name == "NSGA-II"
        N = get_algorithm_param(algorithm_desc, "N")
        algorithm_desc.algorithm = NSGAII(N = N)
    elseif algorithm_desc.name == "MOEA/D"
        N = get_algorithm_param(algorithm_desc, "N")
        T = get_algorithm_param(algorithm_desc, "T")
        algorithm_desc.algorithm = MOEAD(N = N, T = T)
    elseif algorithm_desc.name == "MOEA/D-EP"
        N = get_algorithm_param(algorithm_desc, "N")
        T = get_algorithm_param(algorithm_desc, "T")
        algorithm_desc.algorithm = MOEADEP(N = N, T = T)
    elseif algorithm_desc.name == "MOEA/D-EPHV"
        N = get_algorithm_param(algorithm_desc, "N")
        T = get_algorithm_param(algorithm_desc, "T")
        algorithm_desc.algorithm = MOEADEPHV(N = N, T = T)
    elseif algorithm_desc.name == "MOEA/D-DE"
        N = get_algorithm_param(algorithm_desc, "N")
        T = get_algorithm_param(algorithm_desc, "T")
        algorithm_desc.algorithm = MOEADDE(N = N, T = T)
    elseif algorithm_desc.name == "MOEA/D-DRA"
        N = get_algorithm_param(algorithm_desc, "N")
        T = get_algorithm_param(algorithm_desc, "T")
        algorithm_desc.algorithm = MOEADDRA(N = N, T = T)
    elseif algorithm_desc.name == "SMS-EMOA"
        N = get_algorithm_param(algorithm_desc, "N")
        algorithm_desc.algorithm = SMSEMOA(N = N)
    elseif algorithm_desc.name == "MOEA/D-MY"
        N = get_algorithm_param(algorithm_desc, "N")
        T = get_algorithm_param(algorithm_desc, "T")
        algorithm_desc.algorithm = MOEADMY(N = N, T = T)
    end
    algorithm_desc
end
