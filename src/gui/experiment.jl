function experiment_window()

    algorithms = [get_algorithm_desc("MOEA/D")]
    problems = [get_problem_desc("ZDT1")]

    function()
        CImGui.Begin("Experiment")



        if CImGui.Button("+ Algorithm")

        end

        if CImGui.Button("+ Problem")

        end

        CImGui.End()
    end
end
