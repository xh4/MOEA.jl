Base.@kwdef mutable struct ParamDesc
    name::String
    type::String
    default_value
    value
    desc::String = ""
    min = nothing
    max = nothing
end

function param_gui(desc, param_desc)
    if param_desc.type == "Int"
        @cstatic v = Cint(0) begin
            v = Int32(param_desc.value)
            @c CImGui.InputInt(param_desc.name, &v)
            # if $param_desc.value === nothing
            #     v = Int32($param_desc.default_value)
            # end
            if isa(desc, ProblemDesc)
                update_problem_param(desc, nothing, param_desc.name, v)
            elseif isa(desc, AlgorithmDesc)
                update_algorithm_param(desc, nothing, param_desc.name, v)
            end
            CImGui.SameLine()
            HelpMarker(param_desc.desc)
        end
    end
end
