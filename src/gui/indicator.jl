Base.@kwdef mutable struct IndicatorDesc
    name::String
    better = <
    mean::Bool = true
    min::Bool = true
    max::Bool = true
    std::Bool = true
end

function igd_indicator_desc()
    IndicatorDesc(name="IGD")
end

function hv_indicator_desc()
    IndicatorDesc(name="HV", better= >)
end
