Base.@kwdef mutable struct Parameter
    type::DataType
    name::String
    value::Any
end
