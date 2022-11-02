abstract type AbstractAlgorithm end

show(io::IO, method::AbstractAlgorithm) = print(io, identifier(method))

function identifier(a::AbstractAlgorithm)
    join([string(typeof(a).name.name), "(", properties_string(a), ")"])
end

function properties_string(s)
    properties_string(s, propertynames(s))
end

function properties_string(s, propertynames)
    join(map(name->"$(name)=$(getproperty(s, name))", propertynames), ",")
end

population_size(method::AbstractAlgorithm) = method.N
