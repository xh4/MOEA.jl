abstract type AbstractOptimizer end

function print_header(method::AbstractOptimizer)
    println("Iter     Function value")
end
population_size(method::AbstractOptimizer) =
    error("`population_size` is not implemented for $(summary(method)).")
metrics(method::AbstractOptimizer) = method.metrics
