abstract type Optimizer end

function print_header(method::Optimizer)
    println("Iter     Function value")
end
population_size(method::Optimizer) =
    error("`population_size` is not implemented for $(summary(method)).")
metrics(method::Optimizer) = method.metrics
