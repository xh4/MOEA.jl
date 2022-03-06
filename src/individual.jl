mutable struct Individual
    variable::AbstractVector
    objective::Union{AbstractVector, Number, Nothing}
    constraint::Union{AbstractVector, Nothing}

    Individual(var) = new(var, nothing, nothing)
    Individual(var, obj) = new(var, obj, nothing)
    Individual(var, obj, cst) = new(var, obj, cst)
end

variable(i::Individual) = i.variable

objective(i::Individual) = i.objective

constraint(i::Individual) = i.constraint

copy(i::Individual) = Individual(i.variable, i.objective, i.constraint)