abstract type AbstractIndividual end

variables(i::AbstractIndividual) = error("Function `variables` not implemented")

objectives(i::AbstractIndividual) = error("Function `objectives` not implemented")

mutable struct Individual <: AbstractIndividual
    variables::AbstractVector
    objectives::Union{AbstractVector, Number, Nothing}

    Individual(var) = new(var, nothing)
    Individual(var, obj) = new(var, obj)
end

variables(i::Individual) = i.variables

objectives(i::Individual) = i.objectives

copy(n::Nothing) = nothing

copy(i::Individual) = Individual(copy(i.variables), copy(i.objectives))