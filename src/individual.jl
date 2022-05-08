abstract type AbstractIndividual end

variables(i::AbstractIndividual) = error("Function `variables` not implemented")

objectives(i::AbstractIndividual) = error("Function `objectives` not implemented")

constraints(i::AbstractIndividual) = error("Function `constraints` not implemented")

mutable struct Individual <: AbstractIndividual
    variables::AbstractVector
    objectives::Union{AbstractVector, Number, Nothing}
    constraints::Union{AbstractVector, Nothing}

    Individual(var) = new(var, nothing, nothing)
    Individual(var, obj) = new(var, obj, nothing)
    Individual(var, obj, cst) = new(var, obj, cst)
end

variables(i::Individual) = i.variables

objectives(i::Individual) = i.objectives

constraints(i::Individual) = i.constraints

copy(n::Nothing) = nothing

copy(i::Individual) = Individual(copy(i.variables), copy(i.objectives), copy(i.constraints))