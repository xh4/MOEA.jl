const Population = Vector{<:AbstractIndividual}

variables(pop::Population) = Matrix(mapreduce(variables, hcat, pop)')

objectives(pop::Population) = mapreduce(objectives, hcat, pop)'

constraints(pop::Population) = map(constraints, pop)

copy(pop::Population) = map(copy, pop)

function initial_population(method, objective)
end