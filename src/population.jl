const Population = Vector{<:AbstractIndividual}

variables(pop::Population) = Matrix(mapreduce(variables, hcat, pop)')

objectives(pop::Population) = mapreduce(objectives, hcat, pop)'

constraints(pop::Population) = map(constraints, pop)

copy(pop::Population) = map(copy, pop)

function initial_population(algorithm, problem)
    if encoding(problem) == :real
        [Individual(rand(problem.D)) for i in 1:algorithm.N]
    elseif encoding(problem) == :binary
        [Individual(Bool.(round.(rand(problem.D)))) for i in 1:algorithm.N]
    elseif encoding(problem) == :permutation
        [Individual(shuffle(collect(1:problem.D))) for i in 1:algorithm.N]
    else 
        error("unknown problem encoding $(encoding(problem))")
    end
end