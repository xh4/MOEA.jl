abstract type AbstractProblem end

show(io::IO, problem::AbstractProblem) = print(io, identifier(problem))

function identifier(p::AbstractProblem)
    join([string(typeof(p).name.name), "(", properties_string(p, parameter_properties(p)), ")"])
end

function lower(problem::AbstractProblem)
    reshape(problem.constraints.bounds.bx, 2, problem.D)[1,:]
end

function upper(problem::AbstractProblem)
    reshape(problem.constraints.bounds.bx, 2, problem.D)[2,:]
end