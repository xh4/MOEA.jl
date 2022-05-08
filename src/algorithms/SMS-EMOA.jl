Base.@kwdef struct SMSEMOA <: AbstractOptimizer
    # number of the subproblems
    N::Int = 100
    # number of the weight vectors in the neighborhoor of each weight vector
    T::Int = ceil(N/10)
    crossover = BINX(0.5)
    mutation = PLM()
    metrics::ConvergenceMetrics = ConvergenceMetric[GD(), GD(true)]
end

population_size(method::SMSEMOA) = method.N
default_options(method::SMSEMOA) = (iterations = 1000,)
summary(m::SMSEMOA) =
    "SMS_EMOA[P=$(population_size(m)),]"
show(io::IO, m::SMSEMOA) = print(io, summary(m))

mutable struct SMSEMOAState <: AbstractOptimizerState
    iteration
    start_time
    stop_time
    metrics::ConvergenceMetrics # collection of convergence metrics

    population                  # population
    m                           # number of objectives
    W                           # weight vectors
    B                           # neighbours of each solution
    z                           # the best value found so far
    EP                          # nondominated solutions
end
pfront(s::SMSEMOAState) = s.EP
value(s::SMSEMOAState) = objectives(pfront(s))
minimizer(s::SMSEMOAState) = objectives(pfront(s))
copy(s::SMSEMOAState) = SMSEMOAState(s.iteration, s.start_time, s.stop_time, copy(s.metrics),
                                 copy(s.population),
                                 copy(s.m),
                                 copy(s.W),
                                 copy(s.B),
                                 copy(s.z),
                                         copy(s.EP))

function initial_state(method::SMSEMOA, objective, population, options)
    value!(objective, population)

    return nothing
end

function update_state!(
    method::SMSEMOA,
    state,
    objective,
    constraints,
    options
)

    return false
end
