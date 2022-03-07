module MOEA

using Random, LinearAlgebra, Statistics
using Base: @kwdef
using TOML
using JLD
using UnPack: @unpack
using StackViews
using Plots
using Distributions
using NLSolversBase:
    NLSolversBase,
    AbstractObjective,
    ConstraintBounds,
    AbstractConstraints,
    nconstraints_x,
    nconstraints
using RedefStructs: @redef
using Infiltrator: @infiltrate
import NLSolversBase: f_calls, value, value!
import Base:
    show,
    copy,
    minimum,
    summary,
    identity,
    getproperty,
    rand,
    getindex,
    length,
    copyto!,
    setindex!,
    replace,
    print

export optimize,
    ### TYPE
    Individual,
    Population,
    Objective,
    variable,
    objective,
    ### Constraint
    BoxConstraints,
    ### TEST PROBLEM
    ZDT1,
    ZDT2,
    ZDT3,
    ZDT4,
    ZDT6,
    truepf,
    sphere,
    ackley,
    ### ALGORITHM
    GA,
    NSGA2,
    ### SELECTION
    ranklinear,
    uniformranking,
    roulette,
    rouletteinv,
    sus,
    susinv,
    tournament,
    truncation,
    uniformranking,
    ### CROSSOVER
    # binary
    SPX,  # single point
    TPX,  # two points
    SHFX, # shuffle
    UX,   # uniform
    BINX, # binary
    EXPX, # exponential
    BSX,  # binary subset
    # real value
    DC,   # discrete
    AX,   # average
    WAX,  # weighted average
    IC,   # intermediate
    LC,   # line
    HX,   # heuristic
    LX,   # Laplace
    MILX, # mixed integer Laplace
    SBX,  # simulated binary
    ### MUTATION
    ## GA
    # binary
    flip,
    inversion,
    # real value
    uniform,
    BGA,  # domain range
    PM,   # power
    MIPM, # mixed-integer power
    PLM,  # polynomial
    # combinatorial
    insertion,
    swap2,
    scramble,
    shifting

include("state.jl")
include("optimizer.jl")
include("individual.jl")
include("population.jl")
include("objective.jl")
include("options.jl")
include("termination.jl")
include("result.jl")
include("utilities.jl")
include("constraints.jl")
include("optimize.jl")

include("algorithms/ga.jl")
include("algorithms/nsga2.jl")

include("algorithms/selection.jl")
include("algorithms/crossover.jl")
include("algorithms/mutation.jl")
include("algorithms/moea.jl")

include("problems/type.jl")
include("problems/ZDT.jl")
include("problems/CEC2018.jl")
include("problems/Sphere.jl")
include("problems/Ackley.jl")

end
