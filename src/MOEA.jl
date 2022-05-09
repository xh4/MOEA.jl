module MOEA

using Random, LinearAlgebra, Statistics
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
using Infiltrator: @infiltrate
using Colors
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
    print,
    sortperm,
    sortperm!
import Distances:
    pairwise,
    Euclidean
using CImGui
using CImGui.GLFWBackend
using CImGui.OpenGLBackend
using CImGui.GLFWBackend.GLFW
using CImGui.OpenGLBackend.ModernGL
using CImGui.CSyntax
using CImGui.CSyntax.CStatic
using ImGuiGLFWBackend
using ImGuiOpenGLBackend
using LibCImGui
using GLFW
using Images
using TestImages
using ColorTypes

export optimize,
    ### TYPE
    Individual,
    Population,
    Objective,
    variables,
    objectives,
    constraints,
    ### CONSTRAINT
    BoxConstraints,
    ### METRIC
    gd,
    igd,
    FPL,
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
    DE,
    NSGAII,
    MOEAD,
    MOEADDE,
    SMSEMOA,
    pfront,
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
include("constraints.jl")
include("optimize.jl")
include("FPL.jl")

include("algorithms/selection.jl")
include("algorithms/crossover.jl")
include("algorithms/mutation.jl")
include("algorithms/moea.jl")

include("algorithms/GA.jl")
include("algorithms/DE.jl")
include("algorithms/NSGA-II.jl")
include("algorithms/MOEA-D.jl")
include("algorithms/MOEA-D-DE.jl")
include("algorithms/SMS-EMOA.jl")

include("problems/type.jl")
include("problems/ZDT.jl")
include("problems/DTLZ.jl")
include("problems/CEC2018.jl")
include("problems/sphere.jl")
include("problems/ackley.jl")

include("utilities.jl")

include("gui/utils.jl")
include("gui/image.jl")
include("gui/parameter.jl")
include("gui/problem.jl")
include("gui/algorithm.jl")
include("gui/playground.jl")
include("gui/experiment.jl")
include("gui/gui.jl")

end
