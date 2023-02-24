module MOEA

using Random, LinearAlgebra, Statistics
using TOML
using JLD2
using UnPack: @unpack
using StackViews
using Plots
using Printf
using Distributions
using Combinatorics
using NLSolversBase:
    NLSolversBase,
    AbstractObjective,
    ConstraintBounds,
    AbstractConstraints,
    nconstraints_x,
    nconstraints
using Infiltrator: @infiltrate
using TimerOutputs: @timeit, TimerOutput
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
    Euclidean,
    CosineDist
# using Colors
# using CImGui
# using CImGui.ImGuiGLFWBackend
# using CImGui.ImGuiGLFWBackend.LibGLFW
# using CImGui.ImGuiGLFWBackend.LibCImGui
# using CImGui.ImGuiOpenGLBackend
# using CImGui.ImGuiOpenGLBackend.ModernGL
# using CImGui.CSyntax
# using CImGui.CSyntax.CStatic
# using GLFW
# using Images
# using TestImages
# using ColorTypes
using SQLite
using DBInterface
using DataFrames
using Serialization
using ConcurrentCollections
using ProgressMeter
using PrettyTables
using HypothesisTests
using Hyperscript: @tags, savehtml
import DefaultApplication
using StatsBase

export optimize,
    benchmark,
    experiment,
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
    GD,
    IGD,
    HV,
    FPL,
    ### TEST PROBLEM
    ZDT,
    ZDT1,
    ZDT2,
    ZDT3,
    ZDT4,
    ZDT6,
    DTLZ,
    DTLZ1,
    DTLZ2,
    DTLZ3,
    DTLZ4,
    DTLZ5,
    DTLZ6,
    WFG,
    WFG1,
    WFG2,
    WFG3,
    WFG4,
    WFG5,
    WFG6,
    WFG7,
    WFG8,
    WFG9,
    MOKP,
    MOTSP,
    truepf,
    sphere,
    ackley,
    ### ALGORITHM
    NSGAII,
    MOEAD,
    MOEADDE,
    MOEADDRA,
    MOEADAWA,
    MOEADURAW,
    MOEADUR,
    MOEADMY,
    SMSEMOA,
    SMSEMOAdp,
    RVEA,
    NSGAIII,
    ANSGAIII,
    ARMOEA,
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
    Flip,
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

include("optimizer.jl")
include("individual.jl")
include("population.jl")
include("objective.jl")
include("termination.jl")
include("result.jl")
include("constraints.jl")
include("parameter.jl")
include("problem.jl")
include("algorithm.jl")
include("optimize.jl")
include("benchmark.jl")
include("experiment.jl")

include("indicators/GD.jl")
include("indicators/IGD.jl")
include("indicators/hypervolume/FPL.jl")
include("indicators/hypervolume/MonteCarlo.jl")
include("indicators/HV.jl")
include("indicators/Spread.jl")

include("algorithms/Selection.jl")
include("algorithms/Crossover.jl")
include("algorithms/Mutation.jl")
include("algorithms/Evolution.jl")
include("algorithms/UniformPoints.jl")
include("algorithms/GA.jl")
include("algorithms/DE.jl")
include("algorithms/MOEA.jl")
include("algorithms/NSGA-II.jl")
include("algorithms/NSGA-III.jl")
include("algorithms/A-NSGA-III.jl")
include("algorithms/MOEA-D.jl")
include("algorithms/MOEA-D-DE.jl")
include("algorithms/MOEA-D-DRA.jl")
include("algorithms/MOEA-D-AWA.jl")
include("algorithms/MOEA-D-URAW.jl")
include("algorithms/MOEA-D-UR.jl")
include("algorithms/MOEA-D-MY.jl")
include("algorithms/MOEA-D-MY2.jl")
include("algorithms/MOEA-D-MY3.jl")
include("algorithms/SMS-EMOA.jl")
include("algorithms/SMS-EMOA-dp.jl")
include("algorithms/RVEA.jl")
include("algorithms/AR-MOEA.jl")

include("problems/ZDT.jl")
include("problems/DTLZ.jl")
include("problems/WFG.jl")
include("problems/CEC2018.jl")
include("problems/Sphere.jl")
include("problems/Ackley.jl")
include("problems/MOKP.jl")
include("problems/MOTSP.jl")

include("database.jl")
include("utilities.jl")

# include("gui/utils.jl")
# include("gui/image.jl")
# include("gui/parameter.jl")
# include("gui/problem.jl")
# include("gui/algorithm.jl")
# include("gui/indicator.jl")
# include("gui/playground.jl")
# include("gui/experiment.jl")
# include("gui/gui.jl")

end
