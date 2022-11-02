using MOEA
using Test
import Random: seed!
seed!(42)

for tests in [
    "DTLZ.jl"
]
    include(tests)
end