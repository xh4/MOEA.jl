using MOEA
using Test
import Random: seed!
seed!(42)

for tests in [
    "DE.jl"
]
    include(tests)
end