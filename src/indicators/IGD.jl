"""
    IGD(S,R)

Calculate an inverted generational distance, [`gd`](@ref), between set `S` and the reference set `R`.
Parameters are column-major matrices.
"""
IGD(S, R) = GD(R, S)

IGD(S, R::AbstractProblem) = IGD(S, R.truepf)