
#=
####################################################################################
####################################################################################
####################################################################################
#                Constrained DTLZ
####################################################################################
####################################################################################
####################################################################################
=#

function C1_DTLZ1_f(x, m = 3)
    g = DTLZ_g1(x, m)
    D = length(x)

    fx = fill(0.5*(1 + g), m)
    DTLZ_hyperplane!(fx, x)


    c = fx[end]/0.6 + sum(fx[1:end-1]/0.5) - 1;
    return fx, [c], [0.0]
end


"""
        C1_DTLZ1(m = 3, ref_dirs = gen_ref_dirs(m, 12))

    C1_DTLZ1 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
    where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
    with `n_solutions`.

    ### Parameters
    - `m` number of objective functions
    - `ref_dirs` number of Pareto solutions (default: Das and Dennis' method).

    Main properties:
    - convex
    - multifrontal
    - constraints type 1
    """
function C1_DTLZ1(m=3, ref_dirs = gen_ref_dirs(m, 12))

    D = 10 + m - 1

    bounds = Array([zeros(D) ones(D)]')

    pf = 0.5ref_dirs
    pareto_set = [ generateChild(zeros(0), (fx, [0.0], [0.0])) for fx in pf ]

    return C1_DTLZ1_f, bounds, pareto_set
end


function C1_DTLZ3_f(x,m=3)
    g = DTLZ_g1(x, m)
    fx = fill(1 + g, m)
    DTLZ_hypersphere!(fx, x)

    if m == 2
        r = 6
    elseif m <= 3
        r = 9
    elseif m <= 8
        r = 12.5
    else
        r = 15
    end
    c = -(sum(fx .^ 2) - 16) * (sum(fx .^ 2) - r^2);

    return fx, [c], [0.0]
end

"""
        C1_DTLZ3(m = 3, ref_dirs = gen_ref_dirs(m, 12))

    DTLZ3 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
    where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
    with `n_solutions`.

    ### Parameters
    - `m` number of objective functions
    - `ref_dirs` number of Pareto solutions (default: Das and Dennis' method).

    Main properties:
    - nonconvex
    - multifrontal
    - constraints type 1
    """
function C1_DTLZ3(m = 3, ref_dirs = gen_ref_dirs(m, 12))

    D = 10 + m - 1

    bounds = Array([zeros(D) ones(D)]')

    pf = generic_sphere(ref_dirs)
    pareto_set = [ generateChild(zeros(0), (fx, [0.0], [0.0])) for fx in pf ]

    return C1_DTLZ3_f, bounds, pareto_set
end


function C2_DTLZ2_f(x, m = 3)
    g = DTLZ_g2(x, m)
    fx = fill(1.0 + g, m)
    DTLZ_hypersphere!(fx, x)

    if m == 3
        r = 0.4
    else
        r = 0.5
    end


    p1 = maximum( i -> (fx[i] .- 1).^2 + sum(fx[1:m].^2 .- fx[i]^2) - r^2,1:m)
    p2 = sum((fx .- 1 / sqrt(m)).^2) - r^2
    return fx, [-max(p1, p2)], [0.0]
end

"""
        C2_DTLZ2(m = 3, ref_dirs = gen_ref_dirs(m, 12))

    DTLZ2 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
    where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
    with `n_solutions`.

    ### Parameters
    - `m` number of objective functions
    - `ref_dirs` number of Pareto solutions (default: Das and Dennis' method).

    Main properties:
    - nonconvex
    - unifrontal
    - contraints type 2
    """
function C2_DTLZ2(m=3, ref_dirs = gen_ref_dirs(m, 12))

    D = 10 + m - 1

    bounds = Array([zeros(D) ones(D)]')

    pf = generic_sphere(ref_dirs)
    pareto_set = [ generateChild(zeros(0), (fx, [0.0], [0.0])) for fx in pf ]

    return C2_DTLZ2_f, bounds, pareto_set
end



function C3_DTLZ1_f(x, m = 3)
    g = DTLZ_g1(x, m)
    D = length(x)

    fx = fill(0.5*(1 + g), m)
    DTLZ_hyperplane!(fx, x)


    c = [sum(fx[j] .+ fx/0.5 .- fx[j]/0.5) - 1 for j in 1:m]
    return fx, -c, [0.0]
end

#=
"""
C3_DTLZ1(m = 3, ref_dirs = gen_ref_dirs(m, 12))

C3_DTLZ1 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
with `n_solutions`.

### Parameters
- `m` number of objective functions
- `ref_dirs` number of Pareto solutions (default: Das and Dennis' method).

Main properties:
- convex
- multifrontal
- constraints type 3
"""
function C3_DTLZ1(m=3, ref_dirs = gen_ref_dirs(m, 12))

D = m + 4

bounds = Array([zeros(D) ones(D)]')

@warn "C3_DTLZ1 is under development"

pf = ref_dirs./([(0.5r -6*maximum(r)) for r in ref_dirs])
pareto_set = [ generateChild(zeros(0), (fx, [0.0], [0.0])) for fx in pf ]

return C3_DTLZ1_f, bounds, pareto_set
end
=#




function C3_DTLZ4_f(x,m=3)
    g = DTLZ_g2(x, m)
    fx = fill(1.0 + g, m)
    DTLZ_hypersphere!(fx, x; α = 100)

    c = [fx[j]^2/4  + sum(fx.^2 .- fx[j]^2) - 1 for j in 1:m]
    return fx, -c, [0.0]
end

"""
        C3_DTLZ4(m = 3, ref_dirs = gen_ref_dirs(m, 12))

    C3_DTLZ4 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
    where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
    with `n_solutions`.

    ### Parameters
    - `m` number of objective functions
    - `ref_dirs` number of Pareto solutions (default: Das and Dennis' method).

    Main properties:
    - nonconvex
    - unifrontal
    - constraints type 3
    """
function C3_DTLZ4(m = 3, ref_dirs = gen_ref_dirs(m, 12))

    D = m + 4

    bounds = Array([zeros(D) ones(D)]')

    pf = ref_dirs./([sqrt.(sum(r.^2)-3/4*maximum(r.^2)) for r in ref_dirs])

    pareto_set = [ generateChild(zeros(0), (fx, [0.0], [0.0])) for fx in pf ]

    return C3_DTLZ4_f, bounds, pareto_set
end
