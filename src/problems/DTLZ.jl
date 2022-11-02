abstract type DTLZ <: AbstractProblem end

parameter_properties(::DTLZ) = [:M, :D, :maxFE]

generic_sphere(ref_dirs) = ref_dirs ./ ([norm(v) for v in eachrow(ref_dirs)])

function DTLZ_hypersphere!(fx, x; α = 1)
    m = length(fx)
    for i = 1:m
        if i < m
            fx[i] *= prod(cos.((π*0.5) * x[1:m-i] .^ α ))
        end
        if i > 1
            fx[i] *= sin((π*0.5) * x[1+m - i] .^ α)
        end
    end
    fx
end

function DTLZ_hyperplane!(fx, x)
    m = length(fx)
    for i in 1:m
        fx[i] *= prod(x[1:m-i])
        if i > 1
            fx[i] *= 1 - x[1+m - i]
        end
    end
    fx
end

function DTLZ_g1(x, m)
    y = view(x, m:length(x)) .- 0.5
    100*( length(y) + sum( y.^2 - cos.(20π*y) ))
end

function DTLZ_g2(x, m)
    sum( (view(x, m:length(x)) .- 0.5).^2 )
end


"""
        DTLZ1(M = 3, ref_dirs = gen_ref_dirs(m, 12))

    DTLZ1 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
    where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
    with `n_solutions`.

    ### Parameters
    - `m` number of objective functions
    - `ref_dirs` number of Pareto solutions (default: Das and Dennis' method).

    Main properties:
    - convex
    - multifrontal
    """
struct DTLZ1 <: DTLZ
    M 
    D
    maxFE
    fn
    constraints
    truepf
end

function DTLZ1(; M = 3, D = M + 4, maxFE = 10_000, n = 10_000)
    function dtlz1(x)
        g = DTLZ_g1(x, M)
        fx = fill(0.5*(1 + g), M)
        DTLZ_hyperplane!(fx, x)
    end

    bounds = Array([zeros(D) ones(D)]')
    constraints = BoxConstraints(bounds)

    pf, _ = NBI(n, M)
    pf = 0.5pf
    truepf = [ Individual(zeros(D), fx) for fx in eachrow(pf) ]

    return DTLZ1(M, D, maxFE, dtlz1, constraints, truepf)
end

"""
        DTLZ2(M = 3, ref_dirs = gen_ref_dirs(m, 12))

    DTLZ2 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
    where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
    with `n_solutions`.

    ### Parameters
    - `m` number of objective functions
    - `ref_dirs` number of Pareto solutions (default: Das and Dennis' method).

    Main properties:
    - nonconvex
    - unifrontal
    """
struct DTLZ2 <: DTLZ
    M 
    D
    maxFE
    fn
    constraints
    truepf
end

function DTLZ2(; M = 3, D = M + 9, maxFE = 10_000, n = 10_000)
    function dtlz2(x)
        g = DTLZ_g2(x, M)
        fx = fill(1.0 + g, M)
        DTLZ_hypersphere!(fx, x)
    end

    bounds = Array([zeros(D) ones(D)]')
    constraints = BoxConstraints(bounds)

    ref_dirs, _ = NBI(n, M)
    pf = generic_sphere(ref_dirs)
    truepf = [ Individual(zeros(D), fx) for fx in eachrow(pf) ]

    return DTLZ2(M, D, maxFE, dtlz2, constraints, truepf)
end

"""
        DTLZ3(M = 3, ref_dirs = gen_ref_dirs(m, 12))

    DTLZ3 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
    where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
    with `n_solutions`.

    ### Parameters
    - `m` number of objective functions
    - `ref_dirs` number of Pareto solutions (default: Das and Dennis' method).

    Main properties:
    - nonconvex
    - multifrontal
    """
struct DTLZ3 <: DTLZ
    M 
    D
    maxFE
    fn
    constraints
    truepf
end

function DTLZ3(; M = 3, D = M + 9, maxFE = 10_000, n = 10_000)
    function dtlz3(x)
        g = DTLZ_g1(x, M)
        fx = fill(1 + g, M)
        DTLZ_hypersphere!(fx, x)
    end

    bounds = Array([zeros(D) ones(D)]')
    constraints = BoxConstraints(bounds)

    ref_dirs, _ = NBI(n, M)
    pf = generic_sphere(ref_dirs)
    truepf = [ Individual(zeros(D), fx) for fx in eachrow(pf) ]

    return DTLZ3(M, D, maxFE, dtlz3, constraints, truepf)
end


"""
        DTLZ4(M = 3, ref_dirs = gen_ref_dirs(m, 12))

    DTLZ4 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
    where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
    with `n_solutions`.

    ### Parameters
    - `m` number of objective functions
    - `ref_dirs` number of Pareto solutions (default: Das and Dennis' method).

    Main properties:
    - nonconvex
    - unifrontal
    """
struct DTLZ4 <: DTLZ
    M 
    D
    maxFE
    fn
    constraints
    truepf
end

function DTLZ4(; M = 3, D = M + 9, maxFE = 10_000, n = 10_000)
    function dtlz4(x)
        g = DTLZ_g2(x, M)
        fx = fill(1.0 + g, M)
        DTLZ_hypersphere!(fx, x; α = 100)
    end

    bounds = Array([zeros(D) ones(D)]')
    constraints = BoxConstraints(bounds)

    ref_dirs, _ = NBI(n, M)
    pf = generic_sphere(ref_dirs)
    truepf = [ Individual(zeros(D), fx) for fx in eachrow(pf) ]

    return DTLZ4(M, D, maxFE, dtlz4, constraints, truepf)
end

"""
        DTLZ5(M = 3)

    DTLZ5 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
    where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
    with `n_solutions`.

    ### Parameters
    - `m` number of objective functions

    """
struct DTLZ5 <: DTLZ
    M 
    D
    maxFE
    fn
    constraints
    truepf
end

function DTLZ5(; M = 3, D = M + 9, maxFE = 10_000, n = 10_000)
    function dtlz5(x)
        g = DTLZ_g2(x, M)
        fθ = fill(1.0 + g, M)
        θ = @. 1 / (2*(1 + g)) * (1 + 2g * x[1:M-1] )
        θ[1] = x[1]
    
        DTLZ_hypersphere!(fθ, θ)
    end

    bounds = Array([zeros(D) ones(D)]')
    constraints = BoxConstraints(bounds)

    if M == 3
        ref_dirs, _ = NBI(n, 2)
        X = fill(0.5, size(ref_dirs)[1], D)

        for i in 1:size(ref_dirs)[1]
            X[i,1:2] = ref_dirs[i,:]
        end
        truepf = [ Individual(X[i,:], dtlz5(X[i,:])) for i in 1:size(ref_dirs)[1] ]
    else
        truepf = []
    end
    return DTLZ5(M, D, maxFE, dtlz5, constraints, truepf)
end


"""
        DTLZ6(M = 3)

    DTLZ6 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
    where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
    with `n_solutions`.

    ### Parameters
    - `m` number of objective functions

    """
struct DTLZ6 <: DTLZ
    M 
    D
    maxFE
    fn
    constraints
    truepf
end

function DTLZ6(; M = 3, D = M + 9, maxFE = 10_000, n = 10_000)
    function dtlz6(x)
        g = sum(x[M:end] .^ 0.1)
        fθ = fill(1.0 + g, M)
        θ = @. 1 / (2*(1 + g)) * (1 + 2g * x[1:M-1] )
        θ[1] = x[1]
    
        DTLZ_hypersphere!(fθ, θ)   
    end

    bounds = Array([zeros(D) ones(D)]')
    constraints = BoxConstraints(bounds)

    if M == 3
        ref_dirs, _ = NBI(n, 2)
        X = fill(0.0, size(ref_dirs)[1], D)

        for i in 1:size(ref_dirs)[1]
            X[i,1:2] = ref_dirs[i,:]
        end
        truepf = [ Individual(X[i,:], dtlz6(X[i,:])) for i in 1:size(ref_dirs)[1] ]
    else
        truepf = []
    end

    return DTLZ6(M, D, maxFE, dtlz6, constraints, truepf)
end

function DTLZ(; params...)
    [
        DTLZ1(; params...),
        DTLZ2(; params...),
        DTLZ3(; params...),
        DTLZ4(; params...),
        DTLZ5(; params...),
        DTLZ6(; params...)
    ]
end