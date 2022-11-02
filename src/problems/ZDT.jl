
"""
@article{zitzler2000comparison,
  title={Comparison of multiobjective evolutionary algorithms: Empirical results},
  author={Zitzler, Eckart and Deb, Kalyanmoy and Thiele, Lothar},
  journal={Evolutionary computation},
  volume={8},
  number={2},
  pages={173--195},
  year={2000},
  publisher={MIT Press}
}

ZDT 测试函数的特点:
- 有固定的 2 个目标
- 有确定的 True Pareto Front
"""

abstract type ZDT <: AbstractProblem end

parameter_properties(::ZDT) = [:D, :maxFE]

struct ZDT1 <: ZDT
    M
    D
    maxFE
    fn
    constraints
    truepf
end

function ZDT1(; D = 30, maxFE = 10_000)
    function zdt1(x)
        g = 1.0 + 9.0 * ( sum(x[2:end]) / (length(x)-1) )
        h = 1 - sqrt(x[1] / g)
        [ x[1], g*h ]
    end

    bounds = Array([zeros(D) ones(D)]')
    constraints = BoxConstraints(bounds)

    n_solutions = 10_000
    x = range(0, 1, length=n_solutions)
    y = 1 .- x .^ 0.5
    truepf = [ Individual(zeros(D), [x[i], y[i]]) for i in 1:n_solutions ]

    return ZDT1(2, D, maxFE, zdt1, constraints, truepf)
end

struct ZDT2 <: ZDT
    M 
    D
    maxFE
    fn
    constraints
    truepf
end

function ZDT2(; D = 30, maxFE = 10_000)
    function zdt2(x)
        gx = 1.0 + 9.0 * ( sum(x[2:end]) / (length(x)-1) )
        [x[1], gx*(1 - (x[1] / gx)^2) ]
    end

    bounds = Array([zeros(D) ones(D)]')
    constraints = BoxConstraints(bounds)

    n_solutions = 10_000
    x = range(0, 1, length=n_solutions)
    x = [vcat(x[i], zeros(D - 1)) for i in 1:n_solutions]
    truepf = [ Individual(xx, zdt2(xx)) for xx in x ]

    return ZDT2(2, D, maxFE, zdt2, constraints, truepf)
end

struct ZDT3 <: ZDT
    M 
    D
    maxFE
    fn
    constraints
    truepf
end

function ZDT3(; D = 30, maxFE = 10_000)
    function zdt3(x)
        gx = 1.0 + 9.0 * ( sum(x[2:end]) / (length(x)-1) )
        a = x[1] / gx
        [x[1], gx*(1 - sqrt(a) - a*sin(10π*x[1])) ]
    end

    bounds = Array([zeros(D) ones(D)]')
    constraints = BoxConstraints(bounds)

    n_solutions = 10_000
    if n_solutions < 6
        n_solutions = 6
    end
    regions = [ 0 0.0830015349
                0.182228780 0.2577623634
                0.4093136748 0.4538821041
                0.6183967944 0.6525117038
                0.8233317983 0.8518328654]
    n = Int(ceil(n_solutions / size(regions, 1)))
    x = Float64[]
    for i in 1:size(regions, 1)
        x = vcat(x, range(regions[i,1], regions[i,2], length=n))
    end
    x = [vcat(x[i], zeros(D - 1)) for i in 1:n_solutions]
    truepf = get_non_dominated_solutions([ Individual(xx, zdt3(xx)) for xx in x ])

    return ZDT3(2, D, maxFE, zdt3, constraints, truepf)
end

struct ZDT4 <: ZDT
    M 
    D
    maxFE
    fn
    constraints
    truepf
end

function ZDT4(; D = 10, maxFE = 10_000)
    function zdt4(x)
        gx = 1.0 + 10*(length(x)-1) + sum( x[2:end].^2 - 10cos.(4π*x[2:end]))
        [x[1], gx*(1 - sqrt(x[1] / gx)) ]
    end

    bounds = Array([-5ones(D) 5ones(D)]')
    bounds[:,1] = [0, 1.0]
    constraints = BoxConstraints(bounds)

    n_solutions = 10_000
    x = range(0, 1, length=n_solutions)
    x = [vcat(x[i], zeros(D - 1)) for i in 1:n_solutions]
    truepf = [ Individual(xx, zdt4(xx)) for xx in x ]

    return ZDT4(2, D, maxFE, zdt4, constraints, truepf)
end

struct ZDT6 <: ZDT
    M 
    D
    maxFE
    fn
    constraints
    truepf
end

function ZDT6(; D = 10, maxFE = 10_000)
    function zdt6(x)
        gx = 1.0 + 9.0 * ( sum(x[2:end]) / (length(x)-1.0) )^(0.25)
        ff1 = 1.0 - exp(-4.0x[1])*sin(6.0π*x[1])^6
        [ ff1 , gx*(1.0 - (ff1 / gx)^2) ]
    end

    bounds = Array([zeros(D) ones(D)]')
    constraints = BoxConstraints(bounds)

    n_solutions = 10_000
    #x = range(0, 1, length=n_solutions)
    #x = [vcat(x[i], zeros(D - 1)) for i in 1:n_solutions]
    xx = range(0.2807753191, 1, length=n_solutions)
    yy = 1 .- (xx).^2
    truepf = [ Individual(zeros(0), [xx[i], yy[i]]) for i in 1:length(xx) ]

    return ZDT6(2, D, maxFE, zdt6, constraints, truepf)
end

function ZDT(; params...)
    [
        ZDT1(; params...),
        ZDT2(; params...),
        ZDT3(; params...),
        ZDT4(; params...),
        ZDT6(; params...)
    ]
end