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

abstract type ZDT <: TestProblem end

abstract type ZDT1 <: ZDT end

function ZDT1(m = 30, n_solutions = 100)
    f(x) = begin
        g = 1.0 + 9.0 * ( sum(x[2:end]) / (length(x)-1) )
        h = 1 - sqrt(x[1] / g)
        return [ x[1], g*h ]
    end

    bounds = Array([zeros(m) ones(m)]')

    x = range(0, 1, length=n_solutions)

    X = [vcat(x[i], zeros(m - 1)) for i in 1:n_solutions]
    pareto_set = [ Individual(x, f(x)) for x in X ]

    # @infiltrate

    return f, bounds, pareto_set
end

function truepf(x)
    "truepf"
end

thename(::Type{ZDT1}) = "ZDT1"

abstract type ZDT2 <: ZDT end

function ZDT2(D = 30, n_solutions = 100)
    f(x) = begin
        gx = 1.0 + 9.0 * ( sum(x[2:end]) / (length(x)-1) )
        return [x[1], gx*(1 - (x[1] / gx)^2) ]
    end
    bounds = Array([zeros(D) ones(D)]')

    x = range(0, 1, length=n_solutions)

    x = [vcat(x[i], zeros(D - 1)) for i in 1:n_solutions]
    pareto_set = [ Individual(xx, f(xx)) for xx in x ]

    return f, bounds, pareto_set
end

abstract type ZDT3 <: ZDT end

function ZDT3(D = 30, n_solutions = 100)
    f(x) = begin
        gx = 1.0 + 9.0 * ( sum(x[2:end]) / (length(x)-1) )
        a = x[1] / gx
        return [x[1], gx*(1 - sqrt(a) - a*sin(10π*x[1])) ]
    end
    bounds = Array([zeros(D) ones(D)]')

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
    pareto_set = [ Individual(xx, f(xx)) for xx in x ]

    return f, bounds, get_non_dominated_solutions(pareto_set)
end

abstract type ZDT4 <: ZDT end

function ZDT4(D = 10, n_solutions = 100)
    f(x) = begin
        gx = 1.0 + 10*(length(x)-1) + sum( x[2:end].^2 - 10cos.(4π*x[2:end]))
        return [x[1], gx*(1 - sqrt(x[1] / gx)) ]
    end
    bounds = Array([-5zeros(D) 5ones(D)]')
    bounds[:,1] = [0, 1.0]

    x = range(0, 1, length=n_solutions)

    x = [vcat(x[i], zeros(D - 1)) for i in 1:n_solutions]
    pareto_set = [ Individual(xx, f(xx)) for xx in x ]

    return f, bounds, pareto_set
end

abstract type ZDT6 <: ZDT end

function ZDT6(D = 10, n_solutions = 100)
    f(x) = begin
        gx = 1.0 + 9.0 * ( sum(x[2:end]) / (length(x)-1.0) )^(0.25)
        ff1 = 1.0 - exp(-4.0x[1])*sin(6.0π*x[1])^6
        return [ ff1 , gx*(1.0 - (ff1 / gx)^2) ]
    end

    bounds = Array([zeros(D) ones(D)]')

    #x = range(0, 1, length=n_solutions)

    #x = [vcat(x[i], zeros(D - 1)) for i in 1:n_solutions]
    xx = range(0.2807753191, 1, length=100)
    yy = 1 .- (xx).^2
    pareto_set = [ Individual(zeros(0), [xx[i], yy[i]]) for i in 1:length(xx) ]

    return f, bounds, pareto_set
end
