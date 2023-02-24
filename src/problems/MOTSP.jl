struct MOTSP <: AbstractProblem
    M 
    D
    C # Adjacency matrix of each map
    maxFE
    fn
    truepf
end

encoding(::MOTSP) = :permutation
parameter_properties(::MOTSP) = [:M, :D, :maxFE]

IGD(_, ::MOTSP) = NaN

# c: Correlation parameter
function MOTSP(; M = 2, D = 30, c = 0, maxFE = 10_000)
    path = "settings/MOTSP-M$(M)-D$(D)-c$(c).jld2"
    if isfile(path)
        data = load(path)
        C = data["C"]
    else
        C = Array{Any}(undef,M)
        C[1] = rand(D, D)
        for i = 2:M
            C[i] = c * C[i-1] + (1-c)*rand(D, D)
        end
        for i = 1:M
            C[i] = tril(C[i], -1) + triu(C[i]', 1)
        end
        save(path, "C", C)
    end

    function motsp(x)
        D = length(x)
        v = zeros(M)
        for i = 1:M
            for k = 1 : D-1
                v[i] = v[i] + C[i][x[k], x[k+1]]
            end
            v[i] = v[i] + C[i][D, 1]
        end
        v
    end

    truepf = zeros(M) .+ D

    return MOTSP(M, D, C, maxFE, motsp, truepf)
end
