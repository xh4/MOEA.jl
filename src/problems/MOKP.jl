struct MOKP <: AbstractProblem
    M 
    D
    P # Profit of each item according to each knapsack
    W # Weight of each item according to each knapsack
    maxFE
    fn
    truepf
end

encoding(::MOKP) = :binary
parameter_properties(::MOKP) = [:M, :D, :maxFE]

IGD(_, ::MOKP) = NaN

function MOKP(; M = 2, D = 250, maxFE = 10_000)
    path = "settings/MOKP-M$(M)-D$(D).jld2"
    if isfile(path)
        data = load(path)
        P = data["P"]
        W = data["W"]
    else
        P = Int.(round.(rand(M, D) .* 90)) .+ 10
        W = Int.(round.(rand(M, D) .* 90)) .+ 10
        save(path, "P", P, "W", W)
    end

    function mokp(x)        
        # Repair invalid solutions
        C = sum(W, dims=2) ./ 2
        rank = sortperm(maximum(P ./ W, dims=1))
        while any(W * x .> C)
            k = findfirst(x[rank])
            x[rank[k]] = 0
        end

        # Calculate objective values
        (sum(P, dims=2)' - x' * P')[1,:]
    end

    truepf = sum(P, dims=2)[:,1]

    return MOKP(M, D, P, W, maxFE, mokp, truepf)
end