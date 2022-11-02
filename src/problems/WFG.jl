abstract type WFG <: AbstractProblem end

parameter_properties(::WFG) = [:M, :D, :K, :maxFE]

struct WFG1 <: WFG
    M 
    D
    K
    L
    maxFE
    fn
    constraints
    truepf
end

function WFG1(; M=3, D=M+9, K=M-1, maxFE = 10_000)
    N = D
    lower = zeros(N)
    upper = 2 : 2 : 2N
    L = D - K
    D = 1
    S = 2 : 2 : 2M
    A = ones(M-1)

    function wfg1(z)
        z₀₁ = z ./ upper

        t1 = zeros(N)
        t1[1:K] .= z₀₁[1:K]
        t1[K+1:end] .= s_linear.(z₀₁[K+1:end], 0.35)

        t2 = zeros(N)
        t2[1:K] .= t1[1:K]
        t2[K+1:end] .= b_flat.(t1[K+1:end],0.8,0.75,0.85)

        t3 = zeros(N)
        t3 .= b_poly.(t2, 0.02)

        t4 = zeros(M)
        for i = 1 : M-1
            t4[i] = r_sum(rangef(t3, (i-1)*K/(M-1)+1:i*K/(M-1)), 
                          2*((i-1)*K/(M-1)+1) : 2 : 2*i*K/(M-1))
        end
        t4[M] = r_sum(t3[K+1:K+L], 2*(K+1) : 2 : 2*(K+L))

        x = zeros(M)
        for i = 1 : M-1
            x[i] = max(t4[M],A[i])*(t4[i]-0.5)+0.5
        end
        x[M] = t4[M]

        h = convex(x)
        h[M] = mixed(x, 1, 5)
        fill(D*x[M],M) + S.*h
    end

    bounds = Array([lower upper]')
    constraints = BoxConstraints(bounds)

    truepf = wfg_truepf(1, wfg1, K, L, 10000)

    return WFG1(M, N, K, L, maxFE, wfg1, constraints, truepf)
end

struct WFG2 <: WFG
    M 
    D
    K
    L
    maxFE
    fn
    constraints
    truepf
end

function WFG2(; M=3, D=M+9, K=M-1, maxFE = 10_000)
    N = D
    lower = zeros(N)
    upper = 2 : 2 : 2N
    L = N - K
    D = 1
    S = 2 : 2 : 2M
    A = ones(M-1)

    function wfg2(z)
        z₀₁ = z ./ upper

        t1 = zeros(N)
        t1[1:K] .= z₀₁[1:K]
        t1[K+1:end] .= s_linear.(z₀₁[K+1:end], 0.35)

        t2 = zeros(N)
        t2[1:K] .= t1[1:K]
        for i = K+1:Int(K+L/2)
            t2[i] = r_nonsep(t1[Int(K+2*(i-K)-1):Int(K+2*(i-K))], 2)
        end
        
        t3 = zeros(M)
        for i = 1:M-1
            t3[i] = r_sum(rangef(t2, (i-1)*K/(M-1)+1:i*K/(M-1)), ones(Int(K/(M-1))))
        end
        t3[M] = r_sum(rangef(t2, K+1:K+L/2), ones(Int(L/2)))

        x = zeros(M)
        for i = 1 : M-1
            x[i] = max(t3[M],A[i])*(t3[i]-0.5)+0.5
        end
        x[M] = t3[M]

        h = convex(x)
        h[M] = disc(x, 1, 1, 5)
        fill(D*x[M],M) + S.*h
    end

    bounds = Array([lower upper]')
    constraints = BoxConstraints(bounds)

    truepf = wfg_truepf(2, wfg2, K, L, 10000)

    return WFG2(M, N, K, L, maxFE, wfg2, constraints, truepf)
end

struct WFG3 <: WFG
    M 
    D
    K
    L
    maxFE
    fn
    constraints
    truepf
end

function WFG3(; M=3, D=M+9, K=M-1, maxFE = 10_000)
    N = D
    lower = zeros(N)
    upper = 2 : 2 : 2N
    L = D - K
    D = 1
    S = 2 : 2 : 2M
    A = [1; zeros(M-2)]

    function wfg3(z)
        z₀₁ = z ./ upper

        t1 = zeros(N)
        t1[1:K] .= z₀₁[1:K]
        t1[K+1:end] .= s_linear.(z₀₁[K+1:end], 0.35)

        t2 = zeros(N)
        t2[1:K] .= t1[1:K]
        for i = K+1:Int(K+L/2)
            t2[i] = r_nonsep(t1[Int(K+2*(i-K)-1):Int(K+2*(i-K))], 2)
        end
        
        t3 = zeros(M)
        for i = 1:M-1
            t3[i] = r_sum(rangef(t2, (i-1)*K/(M-1)+1:i*K/(M-1)), ones(Int(K/(M-1))))
        end
        t3[M] = r_sum(rangef(t2, K+1:K+L/2), ones(Int(L/2)))

        x = zeros(M)
        for i = 1 : M-1
            x[i] = max(t3[M],A[i])*(t3[i]-0.5)+0.5
        end
        x[M] = t3[M]

        h = linear(x)
        fill(D*x[M],M) + S.*h
    end

    bounds = Array([lower upper]')
    constraints = BoxConstraints(bounds)

    truepf = wfg_truepf(3, wfg3, K, L, 10000)

    return WFG3(M, N, K, L, maxFE, wfg3, constraints, truepf)
end

struct WFG4 <: WFG
    M 
    D
    K
    L
    maxFE
    fn
    constraints
    truepf
end

function WFG4(; M=3, D=M+9, K=M-1, maxFE = 10_000)
    N = D
    lower = zeros(N)
    upper = 2 : 2 : 2N
    L = D - K
    D = 1
    S = 2 : 2 : 2M
    A = ones(M-1)

    function wfg4(z)
        z₀₁ = z ./ upper

        t1 = s_multi.(z₀₁, 30, 10, 0.35)

        t2 = zeros(M)
        for i = 1:M-1
            t2[i] = r_sum(rangef(t1, (i-1)*K/(M-1)+1:i*K/(M-1)), ones(Int(K/(M-1))))
        end
        t2[M] = r_sum(t1[K+1:N], ones(L))

        x = zeros(M)
        for i = 1 : M-1
            x[i] = max(t2[M],A[i])*(t2[i]-0.5)+0.5
        end
        x[M] = t2[M]

        h = concave(x)
        fill(D*x[M],M) + S.*h
    end

    bounds = Array([lower upper]')
    constraints = BoxConstraints(bounds)

    truepf = wfg_truepf(4, wfg4, K, L, 10000)

    return WFG4(M, N, K, L, maxFE, wfg4, constraints, truepf)
end

struct WFG5 <: WFG
    M 
    D
    K
    L
    maxFE
    fn
    constraints
    truepf
end

function WFG5(; M=3, D=M+9, K=M-1, maxFE = 10_000)
    N = D
    lower = zeros(N)
    upper = 2 : 2 : 2N
    L = D - K
    D = 1
    S = 2 : 2 : 2M
    A = ones(M-1)

    function wfg5(z)
        z₀₁ = z ./ upper

        t1 = s_decept.(z₀₁, 0.35, 0.001, 0.05)

        t2 = zeros(M)
        for i = 1:M-1
            t2[i] = r_sum(rangef(t1, (i-1)*K/(M-1)+1:i*K/(M-1)), ones(Int(K/(M-1))))
        end
        t2[M] = r_sum(t1[K+1:N], ones(L))

        x = zeros(M)
        for i = 1 : M-1
            x[i] = max(t2[M],A[i])*(t2[i]-0.5)+0.5
        end
        x[M] = t2[M]

        h = concave(x)
        fill(D*x[M],M) + S.*h
    end

    bounds = Array([lower upper]')
    constraints = BoxConstraints(bounds)

    truepf = wfg_truepf(5, wfg5, K, L, 10000)

    return WFG5(M, N, K, L, maxFE, wfg5, constraints, truepf)
end

struct WFG6 <: WFG
    M 
    D
    K
    L
    maxFE
    fn
    constraints
    truepf
end

function WFG6(; M=3, D=M+9, K=M-1, maxFE = 10_000)
    N = D
    lower = zeros(N)
    upper = 2 : 2 : 2N
    L = D - K
    D = 1
    S = 2 : 2 : 2M
    A = ones(M-1)

    function wfg6(z)
        z₀₁ = z ./ upper

        t1 = zeros(N)
        t1[1:K] .= z₀₁[1:K]
        t1[K+1:end] .= s_linear.(z₀₁[K+1:end], 0.35)

        t2 = zeros(M)
        for i = 1:M-1
            t2[i] = r_nonsep(rangef(t1, (i-1)*K/(M-1)+1:i*K/(M-1)), K/(M-1))
        end
        t2[M] = r_nonsep(t1[K+1:N], L)

        x = zeros(M)
        for i = 1 : M-1
            x[i] = max(t2[M],A[i])*(t2[i]-0.5)+0.5
        end
        x[M] = t2[M]

        h = concave(x)
        fill(D*x[M],M) + S.*h
    end

    bounds = Array([lower upper]')
    constraints = BoxConstraints(bounds)

    truepf = wfg_truepf(6, wfg6, K, L, 10000)

    return WFG6(M, N, K, L, maxFE, wfg6, constraints, truepf)
end

struct WFG7 <: WFG
    M 
    D
    K
    L
    maxFE
    fn
    constraints
    truepf
end

function WFG7(; M=3, D=M+9, K=M-1, maxFE = 10_000)
    N = D
    lower = zeros(N)
    upper = 2 : 2 : 2N
    L = D - K
    D = 1
    S = 2 : 2 : 2M
    A = ones(M-1)

    function wfg7(z)
        z₀₁ = z ./ upper

        t1 = zeros(N)
        for i = 1:K
            t1[i] = b_param(z₀₁[i], r_sum(z₀₁[i+1:N], ones(K+L-i)), 0.98/49.98, 0.02, 50)
        end
        t1[K+1:N] = z₀₁[K+1:N]

        t2 = zeros(N)
        t2[1:K] .= t1[1:K]
        t2[K+1:end] .= s_linear.(t1[K+1:end], 0.35)

        t3 = zeros(M)
        for i = 1:M-1
            t3[i] = r_sum(rangef(t2, (i-1)*K/(M-1)+1:i*K/(M-1)), ones(Int(K/(M-1))))
        end
        t3[M] = r_sum(t2[K+1:N], ones(L))

        x = zeros(M)
        for i = 1 : M-1
            x[i] = max(t3[M],A[i])*(t3[i]-0.5)+0.5
        end
        x[M] = t3[M]

        h = concave(x)
        fill(D*x[M],M) + S.*h
    end

    bounds = Array([lower upper]')
    constraints = BoxConstraints(bounds)

    truepf = wfg_truepf(7, wfg7, K, L, 10000)

    return WFG7(M, N, K, L, maxFE, wfg7, constraints, truepf)
end

struct WFG8 <: WFG
    M 
    D
    K
    L
    maxFE
    fn
    constraints
    truepf
end

function WFG8(; M=3, D=M+9, K=M-1, maxFE = 10_000)
    N = D
    lower = zeros(N)
    upper = 2 : 2 : 2N
    L = D - K
    D = 1
    S = 2 : 2 : 2M
    A = ones(M-1)

    function wfg8(z)
        z₀₁ = z ./ upper

        t1 = zeros(N)
        t1[1:K] .= z₀₁[1:K]
        for i = K+1:N
            t1[i] = b_param(z₀₁[i], r_sum(z₀₁[1:i-1], ones(i-1)), 0.98/49.98, 0.02, 50)
        end

        t2 = zeros(N)
        t2[1:K] .= t1[1:K]
        t2[K+1:end] .= s_linear.(t1[K+1:end], 0.35)

        t3 = zeros(M)
        for i = 1:M-1
            t3[i] = r_sum(rangef(t2, (i-1)*K/(M-1)+1:i*K/(M-1)), ones(Int(K/(M-1))))
        end
        t3[M] = r_sum(t2[K+1:N], ones(L))

        x = zeros(M)
        for i = 1 : M-1
            x[i] = max(t3[M],A[i])*(t3[i]-0.5)+0.5
        end
        x[M] = t3[M]

        h = concave(x)
        fill(D*x[M],M) + S.*h
    end

    bounds = Array([lower upper]')
    constraints = BoxConstraints(bounds)

    truepf = wfg_truepf(8, wfg8, K, L, 10000)

    return WFG8(M, N, K, L, maxFE, wfg8, constraints, truepf)
end

struct WFG9 <: WFG
    M 
    D
    K
    L
    maxFE
    fn
    constraints
    truepf
end

function WFG9(; M=3, D=M+9, K=M-1, maxFE = 10_000)
    N = D
    lower = zeros(D)
    upper = 2 : 2 : 2D
    L = D - K
    D = 1
    S = 2 : 2 : 2M
    A = ones(M-1)

    function wfg9(z)
        z₀₁ = z ./ upper

        t1 = zeros(N)
        for i = 1:N-1
            t1[i] = b_param(z₀₁[i], r_sum(z₀₁[i+1:N], ones(K+L-i)), 0.98/49.98, 0.02, 50)
        end
        t1[N] = z₀₁[N]

        t2 = zeros(N)
        t2[1:K] .= s_decept.(t1[1:K], 0.35, 0.001, 0.05)
        t2[K+1:N] .= s_multi.(t1[K+1:N], 30, 95, 0.35)

        t3 = zeros(M)
        for i = 1:M-1
            t3[i] = r_nonsep(rangef(t2, (i-1)*K/(M-1)+1:i*K/(M-1)), K/(M-1))
        end
        t3[M] = r_nonsep(t2[K+1:N], L)

        x = zeros(M)
        for i = 1 : M-1
            x[i] = max(t3[M],A[i])*(t3[i]-0.5)+0.5
        end
        x[M] = t3[M]

        h = concave(x)
        fill(D*x[M],M) + S.*h
    end

    bounds = Array([lower upper]')
    constraints = BoxConstraints(bounds)

    truepf = wfg_truepf(9, wfg9, K, L, 10000)

    return WFG9(M, N, K, L, maxFE, wfg9, constraints, truepf)
end

function linear(x)
    M = length(x)
    v = zeros(M)
    v[1] = prod(x[1:M-1])
    for m = 2:M-1
        v[m] = prod(x[1:M-m]) * (1 - x[M-m+1])
    end
    v[M] = 1 - x[1]
    v
end

function convex(x)
    M = length(x)

    v = zeros(M)
    v[1] = prod(1 .- cos.(x[1:M-1]*π/2))
    for m = 2:M-1
        v[m] = prod(1 .- cos.(x[1:M-m]*π/2)) * (1 .- sin.(x[M-m+1]*π/2))
    end
    v[M] = 1 - sin(x[1]*π/2)
    v
end

function concave(x)
    M = length(x)
    v = zeros(M)
    v[1] = prod(sin.(x[1:M-1]*π/2))
    for m = 2:M-1
        v[m] = prod(sin.(x[1:M-m]*π/2)) * cos.(x[M-m+1]*π/2)
    end
    v[M] = cos(x[1]*π/2)
    v
end

function mixed(x, α, A)
    (1 - x[1] - cos(2A*π*x[1] + π/2) / (2A*π)) ^ α
end

function disc(x, α, β, A)
    1 - x[1]^α * cos(A * x[1]^β * π) ^ 2
end

function b_poly(y, α)
    y ^ α
end

function b_flat(y, A, B, C)
    v = A + min(0, floor(y-B)) * A * (B-y) / B - min(0, floor(C-y)) * (1-A) * (y-C) / (1-C)
    round(v*1e4)/1e4
end

function b_param(y, Y, A, B, C)
    y^(B+(C-B)*(A-(1-2*Y)*abs(floor(0.5-Y)+A)))
end

function s_linear(y, A)
    abs(y-A) / abs(floor(A-y)+A)
end

function s_decept(y, A, B, C)
    1+(abs(y-A)-B)*(floor(y-A+B)*(1-C+(A-B)/B)/(A-B)+floor(A+B-y)*(1-C+(1-A-B)/B)/(1-A-B)+1/B)
end

function s_multi(y, A, B, C)
    (1+cos((4*A+2)*pi*(0.5-abs(y-C)/2/(floor(C-y)+C)))+4*B*(abs(y-C)/2/(floor(C-y)+C))^2)/(B+2)
end

function r_sum(y, w)
    sum(y.*w) / sum(w)
end

function r_nonsep(y, A)
    v = 0
    for j = 1:length(y)
        tmp = 0
        for k = 0:A-2
            tmp += abs(y[j]-y[1+mod(j+k,length(y))])
        end
        v += y[j] + tmp
    end
    v / (length(y)/A) / ceil(A/2) / (1+2A-2*ceil(A/2))
end

function rangef(v, r)
    n = length(r)
    v[Int(round(r[1])):Int(round(r[1]))+n-1]
end

function wfg1_random_z(k, l)
    n = k+l
    z = zeros(n)
    for i = 1:k
        z[i] = rand() ^ 50.0
    end
    for i = k+1:n
        z[i] = 0.35
    end
    for i = 1:n
        z[i] *= 2.0i
    end
    z
end

function wfg2_to_7_random_z(k, l)
    n = k+l
    z = zeros(n)
    for i = 1:k
        z[i] = rand()
    end
    for i = k+1:n
        z[i] = 0.35
    end
    for i = 1:n
        z[i] *= 2.0i
    end
    z
end

function wfg8_random_z(k, l)
    n = k+l
    z = zeros(n)
    for i = 1:k
        z[i] = rand()
    end
    for i = k+1:n
        w = ones(i-1)
        u = r_sum(z[1:i-1], w)
        tmp1 = abs(floor(0.5-u) + 0.98/49.98)
        tmp2 = 0.02 + 49.98 * (0.98 / 49.98 - (1.0 - 2.0 * u) * tmp1)
        z[i] = 0.35 ^ (tmp2 ^ -1.0)
    end
    for i = 1:n
        z[i] *= 2.0i
    end
    z
end

function wfg9_random_z(k, l)
    n = k+l
    z = zeros(n)
    for i = 1:k
        z[i] = rand()
    end
    z[n] = 0.35
    for i = n-1:-1:k+1
        u = r_sum(z[i+1:n], ones(n-i))
        z[i] = 0.35 ^ ((0.02 + 1.96u) ^ -1)
    end
    for i = 1:n
        z[i] *= 2.0i
    end
    z
end

function wfg_truepf(which, wfg_fn, K, L, n)
    random_z_fn = if which == 1
        wfg1_random_z
    elseif 2 <= which <= 7
        wfg2_to_7_random_z
    elseif which == 8
        wfg8_random_z
    elseif which == 9
        wfg9_random_z
    else
        error("No WFG$which problem")
    end
    individuals = Vector{Individual}()
    for i = 1:n
        z = random_z_fn(K, L)
        obj = wfg_fn(z)
        individual = Individual(z, obj)
        push!(individuals, individual)
    end
    individuals
end

function WFG(; params...)
    [
        WFG1(; params...),
        WFG2(; params...),
        WFG3(; params...),
        WFG4(; params...),
        WFG5(; params...),
        WFG6(; params...),
        WFG7(; params...),
        WFG8(; params...),
        WFG9(; params...)
    ]
end