function DasDennis(N, M)
    # gen_ref_dirs 会生成 N+1 个...，需修复
    W = reduce(vcat, gen_ref_dirs(M, N-1)')
    # gen_ref_dirs(M, N-1)
    W = broadcast((v) -> max(v, 1e-6), W)
    W, size(W)[1]
end

function TwoLayer(N, M)
    H1 = 1
    while binomial(H1+M,M-1) <= N
        H1 = H1 + 1
    end
    W = reduce(vcat, transpose.(combinations(1:H1+M-1,M-1))) .- repeat(0:M-2, 1, binomial(H1+M-1,M-1))' .- 1
    W = (hcat(W,zeros(size(W,1),1).+H1)-hcat(zeros(size(W,1),1),W))./H1
    if H1 < M
        H2 = 0
        while binomial(H1+M-1,M-1)+binomial(H2+M,M-1) <= N
            H2 = H2 + 1
        end
        if H2 > 0
            W2 = reduce(hcat,collect(combinations(1:H2+M-1,M-1)))' - repeat(collect(0:M-2)',binomial(H2+M-1,M-1)) .- 1
            W2 = (hcat(W2,zeros(size(W2,1),1).+H2)-hcat(zeros(size(W2,1),1),W2))./H2
            W  = vcat(W, W2./2 .+ 1/(2*M))
        end
    end
    W = broadcast((v) -> max(v, 1e-6), W)
    N = size(W, 1)
    W, N
end