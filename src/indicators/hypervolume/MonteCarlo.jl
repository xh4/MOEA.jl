function hv_monte_carlo(front, referencePoint)
    obj = reduce(hcat, front)'
    nsamples = 1_000_000
    maxvalue = referencePoint'
    minvalue = minimum(obj, dims=1)
    samples = reduce(hcat, map((max, min) -> rand(Uniform(min, max), nsamples), maxvalue, minvalue))
    M = length(referencePoint)
    for i = 1:size(obj)[1]
        domi = fill(true, size(samples)[1])
        m = 1
        while m <= M && any(domi)
            domi = domi .&& obj[i,m] .<= samples[:,m]
            m = m + 1
        end
        samples = samples[Not(domi),:]
    end
    prod(maxvalue-minvalue)*(1-size(samples)[1]/nsamples)
end