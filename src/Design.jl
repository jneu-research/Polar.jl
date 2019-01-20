__precompile__()

module Design

    # design as Reed-Muller code

    export design_reedmuller

    function design_reedmuller(m, k)
        n = 2^m

        F = [1 1; 0 1]
        G = reduce(kron, [ F for i in 1:m ])
        G_weights = sum(G, dims=1)[:]
        G_weights_sorted = sort(G_weights; rev=false)
        threshold = G_weights_sorted[n-k]

        code = G_weights .<= threshold
        @assert count(code) == n-k
        return code
    end


    # design using density evolution

    export design_densityevolution
    using Polar.DensityEvolution

    function design_densityevolution(m, k, D)
        n = 2^m

        resde = densityevolution(m, D)
        code = freeze_by_errorprobability(resde, n-k)
        @assert count(code) == n-k
        return code
    end

end
