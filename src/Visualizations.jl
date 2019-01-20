__precompile__()

module Visualizations

    using Plots
    using Polar.DensityEvolution


    export plot_errorprobability, plot_errorprobability_annotated

    function plot_errorprobability(res_de)
        N = length(res_de)

        p = plot(xlabel="Channel", ylabel="Bit Error Probability", legend=:topright, size=(900, 500), xticks=N*(0.0:1/8:1.0))
        plot!(0:(N-1), errorprobability.(res_de), label="\$P_{error}\$", linealpha=0, marker=:auto, markerstrokewidth=0, markerstrokealpha=0)
        return p
    end

    function plot_errorprobability_annotated(res_de, frozen_bits)
        N = length(res_de)
        idx_frozen = findall(frozen_bits)
        idx_unfrozen = findall(.! frozen_bits)
        ch_frozen = res_de[idx_frozen]
        ch_unfrozen = res_de[idx_unfrozen]
        min_frozen = minimum(errorprobability.(ch_frozen))
        max_unfrozen = maximum(errorprobability.(ch_unfrozen))

        p = plot(xlabel="Channel", ylabel="Bit Error Probability \$P_{error}\$", legend=:topright, size=(900, 500), xticks=N*(0.0:1/8:1.0))
        plot!(idx_frozen .- 1, errorprobability.(ch_frozen), label="Frozen channels", linealpha=0, marker=:auto, markerstrokewidth=0, markerstrokealpha=0)
        plot!(idx_unfrozen .- 1, errorprobability.(ch_unfrozen), label="Unfrozen channels", linealpha=0, marker=:auto, markerstrokewidth=0, markerstrokealpha=0)
        plot!(0:(N-1), x -> min_frozen, label="Best of $(length(idx_frozen)) frozen channels")
        plot!(0:(N-1), x -> max_unfrozen, label="Worst of $(length(idx_unfrozen)) unfrozen channels")
        return p
    end


    export plot_mutualinformation, plot_mutualinformation_annotated

    function plot_mutualinformation(res_de)
        N = length(res_de)

        p = plot(xlabel="Channel", ylabel="Mutual Information", legend=:bottomright, size=(900, 500), xticks=N*(0.0:1/8:1.0))
        plot!(0:(N-1), mutualinformation.(res_de), label="\$I(.)\$", linealpha=0, marker=:auto, markerstrokewidth=0, markerstrokealpha=0)
        return p
    end

    function plot_mutualinformation_annotated(res_de, frozen_bits)
        N = length(res_de)
        idx_frozen = findall(frozen_bits)
        idx_unfrozen = findall(.! frozen_bits)
        ch_frozen = res_de[idx_frozen]
        ch_unfrozen = res_de[idx_unfrozen]
        max_frozen = maximum(mutualinformation.(ch_frozen))
        min_unfrozen = minimum(mutualinformation.(ch_unfrozen))

        p = plot(xlabel="Channel", ylabel="Mutual Information", legend=:bottomright, size=(900, 500), xticks=N*(0.0:1/8:1.0))
        plot!(idx_frozen .- 1, mutualinformation.(ch_frozen), label="Frozen channels", linealpha=0, marker=:auto, markerstrokewidth=0, markerstrokealpha=0)
        plot!(idx_unfrozen .- 1, mutualinformation.(ch_unfrozen), label="Unfrozen channels", linealpha=0, marker=:auto, markerstrokewidth=0, markerstrokealpha=0)
        plot!(0:(N-1), x -> max_frozen, label="Best of $(length(idx_frozen)) frozen channels")
        plot!(0:(N-1), x -> min_unfrozen, label="Worst of $(length(idx_unfrozen)) unfrozen channels")
        return p
    end


    export plot_codecomparison, plot_code

    function plot_codecomparison(c1, c2)
        @assert length(c1) == length(c2)
        # @assert count(c1) == count(c2)
        
        n = length(c1)
        idx1_frozen = findall(c1) .- 1
        idx1_unfrozen = findall(.! c1) .- 1
        idx2_frozen = findall(c2) .- 1
        idx2_unfrozen = findall(.! c2) .- 1
        idx_diff = findall(c1 .!= c2) .- 1
        idx_same = findall(c1 .== c2) .- 1
        idx_tofrozen = findall((c1 .== false) .& (c2 .== true)) .- 1
        idx_tounfrozen = findall((c1 .== true) .& (c2 .== false)) .- 1
        
        p = plot(xlim=(0-3, n-1+3), ylim=(1, 5), xticks=n*(0.0:1/8:1.0), size=(900, 100), legend=false, yticks=nothing)
        plot!(idx1_frozen, x -> 4, linewidth=0, marker=:circle, markercolor=:black, markerstrokewidth=0.3)
        plot!(idx1_unfrozen, x -> 4, linewidth=0, marker=:circle, markercolor=:white, markerstrokewidth=0.3)
        plot!(idx2_frozen, x -> 3, linewidth=0, marker=:circle, markercolor=:black, markerstrokewidth=0.3)
        plot!(idx2_unfrozen, x -> 3, linewidth=0, marker=:circle, markercolor=:white, markerstrokewidth=0.3)
        plot!(idx_tofrozen, x -> 2, linewidth=0, marker=:diamond, markercolor=:red, markerstrokecolor=:red, markerstrokewidth=1, markersize=3)
        plot!(idx_tounfrozen, x -> 2, linewidth=0, marker=:diamond, markercolor=:white, markerstrokecolor=:red, markerstrokewidth=1, markersize=3)
        return p
    end

    function plot_code(c1)
        n = length(c1)
        idx1_frozen = findall(c1) .- 1
        idx1_unfrozen = findall(.! c1) .- 1
        
        p = plot(xlim=(0-3, n-1+3), ylim=(1, 3), xticks=n*(0.0:1/8:1.0), size=(900, 50), legend=false, yticks=nothing)
        plot!(idx1_frozen, x -> 2, linewidth=0, marker=:circle, markercolor=:black, markerstrokewidth=0.3)
        plot!(idx1_unfrozen, x -> 2, linewidth=0, marker=:circle, markercolor=:white, markerstrokewidth=0.3)
        return p
    end

end
