__precompile__()

module BridgeCh2DE

    export channeloutputdistribution

    using Polar.GF2n
    using Polar.Channels
    using Polar.DensityEvolution
    using Polar.CommunicationsUtils
    using Polar.Utils


    function channeloutputdistribution(ch::BiAWGNChannel)
        σ² = ch.ch.σ²
        σ = sqrt(σ²)
        I = J(2/σ)
        return GADEDistribution(I)
    end

    function channeloutputdistribution(ch::Tch, v::Type{PostPPValue{Tllr,f_map}}) where {Tch<:AbstractCommunicationsChannel,Tllr,f_map}
        return DiscreteDEDistribution{PostPPValue{Tllr,f_map}}(Dict( PostPPValue{Tllr,f_map}(val) => prob for (val, prob) in channeloutputdistribution(ch, Tllr) ))
    end

    function channeloutputdistribution(ch::BSChannel)
        @assert !(ch.p ≈ 0.0)
        Δ = ln((1-ch.p)/ch.p)
        return DiscreteDEDistribution{FloatValue}(Dict( FloatValue(+Δ) => (1-ch.p), FloatValue(-Δ) => ch.p ))
    end

end
