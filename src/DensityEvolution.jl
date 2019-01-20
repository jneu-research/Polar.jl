__precompile__()

module DensityEvolution

    using Distributed
    using Polar.Distributions
    using Polar.CommunicationsUtils


    # value types for the density evolution

    export DEValueType
    export f_vn, f_cn
    import Base: isless, isapprox, show

    abstract type DEValueType end


    # float values

    export FloatValue

    struct FloatValue <: DEValueType
        v :: Float64
    end

    show(io::IO, x::FloatValue) = show(io, x.v)
    isless(x::FloatValue, y::FloatValue) = isless(x.v, y.v)
    isless(x::FloatValue, y::Float64) = isless(x.v, y)
    isapprox(x::FloatValue, y::Float64) = isapprox(x.v, y)

    f_vn(v1::FloatValue, v2::FloatValue, i) = FloatValue(v1.v + v2.v)
    f_cn(v1::FloatValue, v2::FloatValue, i) = FloatValue(sign(v1.v) * sign(v2.v) * min(abs(v1.v), abs(v2.v)))


    # "post p(olar) p(rocessed) value"

    export PostPPValue

    struct PostPPValue{T<:DEValueType,f_map} <: DEValueType
        v :: T
    end

    show(io::IO, x::PostPPValue) = show(io, x.v)
    isless(x::PostPPValue, y::PostPPValue) = isless(x.v, y.v)
    isless(x::PostPPValue, y::Float64) = isless(x.v, y)
    isapprox(x::PostPPValue, y::Float64) = isapprox(x.v, y)

    f_vn(v1::PostPPValue{T,f_map}, v2::PostPPValue{T,f_map}, i) where {T,f_map} = PostPPValue{T,f_map}( f_map( f_vn(v1.v, v2.v, i) , i) )
    f_cn(v1::PostPPValue{T,f_map}, v2::PostPPValue{T,f_map}, i) where {T,f_map} = PostPPValue{T,f_map}( f_map( f_cn(v1.v, v2.v, i) , i) )


    # compound values

    f_vn(v1::Tuple, v2::Tuple, i) = Tuple(map((e1, e2) -> f_vn(e1, e2, i), v1, v2))
    f_cn(v1::Tuple, v2::Tuple, i) = Tuple(map((e1, e2) -> f_cn(e1, e2, i), v1, v2))


    # distributions for DE

    export DEDistributionType
    export de_vn, de_cn, errorprobability, mutualinformation

    abstract type DEDistributionType end


    # Gaussian approximation

    export GADEDistribution, J, Jinv

    struct GADEDistribution <: DEDistributionType
        I :: Float64
    end

    # see Brännström et al.
    J_H1 = 0.3073
    J_H2 = 0.8935
    J_H3 = 1.1064

    J(σ) = (1-2^(-J_H1*σ^(2*J_H2)))^J_H3
    Jinv(I) = (-1/J_H1*log2(1-I^(1/J_H3)))^(1/(2*J_H2))
    # σ2I = J
    # I2σ = Jinv

    de_vn(d::GADEDistribution, i) = GADEDistribution( J(sqrt(2)*Jinv(d.I)) )
    de_cn(d::GADEDistribution, i) = GADEDistribution( 1-J(sqrt(2)*Jinv(1-d.I)) )

    function errorprobability(d :: GADEDistribution)
        I = d.I

        if I ≈ 0.0
            return 0.5
        elseif I ≈ 1.0
            return 0.0
        end

        σ = Jinv(I)
        σ² = σ^2
        μ = σ²/2
        F = N_cdf(μ, σ²)

        return F(0)
    end

    function mutualinformation(d :: GADEDistribution)
        return d.I
    end


    # BEC distribution

    export BECDistribution

    struct BECDistribution <: DEDistributionType
        ϵ :: Float64
    end

    de_vn(d::BECDistribution, i) = BECDistribution(d.ϵ^2)
    de_cn(d::BECDistribution, i) = BECDistribution(2*d.ϵ-d.ϵ^2)

    function errorprobability(d :: BECDistribution)
        return 1//2 * d.ϵ
    end

    function mutualinformation(d :: BECDistribution)
        return 1 - d.ϵ
    end


    # discrete distribution

    export DiscreteDEDistribution
    import Base: iterate, length, eltype

    struct DiscreteDEDistribution{T} <: DEDistributionType
        dist :: Dict{T,Float64}
    end

    function DiscreteDEDistribution(d::Distribution{Tuple{T},1}) where {T}
        d_new = DiscreteDEDistribution{T}(Dict(
            v[1] => p for (v, p) in d.data
        ))
        return d_new
    end

    function _transformrawdistribution(f, d::Tdist) where {T,Tdist<:Dict{T,Float64}}
        @assert isdistribution(d)
        d_result = Tdist()
        for (val, prob) in d
            val_new = f(val)
            d_result[val_new] = get(d_result, val_new, 0.0) + prob
        end
        @assert isdistribution(d_result)
        return d_result
    end

    function _transformrawdistribution(f, d1::Tdist, d2::Tdist) where {T,Tdist<:Dict{T,Float64}}
        @assert isdistribution(d1)
        @assert isdistribution(d2)
        d_result = Tdist()
        for (val1, prob1) in d1
            for (val2, prob2) in d2
                val_new = f(val1, val2)
                d_result[val_new] = get(d_result, val_new, 0.0) + prob1*prob2
            end
        end
        @assert isdistribution(d_result)
        return d_result
    end

    de_vn(d::DiscreteDEDistribution{T}, i) where {T} = DiscreteDEDistribution{T}(_transformrawdistribution((e1, e2) -> f_vn(e1, e2, i), d.dist, d.dist))
    de_cn(d::DiscreteDEDistribution{T}, i) where {T} = DiscreteDEDistribution{T}(_transformrawdistribution((e1, e2) -> f_cn(e1, e2, i), d.dist, d.dist))

    function errorprobability(d :: DiscreteDEDistribution)
        p = 0.0
        for (val, prob) in d.dist
            if val < 0.0
                p += prob
            elseif val ≈ 0.0
                p += 0.5 * prob
            end
        end
        return p
    end

    iterate(d::DiscreteDEDistribution) = iterate(d.dist)
    iterate(d::DiscreteDEDistribution, state) = iterate(d.dist, state)
    length(d::DiscreteDEDistribution) = length(d.dist)
    eltype(d::DiscreteDEDistribution) = eltype(d.dist)


    # actual density evolution implementation

    export densityevolution, densityevolution_raw

    function densityevolution_raw_extend(k, ch; debug=false)
        if debug
            println("Extending $(k)...")
        end

        return Dict(
            "$(k)0" => de_cn(ch, length(k)),
            "$(k)1" => de_vn(ch, length(k)),
        )
    end

    function densityevolution_raw(m, d; debug=false)
        channels = Dict{String,typeof(d)}("" => d)
        
        for i in 1:m
            keys_to_extend = collect(keys(channels))
            channels = reduce(merge, pmap(k -> densityevolution_raw_extend(k, channels[k]; debug=debug), keys_to_extend))
        end
        
        return channels
    end

    function densityevolution_raw_full(m, d; debug=false)
        channels = Dict("" => d)
        
        for i in 1:m
            keys_to_extend = collect(Iterators.filter(k->length(k) == i-1, keys(channels)))
            merge!(channels, reduce(merge, pmap(k -> densityevolution_raw_extend(k, channels[k]; debug=debug), keys_to_extend)))
        end
        
        return channels
    end

    function densityevolution(m, d; debug=false)
        result = Vector{typeof(d)}(undef, 2^m)
        channels = densityevolution_raw(m, d; debug=debug)
        for idx in 0:(2^m-1)
            result[idx+1] = channels[binarystring(idx, m)]
        end
        return result
    end


    export freeze_by_errorprobability, freeze_by_mutualinformation

    function freeze_by_errorprobability(res_de, k_)
        n = length(res_de)
        k = n - k_

        res_de_proberrs = errorprobability.(res_de)
        res_de_proberrs_sorted = sort(res_de_proberrs; rev=false)
        thr = res_de_proberrs_sorted[k]
        frozen = res_de_proberrs .> thr

        @assert count(frozen) == k_
        return frozen
    end

    function freeze_by_mutualinformation(res_de, k_)
        n = length(res_de)
        k = n - k_

        res_de_proberrs = mutualinformation.(res_de)
        res_de_proberrs_sorted = sort(res_de_proberrs; rev=true)
        thr = res_de_proberrs_sorted[k]
        frozen = res_de_proberrs .< thr

        @assert count(frozen) == k_
        return frozen
    end

end
