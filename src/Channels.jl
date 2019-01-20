__precompile__()

module Channels

    using SpecialFunctions
    using Polar.GF2n: GF2Element, toint
    using Polar.CommunicationsUtils


    # channel alphabets

    export TernaryMΔ, Ternary0, TernaryPΔ
    @enum TernarySymbol TernaryMΔ=1 Ternary0=2 TernaryPΔ=3

    export ContinuousAlphabet, BinaryAlphabet, TernaryAlphabet
    const ContinuousAlphabet = Float64
    const BinaryAlphabet = GF2Element
    const TernaryAlphabet = TernarySymbol


    # modulation

    export mod_bpsk

    function mod_bpsk(x::GF2Element)
        return 1.0-2.0*toint(x)
    end


    # channels

    export AbstractCommunicationsChannel
    export AWGNChannel, BiAWGNChannel
    export BEChannel
    export intype, outtype, get_w, get_f, get_llr

    abstract type AbstractCommunicationsChannel end


    struct AWGNChannel <: AbstractCommunicationsChannel
        σ²::Float64
        _a::Float64
        _b::Float64
        _c::Float64

        function AWGNChannel(σ²::Float64)
            _a = 1/sqrt(2π*σ²)
            _b = -1/(2*σ²)
            _c = sqrt(σ²)
            return new(σ², _a, _b, _c)
        end
    end
    intype(::Type{AWGNChannel}) = ContinuousAlphabet
    outtype(::Type{AWGNChannel}) = ContinuousAlphabet
    get_w(ch::AWGNChannel, x::ContinuousAlphabet) = y::ContinuousAlphabet -> (ch._a * exp(ch._b*(x-y)^2))::Float64
    get_f(ch::AWGNChannel) = x::ContinuousAlphabet -> (x + ch._c*randn())::ContinuousAlphabet


    struct BiAWGNChannel <: AbstractCommunicationsChannel
        ch::AWGNChannel

        function BiAWGNChannel(σ²::Float64)
            ch = AWGNChannel(σ²)
            return new(ch)
        end
    end
    intype(::Type{BiAWGNChannel}) = BinaryAlphabet
    outtype(::Type{BiAWGNChannel}) = ContinuousAlphabet
    function get_w(ch::BiAWGNChannel, x::BinaryAlphabet)
        return get_w(ch.ch, mod_bpsk(x))
    end
    function get_f(ch::BiAWGNChannel)
        f = get_f(ch.ch)
        return x::BinaryAlphabet -> f(mod_bpsk(x))
    end
    function get_llr(ch::BiAWGNChannel)
        _a = 2/ch.ch.σ²
        function llr(y::ContinuousAlphabet)
            return _a * y
        end
        return llr
    end


    struct BEChannel <: AbstractCommunicationsChannel
        ϵ::Float64

        function BEChannel(ϵ::Float64)
            @assert 0.0 <= ϵ <= 1.0
            return new(ϵ)
        end
    end
    intype(::Type{BEChannel}) = BinaryAlphabet
    outtype(::Type{BEChannel}) = TernaryAlphabet
    function get_w(ch::BEChannel, x::BinaryAlphabet)
        function w(y::TernaryAlphabet)::Float64
            if y == TernaryPΔ
                return x == GF2Element(0) ? 1-ch.ϵ : 0.0
            elseif y == Ternary0
                return ch.ϵ
            else # y == TernaryMΔ
                return x == GF2Element(1) ? 1-ch.ϵ : 0.0
            end
        end
        return w
    end
    function get_f(ch::BEChannel)
        function f(x::BinaryAlphabet)::TernaryAlphabet
            if rand() <= ch.ϵ
                return Ternary0
            else
                return x == GF2Element(0) ? TernaryPΔ : TernaryMΔ
            end
        end
        return f
    end
    function get_llr(ch::BEChannel)
        function llr(y::TernaryAlphabet)
            if y == TernaryMΔ
                return -Inf
            elseif y == TernaryPΔ
                return +Inf
            else
                return 0.0
            end
        end
        return llr
    end

end
