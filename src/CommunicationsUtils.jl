__precompile__()

module CommunicationsUtils

    using Printf


    # power unit conversion

    export fraction2db, db2fraction

    fraction2db(x) = 10*log10(x)
    db2fraction(x) = 10^(x/10)


    # information theory utils

    export Plog2Q, h2, ln
    export binarystring, hexadecimalstring

    Plog2Q(p, q) = p > 1e-300 ? p*log2(q) : 0.0
    h2(p) = -Plog2Q(p, p) -Plog2Q(1-p, 1-p)
    ln = log

    binarystring(num) = reduce(*, [ c == 0 ? "0" : "1" for c in reverse(digits(num, base=2)) ])
    binarystring(num, m) = reduce(*, [ c == 0 ? "0" : "1" for c in reverse(digits(num, base=2, pad=m)) ])
    hexadecimalstring(num) = @sprintf("%x", num)


    # Gaussian distribution

    export N_pdf, N_cdf
    using SpecialFunctions

    function N_pdf(μ, σ²)
        return x -> 1.0/sqrt(2π*σ²) * exp(-(x-μ)^2 / (2*σ²))
    end

    function N_cdf(μ, σ²)
        σ = sqrt(σ²)
        return x -> 1/2*(1+erf( (x-μ)/(σ*sqrt(2)) ))
    end


    # conversion of channel parameters for AWGN-type channels

    export esn02σ², σ²2esn0
    export esn0db2σ², σ²2esn0db

    esn02σ²(esn0) = 1/2 * 1/esn0
    σ²2esn0(σ²) = 1/(2*σ²)
    esn0db2σ² = esn02σ² ∘ db2fraction
    σ²2esn0db = fraction2db ∘ σ²2esn0

    export esn0rates2ebn0

    esn0rates2ebn0(esn0s, rates) = (esn0s./rates, rates)

    export ebn02esn0, esn02ebn0, ebn0db2esn0db, esn0db2ebn0db

    ebn02esn0(ebn0, R) = ebn0 * R
    esn02ebn0(esn0, R) = esn0 * 1/R
    ebn0db2esn0db(ebn0db, R) = fraction2db(ebn02esn0(db2fraction(ebn0db), R))
    esn0db2ebn0db(esn0db, R) = fraction2db(esn02ebn0(db2fraction(esn0db), R))

    export ebn0db2σ², σ²2ebn0db

    ebn0db2σ²(ebn0db, R) = esn0db2σ²(ebn0db2esn0db(ebn0db, R))
    σ²2ebn0db(σ², R) = esn0db2ebn0db(σ²2esn0db(σ²), R)

    # some channel capacities

    export BSC_capacity, BEC_capacity, BiAWGN_capacity, AWGN_capacity

    BSC_capacity(p) = 1 - h2(p)
    BEC_capacity(e) = 1 - e

    function BiAWGN_capacity(σ²; Δ=0.01, lims=(-20, +20))
        esn0db = σ²2esn0db(σ²)
        @assert -15 <= esn0db <= +20
        @assert length(lims) == 2
        (lim_min, lim_max) = lims

        pYgM1 = N_pdf(-1, σ²)
        pYgP1 = N_pdf(+1, σ²)
        pY = y -> 1/2 * (pYgM1(y) + pYgP1(y))
        MIlogM1 = y -> Plog2Q(pYgM1(y), pYgM1(y)/pY(y))
        MIlogP1 = y -> Plog2Q(pYgP1(y), pYgP1(y)/pY(y))
        MIlog = y -> 1/2 * (MIlogM1(y) + MIlogP1(y))

        return sum([ Δ*MIlog(y) for y in lim_min:Δ:lim_max ])
    end

    function AWGN_capacity(σ²)
        esn0 = 1 / (2*σ²)
        return log2(1 + esn0)
    end

end
