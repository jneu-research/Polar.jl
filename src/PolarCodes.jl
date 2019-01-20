__precompile__()

module PolarCodes

    using Polar.GF2n: GF2Element, GF2_0, GF2_1
    using Polar.Utils: @idx0
    using Polar.CommunicationsUtils: binarystring, hexadecimalstring

    import Random: rand!, rand

    export PolarCode
    export n, k, R, m
    export interleave!, interleave
    export deinterleave!, deinterleave
    export polarencode!, polarencode
    export polardecode!, polardecode
    export bitreverse!, bitreverse


    struct PolarCode{n,k}
        is_frozen::NTuple{n,Bool}
        frozen_values::NTuple{n,GF2Element}
        bitrev_permutation::NTuple{n,Int}

        function PolarCode{n,k}(is_frozen::NTuple{n,Bool}, frozen_values::NTuple{n,GF2Element}) where {n,k}
            @idx0 begin
                @assert n == nextpow(2, n)
                @assert count(is_frozen) == n-k

                m = Int(log2(n))
                bitrev_permutation = Vector{Int}(undef, n)
                for i in 0:(n-1)
                    bitrev_permutation[i] = parse(Int, reverse(binarystring(i, m)); base=2)
                end

                return new(is_frozen, frozen_values, Tuple(bitrev_permutation))
            end
        end
    end

    function PolarCode(is_frozen::Union{Vector{Bool},BitArray{1}}, frozen_values::Vector{GF2Element})
        @assert length(is_frozen) == length(frozen_values)
        n = length(is_frozen)
        k = n - count(is_frozen)
        return PolarCode{n,k}(Tuple(is_frozen), Tuple(frozen_values))
    end

    function PolarCode(is_frozen::Union{Vector{Bool},BitArray{1}})
        n = length(is_frozen)
        return PolarCode(is_frozen, [ GF2_0 for i in 1:n ])
    end


    function Base.show(io::IO, m::PolarCode{n,k}) where {n,k}
        frozen_id = join([ (x == GF2_1 ? '1' : '0') for x in m.is_frozen ])
        values_id = join([ (x == GF2_1 ? '1' : '0') for x in m.frozen_values ])
        id = hash("$(frozen_id)-$(values_id)")
        print(io, "$(typeof(m))(fingerprint=$(hexadecimalstring(id)))")
    end

    @inline n(code::PolarCode{n_,k_}) where {n_,k_} = n_
    @inline k(code::PolarCode{n_,k_}) where {n_,k_} = k_
    @inline R(code::PolarCode{n_,k_}) where {n_,k_} = k_/n_
    @inline m(code::PolarCode{n_,k_}) where {n_,k_} = Int(log2(n(code)))


    function interleave!(code::PolarCode, dst::Vector{GF2Element}, src::Vector{GF2Element})
        @idx0 begin
            @assert length(dst) == n(code)
            @assert length(src) == k(code)
            idx_src = 0
            for idx_dst in 0:(n(code)-1)
                if code.is_frozen[idx_dst]
                    dst[idx_dst] = code.frozen_values[idx_dst]
                else
                    dst[idx_dst] = src[idx_src]
                    idx_src += 1
                end
            end
            @assert idx_src == k(code)
            return nothing
        end
    end

    function deinterleave!(code::PolarCode, dst::Vector{GF2Element}, src::Vector{GF2Element})
        @idx0 begin
            @assert length(dst) == k(code)
            @assert length(src) == n(code)
            idx_dst = 0
            for idx_src in 0:(n(code)-1)
                if code.is_frozen[idx_src]
                    @assert src[idx_src] == code.frozen_values[idx_src]
                else
                    dst[idx_dst] = src[idx_src]
                    idx_dst += 1
                end
            end
            @assert idx_dst == k(code)
            return nothing
        end
    end

    function _polarcode!(u::Vector{GF2Element})
        len = length(u)
        @assert len == nextpow(2, len)
        _polarcode!(u, 0, len)
    end

    function _polarcode!(u::Vector{GF2Element}, idx::Int, len::Int)
        @idx0 @inbounds begin
            if len == 2
                u[idx] = u[idx] + u[idx+1]
            else
                len′ = len>>1
                _polarcode!(u, idx, len′)
                _polarcode!(u, idx+len′, len′)
                for i in 0:(len′-1)
                    u[idx+i] = u[idx+i] + u[idx+i+len′]
                end
            end
            return nothing
        end
    end

    function polarencode!(code::PolarCode, dst::Vector{GF2Element}, src::Vector{GF2Element})
        @assert length(dst) == n(code)
        @assert length(src) == k(code)
        interleave!(code, dst, src)
        bitreverse!(code, dst)
        _polarcode!(dst)
        return nothing
    end

    function polardecode!(code::PolarCode, dst::Vector{GF2Element}, src::Vector{GF2Element})
        @assert length(dst) == k(code)
        @assert length(src) == n(code)
        _polarcode!(src)
        bitreverse!(code, src)
        deinterleave!(code, dst, src)
        return nothing
    end

    function bitreverse!(code::PolarCode, dst::Vector{GF2Element})
        @idx0 begin
            for idx_i in 0:(n(code)-1)
                idx_j = code.bitrev_permutation[idx_i]
                if idx_j > idx_i
                    tmp_a = dst[idx_i]
                    tmp_b = dst[idx_j]
                    dst[idx_i] = tmp_b
                    dst[idx_j] = tmp_a
                end
            end
            return nothing
        end
    end

    function rand!(code::PolarCode, dst::Vector{GF2Element})
        @assert length(dst) == k(code)
        for i in 1:k(code)
            dst[i] = rand(GF2Element)
        end
        return nothing
    end


    function interleave(code::PolarCode, src::Vector{GF2Element})
        dst = Vector{GF2Element}(undef, n(code))
        interleave!(code, dst, src)
        return dst
    end

    function deinterleave(code::PolarCode, src::Vector{GF2Element})
        dst = Vector{GF2Element}(undef, k(code))
        deinterleave!(code, dst, src)
        return dst
    end

    function polarencode(code::PolarCode, src::Vector{GF2Element})
        dst = Vector{GF2Element}(undef, n(code))
        polarencode!(code, dst, src)
        return dst
    end

    function polardecode(code::PolarCode, src::Vector{GF2Element})
        dst = Vector{GF2Element}(undef, k(code))
        polardecode!(code, dst, src)
        return dst
    end

    function bitreverse(code::PolarCode, src::Vector{GF2Element})
        dst = copy(src)
        bitreverse!(code, dst)
        return dst
    end

    function rand(code::PolarCode)
        dst = Vector{GF2Element}(undef, k(code))
        rand!(code, dst)
        return dst
    end

end