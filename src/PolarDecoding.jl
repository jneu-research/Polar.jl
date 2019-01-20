__precompile__()

module PolarDecoding

    using StaticArrays
    using Random

    using Polar.Stacks: FixedCapacityStack
    using Polar.Utils: @idx0, prettyshow
    using Polar.PolarCodes: PolarCode
    using Polar.GF2n: GF2Element, toint, GF2_0, GF2_1, ⊻
    using Polar.CommunicationsUtils: ln

    import Base: ⊻
    import Base: iterate, eltype, length
    import Base: isless
    import Base: zero
    import Base: sign, abs
    import Base: -
    

    # data types to be used in the polar code decoder

    export LLRType, UnquantizedLLR
    export cn_update, vn_update, pm_increment

    export PMType, UnquantizedPM
    export update

    abstract type LLRType end
    abstract type PMType end


    # unquantized LLR type

    struct UnquantizedLLR <: LLRType
        v::Float64
    end

    UnquantizedLLR() = UnquantizedLLR(0.00)
    zero(::Type{UnquantizedLLR}) = UnquantizedLLR()

    # min-approximation:
    function cn_update(a::UnquantizedLLR, b::UnquantizedLLR)
        return UnquantizedLLR(sign(a.v)*sign(b.v)*min(abs(a.v), abs(b.v)))
    end
    # # ideal:
    # function cn_update(a::UnquantizedLLR, b::UnquantizedLLR)
    #   return UnquantizedLLR(ln((exp(a.v+b.v)+1)/(exp(a.v)+exp(b.v))))
    # end

    function vn_update(a::UnquantizedLLR, b::UnquantizedLLR, c::GF2Element)
        if c == GF2_0
            return UnquantizedLLR( a.v+b.v)
        else # c == GF2_1
            return UnquantizedLLR(-a.v+b.v)
        end
    end

    # piecewise-linear approximation:
    function pm_increment(a::UnquantizedLLR, i::Int)
        Δ = 2*ln(2)
        λ = a.v
        if λ <= -Δ
            return -λ
        elseif λ <= +Δ
            return -(λ-Δ)/2
        else
            return 0.0
        end
    end
    # # coarser piecewise-linear approximation (problematic around 0 with quantized LLRs!)
    # function pm_increment(a::UnquantizedLLR, i::Int)
    #   if a.v >= 0.0
    #       return 0.0
    #   else # a.v < 0.0
    #       return -a.v
    #   end
    # end

    function -(a::UnquantizedLLR)
        return UnquantizedLLR(-a.v)
    end


    # unquantized PM type

    struct UnquantizedPM <: PMType
        v::Float64
    end

    UnquantizedPM() = UnquantizedPM(0.00)
    zero(::Type{UnquantizedPM}) = UnquantizedPM()

    function update(a::UnquantizedPM, b::T, c::GF2Element, i::Int) where {T<:LLRType}
        if c == GF2_0
            return UnquantizedPM(a.v + pm_increment( b, i))
        else # c == GF2_1
            return UnquantizedPM(a.v + pm_increment(-b, i))
        end
    end

    function isless(x::UnquantizedPM, y::UnquantizedPM)
        return isless(x.v, y.v)
    end


    # logging for the SCL decoder

    export ObservationCounter, observe!, n

    mutable struct ObservationCounter{T}
        samples::Dict{T,Int}

        function ObservationCounter{T}() where {T}
            return new(Dict{T,Int}())
        end
    end

    function observe!(c::ObservationCounter{T}, v::T) where {T}
        c.samples[v] = get(c.samples, v, 0) + 1
        return nothing
    end

    function n(c::ObservationCounter)
        return sum(values(c.samples))
    end

    export AbstractLLRLogger, BasicLLRLogger, CompressedLLRLogger, SingleFrameLLRLogger

    abstract type AbstractLLRLogger{T_llr<:LLRType,n} end

    struct BasicLLRLogger{T_llr<:LLRType,n} <: AbstractLLRLogger{T_llr,n}
        counters::NTuple{n,ObservationCounter{T_llr}}

        function BasicLLRLogger{T_llr,n}() where {T_llr<:LLRType,n}
            @assert isa(n, Int)
            counters_ = [ ObservationCounter{T_llr}() for i in 1:n ]
            return new(Tuple(counters_))
        end
    end

    function log!(l::BasicLLRLogger{T_llr,n}, i::Int, v::T_llr) where {T_llr,n}
        @idx0 begin
            observe!(l.counters[i], v)
            return nothing
        end
    end

    struct CompressedLLRLogger{T_llr<:LLRType,n,f_compress,T_value} <: AbstractLLRLogger{T_llr,n}
        counters::NTuple{n,ObservationCounter{T_value}}

        function CompressedLLRLogger{T_llr,n,f_compress,T_value}() where {T_llr<:LLRType,n,f_compress,T_value}
            @assert isa(n, Int)
            counters_ = [ ObservationCounter{T_value}() for i in 1:n ]
            return new(Tuple(counters_))
        end
    end

    function log!(l::CompressedLLRLogger{T_llr,n,f_compress,T_value}, i::Int, v::T_llr) where {T_llr,n,f_compress,T_value}
        @idx0 begin
            observe!(l.counters[i], f_compress(v))
            return nothing
        end
    end

    struct SingleFrameLLRLogger{T_llr<:LLRType,n} <: AbstractLLRLogger{T_llr,n}
        llrs::Vector{T_llr}

        function SingleFrameLLRLogger{T_llr,n}() where {T_llr<:LLRType,n}
            @assert isa(n, Int)
            return new(Vector{T_llr}(undef, n))
        end

        function SingleFrameLLRLogger{T_llr,n}(v::Vector{T_llr}) where {T_llr<:LLRType,n}
            @assert isa(n, Int)
            @assert length(v) == n
            return new(v)
        end
    end

    function log!(l::SingleFrameLLRLogger{T_llr,n}, i::Int, v::T_llr) where {T_llr,n}
        @idx0 begin
            l.llrs[i] = v
            return nothing
        end
    end


    # the SCL decoder for polar codes

    export SCLDecoder
    export decode!, decode_with_debug!
    export reconstruct_pm!, reconstruct_pm_with_debug!
    export extract_list, extract_list!, extract_list_C, extract_list_C!
    export extract_best_pm, extract_best_pm!, extract_best_pm_C, extract_best_pm_C!
    export llrtype, pmtype

    struct SCLDecoder{m,L,T_llr,T_pm,code}
        stack_unused_pis::FixedCapacityStack{Int,L}
        bitmap_used_pis::Vector{Bool}
        list_pms::Vector{T_pm}
        array_u::Matrix{GF2Element}
        array_p::Matrix{Vector{T_llr}}
        array_c::Matrix{Matrix{GF2Element}}
        lut_pi2ai::Matrix{Int}
        stacks_unused_ais::Vector{FixedCapacityStack{Int,L}}
        array_refcnt::Matrix{Int}

        tmp_potpaths::Matrix{T_pm}
        tmp_potpaths_cont::Matrix{Bool}
        tmp_potpaths_pms::Vector{T_pm}

        function SCLDecoder{m,L,T_llr,T_pm,code}() where {m,L,T_llr<:LLRType,T_pm<:PMType,code}
            @assert isa(m, Int)
            @assert isa(L, Int)
            @assert isa(code, PolarCode{2^m})

            stack_unused_pis = FixedCapacityStack{Int,L}()
            bitmap_used_pis = [ false for ℓ in 0:(L-1) ]
            list_pms = [ T_pm() for ℓ in 0:(L-1) ]
            array_u = [ GF2_0 for i in 0:(2^m-1), ℓ in 0:(L-1) ]
            array_p = [ [ T_llr() for i in 0:(2^(m-λ)-1) ] for λ in 0:m, s in 0:(L-1) ]
            array_c = [ [ GF2_0 for i in 0:(2^(m-λ)-1), u in 0:1 ] for λ in 0:m, s in 0:(L-1) ]
            lut_pi2ai = [ Int(0) for λ in 0:m, s in 0:(L-1) ]
            stacks_unused_ais = [ FixedCapacityStack{Int,L}() for λ in 0:m ]
            array_refcnt = [ Int(0) for λ in 0:m, s in 0:(L-1) ]

            tmp_potpaths = [ T_pm() for ℓ in 0:(L-1), u in 0:1 ]
            tmp_potpaths_cont = [ false for ℓ in 0:(L-1), u in 0:1 ]
            tmp_potpaths_pms = Vector{T_pm}()

            return new(
                stack_unused_pis,
                bitmap_used_pis,
                list_pms,
                array_u,
                array_p,
                array_c,
                lut_pi2ai,
                stacks_unused_ais,
                array_refcnt,

                tmp_potpaths,
                tmp_potpaths_cont,
                tmp_potpaths_pms,
            )
        end
    end

    function Base.show(io::IO, d::SCLDecoder{m,L,T_llr,T_pm,code}) where {m,L,T_llr,T_pm,code}
        print(io, "$(typeof(d))()")
    end

    @inline n(d::SCLDecoder{m,L,T_llr,T_pm,code}) where {m,L,T_llr,T_pm,code} = Int(2^m)
    @inline m(d::SCLDecoder{m_,L_,T_llr_,T_pm_,code_}) where {m_,L_,T_llr_,T_pm_,code_} = Int(m_)
    @inline L(d::SCLDecoder{m_,L_,T_llr_,T_pm_,code_}) where {m_,L_,T_llr_,T_pm_,code_} = Int(L_)
    @inline code(d::SCLDecoder{m_,L_,T_llr_,T_pm_,code_}) where {m_,L_,T_llr_,T_pm_,code_} = code_
    llrtype(d::SCLDecoder{m_,L_,T_llr_,T_pm_,code_}) where {m_,L_,T_llr_,T_pm_,code_} = T_llr_
    pmtype(d::SCLDecoder{m_,L_,T_llr_,T_pm_,code_}) where {m_,L_,T_llr_,T_pm_,code_} = T_pm_


    # initialize the data structures of the decoder

    function reset!(d::SCLDecoder{m,L,T_llr,T_pm,code}) where {m,L,T_llr,T_pm,code}
        @idx0 @inbounds begin
            # # unnecessary reset steps
            # for λ in 0:m, s in 0:(L-1), i in 0:(2^(m-λ)-1)
            #   d.array_p[λ,s][i] = T_llr()
            #   d.array_c[λ,s][i,0] = GF2_0
            #   d.array_c[λ,s][i,1] = GF2_0
            # end
            # for λ in 0:m, s in 0:(L-1)
            #   d.lut_pi2ai[λ,s] = 0
            # end
            # for ℓ in 0:(L-1)
            #   d.list_pms[ℓ] = T_pm()
            #   for i in 0:(2^m-1)
            #       d.array_u[i,ℓ] = GF2_0
            #   end
            # end

            # necessary reset steps
            empty!(d.stack_unused_pis)
            for λ in 0:m
                empty!(d.stacks_unused_ais[λ])
            end
            for ℓ in 0:(L-1)
                d.bitmap_used_pis[ℓ] = false
                push!(d.stack_unused_pis, ℓ)
            end
            for λ in 0:m, s in 0:(L-1)
                d.array_refcnt[λ,s] = 0
                push!(d.stacks_unused_ais[λ], s)
            end
        end
    end


    # low-level allocation/freeing of path and array indices

    @inline current_list_size(d::SCLDecoder{m,L,T_llr,T_pm,code}) where {m,L,T_llr,T_pm,code} = L-length(d.stack_unused_pis)

    function alloc_pi!(d::SCLDecoder)
        @idx0 @inbounds begin
            ℓ = pop!(d.stack_unused_pis)
            # @assert d.bitmap_used_pis[ℓ] == false
            d.bitmap_used_pis[ℓ] = true
            return ℓ
        end
    end

    function free_pi!(d::SCLDecoder, ℓ::Int)
        @idx0 @inbounds begin
            # @assert d.bitmap_used_pis[ℓ] == true
            push!(d.stack_unused_pis, ℓ)
            d.bitmap_used_pis[ℓ] = false
        end
    end

    function alloc_ai!(d::SCLDecoder, λ::Int)
        @idx0 @inbounds begin
            s = pop!(d.stacks_unused_ais[λ])
            # @assert d.array_refcnt[λ,s] == 0
            d.array_refcnt[λ,s] = 1
            return s
        end
    end

    function free_ai!(d::SCLDecoder, λ::Int, s::Int)
        @idx0 @inbounds begin
            # @assert d.array_refcnt[λ,s] == 1
            push!(d.stacks_unused_ais[λ], s)
            d.array_refcnt[λ,s] = 0
        end
    end


    # low-level iteration over active paths

    struct ActivePIsIterator{m,L,T_llr,T_pm,code}
        d::SCLDecoder{m,L,T_llr,T_pm,code}
    end

    @inline active_pis(d::SCLDecoder{m,L,T_llr,T_pm,code}) where {m,L,T_llr,T_pm,code} = ActivePIsIterator(d)

    function _find_next_active_pi(d::SCLDecoder{m,L,T_llr,T_pm,code}, ℓ::Int) where {m,L,T_llr,T_pm,code}
        @idx0 @inbounds begin
            for i in (ℓ+1):(L-1)
                if d.bitmap_used_pis[i]
                    return i::Int
                end
            end
            return L::Int
        end
    end

    @inline eltype(i::ActivePIsIterator{m,L,T_llr,T_pm,code}) where {m,L,T_llr,T_pm,code} = Int
    @inline length(i::ActivePIsIterator{m,L,T_llr,T_pm,code}) where {m,L,T_llr,T_pm,code} = current_list_size(i.d)

    @inline function iterate(i::ActivePIsIterator{m,L,T_llr,T_pm,code}) where {m,L,T_llr,T_pm,code}
        item = _find_next_active_pi(i.d, -1)
        if item == L
            return nothing
        else
            return (item, _find_next_active_pi(i.d, item))
        end
    end
    @inline function iterate(i::ActivePIsIterator{m,L,T_llr,T_pm,code}, state_current) where {m,L,T_llr,T_pm,code}
        if state_current == L
            return nothing
        else
            return (state_current, _find_next_active_pi(i.d, state_current))
        end
    end


    # low-level functionality of the decoder

    function assign_initial_path!(d::SCLDecoder{m,L,T_llr,T_pm,code}) where {m,L,T_llr,T_pm,code}
        @idx0 @inbounds begin
            ℓ = alloc_pi!(d)
            d.list_pms[ℓ] = T_pm()
            for λ in 0:m
                s = alloc_ai!(d, λ)
                d.lut_pi2ai[λ,ℓ] = s
            end
            return ℓ
        end
    end

    function clone_path!(d::SCLDecoder{m,L,T_llr,T_pm,code}, ℓ::Int, φ::Int) where {m,L,T_llr,T_pm,code}
        @idx0 @inbounds begin
            ℓ′ = alloc_pi!(d)
            d.list_pms[ℓ′] = d.list_pms[ℓ]
            for i in 0:(φ-1)
                d.array_u[i,ℓ′] = d.array_u[i,ℓ]
            end
            for λ in 0:m
                s = d.lut_pi2ai[λ,ℓ]
                d.lut_pi2ai[λ,ℓ′] = s
                d.array_refcnt[λ,s] += 1
            end
            return ℓ′
        end
    end

    function kill_path!(d::SCLDecoder{m,L,T_llr,T_pm,code}, ℓ::Int) where {m,L,T_llr,T_pm,code}
        @idx0 @inbounds begin
            free_pi!(d, ℓ)
            for λ in 0:m
                s = d.lut_pi2ai[λ,ℓ]
                if d.array_refcnt[λ,s] == 1
                    free_ai!(d, λ, s)
                else
                    d.array_refcnt[λ,s] -= 1
                end
            end
        end
    end

    function lazy_copy_pc!(d::SCLDecoder{m,L,T_llr,T_pm,code}, λ::Int, s::Int) where {m,L,T_llr,T_pm,code}
        @idx0 @inbounds begin
            s′ = alloc_ai!(d, λ)
            copyto!(d.array_p[λ,s′], d.array_p[λ,s])
            copyto!(d.array_c[λ,s′], d.array_c[λ,s])
            d.array_refcnt[λ,s] -= 1
            return s′
        end
    end

    function cow_mechanism!(d::SCLDecoder{m,L,T_llr,T_pm,code}, λ::Int, ℓ::Int) where {m,L,T_llr,T_pm,code}
        @idx0 @inbounds begin
            s = d.lut_pi2ai[λ,ℓ]
            if d.array_refcnt[λ,s] == 1
                # data is not shared, return it
                s′ = s
            else
                # perform lazy copy, return the copy
                s′ = lazy_copy_pc!(d, λ, s)
                d.lut_pi2ai[λ,ℓ] = s′
            end
            return s′
        end
    end

    function get_p!(d::SCLDecoder{m,L,T_llr,T_pm,code}, λ::Int, ℓ::Int) where {m,L,T_llr,T_pm,code}
        @idx0 @inbounds begin
            s′ = cow_mechanism!(d, λ, ℓ)
            return d.array_p[λ,s′]::Vector{T_llr}
        end
    end

    function get_c!(d::SCLDecoder{m,L,T_llr,T_pm,code}, λ::Int, ℓ::Int) where {m,L,T_llr,T_pm,code}
        @idx0 @inbounds begin
            s′ = cow_mechanism!(d, λ, ℓ)
            return d.array_c[λ,s′]::Matrix{GF2Element}
        end
    end


    # mid-level functions to be used in the decoder

    function recursively_calc_p!(d::SCLDecoder{m,L,T_llr,T_pm,code}, λ::Int, φ::Int) where {m,L,T_llr,T_pm,code}
        @idx0 @inbounds begin
            if λ == 0
                # there is nothing to compute on the 0th layer, this layer contains the channel output
                return
            end

            # const ψ = φ>>1
            ψ = φ>>1
            if φ % 2 == 0
                # we just entered a branch that requires to recompute the layer "below" (λ-1)
                recursively_calc_p!(d, λ-1, ψ)
            end

            for ℓ in active_pis(d)
                P = get_p!(d, λ, ℓ)
                P′ = get_p!(d, λ-1, ℓ)
                C = get_c!(d, λ, ℓ)

                # # reorder loops for better parallelization:
                # for β in 0:(2^(m-λ)-1)
                #   if φ % 2 == 0
                #       P[β] = cn_update(P′[2β], P′[2β+1])
                #   else
                #       u′ = C[β,0]
                #       P[β] = vn_update(P′[2β], P′[2β+1], u′)
                #   end
                # end
                if φ % 2 == 0
                    @simd for β in 0:(2^(m-λ)-1)
                        P[β] = cn_update(P′[2β], P′[2β+1])
                    end
                else
                    @simd for β in 0:(2^(m-λ)-1)
                        u′ = C[β,0]
                        P[β] = vn_update(P′[2β], P′[2β+1], u′)
                    end
                end
            end

            return nothing
        end
    end

    function recursively_update_c!(d::SCLDecoder{m,L,T_llr,T_pm,code}, λ::Int, φ::Int) where {m,L,T_llr,T_pm,code}
        @idx0 @inbounds begin
            # @assert φ % 2 == 1
            # const ψ = φ>>1
            ψ = φ>>1

            for ℓ in active_pis(d)
                C = get_c!(d, λ, ℓ)
                C′ = get_c!(d, λ-1, ℓ)
                @simd for β in 0:(2^(m-λ)-1)
                    C′[2β,ψ%2] = C[β,0] ⊻ C[β,1]
                    C′[2β+1,ψ%2] = C[β,1]
                end
            end

            if ψ % 2 == 1
                recursively_update_c!(d, λ-1, ψ)
            end

            return nothing
        end
    end


    # decoder main functions

    function continue_paths_frozen!(d::SCLDecoder{m,L,T_llr,T_pm,code}, φ::Int) where {m,L,T_llr,T_pm,code}
        @idx0 @inbounds begin
            for ℓ in active_pis(d)
                P = get_p!(d, m, ℓ)
                C = get_c!(d, m, ℓ)
                C[0,φ%2] = code.frozen_values[φ]
                d.array_u[φ,ℓ] = code.frozen_values[φ]
                d.list_pms[ℓ] = update(d.list_pms[ℓ], P[0], code.frozen_values[φ], φ)
            end

            return nothing
        end
    end

    function continue_paths_unfrozen!(d::SCLDecoder{m,L,T_llr,T_pm,code}, φ::Int) where {m,L,T_llr,T_pm,code}
        @idx0 @inbounds begin
            empty!(d.tmp_potpaths_pms)
            sizehint!(d.tmp_potpaths_pms, 2*L)
            d.tmp_potpaths_cont .= false

            for ℓ in active_pis(d)
                P = get_p!(d, m, ℓ)
                pm0 = update(d.list_pms[ℓ], P[0], GF2_0, φ)
                pm1 = update(d.list_pms[ℓ], P[0], GF2_1, φ)
                d.tmp_potpaths[ℓ,0] = pm0
                d.tmp_potpaths[ℓ,1] = pm1
                push!(d.tmp_potpaths_pms, pm0)
                push!(d.tmp_potpaths_pms, pm1)
            end

            ρ = min(2*current_list_size(d), L)
            sort!(d.tmp_potpaths_pms)
            τ = d.tmp_potpaths_pms[ρ-1] # this can be optimized to O(L) according to Tal+Vardy paper

            for ℓ in active_pis(d)
                d.tmp_potpaths_cont[ℓ,0] = (d.tmp_potpaths[ℓ,0] <= τ)
                d.tmp_potpaths_cont[ℓ,1] = (d.tmp_potpaths[ℓ,1] <= τ)
            end

            # specially with discrete channel outputs / discrete LLRs / discrete PMs, tie-breaking (here randomized) might be required to have the exact desired number of surviving paths
            num_potpaths_cont = count(d.tmp_potpaths_cont)
            num_potpaths_cont_tiebreak = num_potpaths_cont - ρ
            @assert num_potpaths_cont_tiebreak >= 0
            while num_potpaths_cont_tiebreak > 0
                for ℓ in shuffle(collect(active_pis(d))), u in (0, 1)
                    if num_potpaths_cont_tiebreak == 0
                        break
                    end
                    if d.tmp_potpaths_cont[ℓ,u] && d.tmp_potpaths[ℓ,u] == τ
                        if rand() < 0.5
                            d.tmp_potpaths_cont[ℓ,u] = false
                            num_potpaths_cont_tiebreak -= 1
                        end
                    end
                end
            end
            @assert count(d.tmp_potpaths_cont) == ρ

            for ℓ in active_pis(d)
                if !d.tmp_potpaths_cont[ℓ,0] && !d.tmp_potpaths_cont[ℓ,1]
                    kill_path!(d, ℓ)
                end
            end

            # TODO: recheck to make sure active_pis(d) is not affected by newly created paths in the loop

            # for ℓ in collect(active_pis(d)) # buffer this, otherwise we also make newly created paths subject to this loop
            #   if !d.tmp_potpaths_cont[ℓ,0] && !d.tmp_potpaths_cont[ℓ,1]
            #       # was already removed, this should not happen!
            #       @assert false
            for ℓ in active_pis(d)
                if d.tmp_potpaths_cont[ℓ,0] && d.tmp_potpaths_cont[ℓ,1]
                    ℓ′ = clone_path!(d, ℓ, φ)
                    C = get_c!(d, m, ℓ)
                    C′ = get_c!(d, m, ℓ′)
                    C[0,φ%2] = GF2_0
                    C′[0,φ%2] = GF2_1
                    d.array_u[φ,ℓ] = GF2_0
                    d.array_u[φ,ℓ′] = GF2_1
                    d.list_pms[ℓ] = d.tmp_potpaths[ℓ,0]
                    d.list_pms[ℓ′] = d.tmp_potpaths[ℓ,1]
                elseif d.tmp_potpaths_cont[ℓ,0] && !d.tmp_potpaths_cont[ℓ,1]
                    C = get_c!(d, m, ℓ)
                    C[0,φ%2] = GF2_0
                    d.array_u[φ,ℓ] = GF2_0
                    d.list_pms[ℓ] = d.tmp_potpaths[ℓ,0]
                elseif !d.tmp_potpaths_cont[ℓ,0] && d.tmp_potpaths_cont[ℓ,1]
                    C = get_c!(d, m, ℓ)
                    C[0,φ%2] = GF2_1
                    d.array_u[φ,ℓ] = GF2_1
                    d.list_pms[ℓ] = d.tmp_potpaths[ℓ,1]
                else # !d.tmp_potpaths_cont[ℓ,0] && !d.tmp_potpaths_cont[ℓ,1]
                    # was already removed
                    continue
                end
            end

            return nothing
        end
    end

    function decode!(d::SCLDecoder{m,L,T_llr,T_pm,code}, llrs::Vector{T_llr}) where {m,L,T_llr,T_pm,code}
        @idx0 begin
            reset!(d)
            ℓ = assign_initial_path!(d)

            P0 = get_p!(d, 0, ℓ)
            copyto!(P0, llrs)

            for φ in 0:(n(d)-1)
                recursively_calc_p!(d, m, φ)

                if code.is_frozen[φ]
                    continue_paths_frozen!(d, φ)
                else
                    continue_paths_unfrozen!(d, φ)
                end

                if φ % 2 == 1
                    recursively_update_c!(d, m, φ)
                end
            end

            return current_list_size(d)
        end
    end

    function decode_with_debug!(d::SCLDecoder{m,L,T_llr,T_pm,code}, llrs::Vector{T_llr}, llr_logger::T_llrlogger) where {m,L,T_llr,T_pm,code,T_llrlogger<:AbstractLLRLogger}
        @idx0 begin
            reset!(d)
            ℓ = assign_initial_path!(d)

            P0 = get_p!(d, 0, ℓ)
            copyto!(P0, llrs)

            for φ in 0:(n(d)-1)
                recursively_calc_p!(d, m, φ)

                # log llr messages
                for ℓ in active_pis(d)
                    P = get_p!(d, m, ℓ)
                    log!(llr_logger, φ, P[0])
                end

                if code.is_frozen[φ]
                    continue_paths_frozen!(d, φ)
                else
                    continue_paths_unfrozen!(d, φ)
                end

                if φ % 2 == 1
                    recursively_update_c!(d, m, φ)
                end
            end

            return current_list_size(d)
        end
    end

    function reconstruct_pm!(d::SCLDecoder{m,L,T_llr,T_pm,code}, llrs::Vector{T_llr}, codeword::Vector{GF2Element}) where {m,L,T_llr,T_pm,code}
        @idx0 begin
            reset!(d)
            ℓ = assign_initial_path!(d)

            P0 = get_p!(d, 0, ℓ)
            copyto!(P0, llrs)

            for φ in 0:(n(d)-1)
                recursively_calc_p!(d, m, φ)

                P = get_p!(d, m, ℓ)
                C = get_c!(d, m, ℓ)
                C[0,φ%2] = codeword[φ]
                d.array_u[φ,ℓ] = codeword[φ]
                d.list_pms[ℓ] = update(d.list_pms[ℓ], P[0], codeword[φ], φ)

                if φ % 2 == 1
                    recursively_update_c!(d, m, φ)
                end
            end

            return d.list_pms[ℓ]
        end
    end

    function reconstruct_pm_with_debug!(d::SCLDecoder{m,L,T_llr,T_pm,code}, llrs::Vector{T_llr}, codeword::Vector{GF2Element}, llr_logger::T_llrlogger) where {m,L,T_llr,T_pm,code,T_llrlogger<:AbstractLLRLogger}
        @idx0 begin
            reset!(d)
            ℓ = assign_initial_path!(d)

            P0 = get_p!(d, 0, ℓ)
            copyto!(P0, llrs)

            for φ in 0:(n(d)-1)
                recursively_calc_p!(d, m, φ)

                # log llr message
                P = get_p!(d, m, ℓ)
                log!(llr_logger, φ, P[0])

                P = get_p!(d, m, ℓ)
                C = get_c!(d, m, ℓ)
                C[0,φ%2] = codeword[φ]
                d.array_u[φ,ℓ] = codeword[φ]
                d.list_pms[ℓ] = update(d.list_pms[ℓ], P[0], codeword[φ], φ)

                if φ % 2 == 1
                    recursively_update_c!(d, m, φ)
                end
            end

            return d.list_pms[ℓ]
        end
    end

    function extract_list_C!(d::SCLDecoder{m,L,T_llr,T_pm,code}, cs::Matrix{GF2Element}, pms::Vector{T_pm}) where {m,L,T_llr,T_pm,code}
        @idx0 begin
            i = 0
            for ℓ in active_pis(d)
                C0 = get_c!(d, 0, ℓ)
                for β in 0:(n(d)-1)
                    cs[β,i] = C0[β,0]
                end
                pms[i] = d.list_pms[ℓ]
                i += 1
            end
            return nothing
        end
    end

    function extract_list!(d::SCLDecoder{m,L,T_llr,T_pm,code}, cs::Matrix{GF2Element}, pms::Vector{T_pm}) where {m,L,T_llr,T_pm,code}
        @idx0 begin
            if current_list_size(d) == L
                copyto!(cs, d.array_u)
                copyto!(pms, d.list_pms)
            else
                i = 0
                for ℓ in active_pis(d)
                    for j in 0:(n(d)-1)
                        cs[j,i] = d.array_u[j,ℓ]
                    end
                    pms[i] = d.list_pms[ℓ]
                    i += 1
                end
            end
            return nothing
        end
    end

    function extract_best_pm_C!(d::SCLDecoder{m,L,T_llr,T_pm,code}, c::Vector{GF2Element}) where {m,L,T_llr,T_pm,code}
        @idx0 begin
            ℓ_opt = -1
            pm_opt = T_pm()
            for ℓ in active_pis(d)
                pm = d.list_pms[ℓ]
                if pm < pm_opt || ℓ_opt == -1
                    ℓ_opt = ℓ
                    pm_opt = pm
                end
            end
            @assert ℓ_opt != -1

            C0 = get_c!(d, 0, ℓ_opt)
            for β in 0:(n(d)-1)
                c[β] = C0[β,0]
            end

            return pm_opt
        end
    end

    function extract_best_pm!(d::SCLDecoder{m,L,T_llr,T_pm,code}, c::Vector{GF2Element}) where {m,L,T_llr,T_pm,code}
        @idx0 begin
            ℓ_opt = -1
            pm_opt = T_pm()
            for ℓ in active_pis(d)
                pm = d.list_pms[ℓ]
                if pm < pm_opt || ℓ_opt == -1
                    ℓ_opt = ℓ
                    pm_opt = pm
                end
            end
            @assert ℓ_opt != -1

            for i in 0:(n(d)-1)
                c[i] = d.array_u[i,ℓ_opt]
            end

            return pm_opt
        end
    end

    function extract_list_C(d::SCLDecoder{m,L,T_llr,T_pm,code}) where {m,L,T_llr,T_pm,code}
        L′ = current_list_size(d)
        cs = Matrix{GF2Element}(undef, n(d), L′)
        pms = Vector{T_pm}(undef, L′)
        extract_list_C!(d, cs, pms)
        return (cs, pms)
    end

    function extract_list(d::SCLDecoder{m,L,T_llr,T_pm,code}) where {m,L,T_llr,T_pm,code}
        L′ = current_list_size(d)
        cs = Matrix{GF2Element}(undef, n(d), L′)
        pms = Vector{T_pm}(undef, L′)
        extract_list!(d, cs, pms)
        return (cs, pms)
    end

    function extract_best_pm_C(d::SCLDecoder{m,L,T_llr,T_pm,code}) where {m,L,T_llr,T_pm,code}
        c = Vector{GF2Element}(undef, n(d))
        pm = extract_best_pm_C!(d, c)
        return (c, pm)
    end

    function extract_best_pm(d::SCLDecoder{m,L,T_llr,T_pm,code}) where {m,L,T_llr,T_pm,code}
        c = Vector{GF2Element}(undef, n(d))
        pm = extract_best_pm!(d, c)
        return (c, pm)
    end

end
