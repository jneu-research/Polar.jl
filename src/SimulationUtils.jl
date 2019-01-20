__precompile__()

module SimulationUtils

    using Distributed
    using Random

    using Polar.PolarCodes
    using Polar.PolarDecoding
    using Polar.Channels
    using Polar.GF2n
    using Polar.CommunicationsUtils

    import Polar.Utils: gnuplotexport


    # confidence intervals

    export confidence_wilson, confidence_wald

    function confidence_wilson(N_hit::Int, N_all::Int)   # better, more precise, smaller!
        z = 1.96
        fer = N_hit / N_all
        return (z / (1+z^2/N_all) * sqrt(fer*(1-fer)/N_all + z^2/(4*N_all^2)))
    end

    function confidence_wald(N_hit::Int, N_all::Int)   # worse, more loose, larger!
        z = 1.96
        fer = N_hit / N_all
        return (z * sqrt(fer*(1-fer)/N_all))
    end


    # simulation result datatype

    export SimulationResult

    struct SimulationResult{events}
        n_events::Dict{Symbol,Int}

        function SimulationResult{events}(n_events::Dict{Symbol,Int}) where events
            events_ = Tuple(sort(collect(events)))
            n_events_ = copy(n_events)
            @assert all(k -> k ∈ events_ || k == :ANY, keys(n_events_))
            for e in events_
                if e ∉ keys(n_events_)
                    n_events_[e] = 0
                end
            end
            return new{events_}(n_events_)
        end
    end

    function SimulationResult(events::Symbol...)
        n_events = Dict{Symbol,Int}(:ANY => 0)
        for e in events
            n_events[e] = 0
        end
        return SimulationResult{events}(n_events)
    end

    function SimulationResult(n_events::Dict{Symbol,Int})
        n_events_ = copy(n_events)
        @assert :ANY in keys(n_events)
        return SimulationResult{Tuple(Iterators.filter(k -> k != :ANY, keys(n_events_)))}(n_events_)
    end

    # events(r::SimulationResult{evs}) where evs = evs


    export rate

    function rate(r::SimulationResult, e::Symbol)
        return r[e]/r[:ANY]
    end


    export ci_wilson, ci_wald, ci_wilson_rel, ci_wald_rel

    function ci_wilson(r::SimulationResult, e::Symbol)
        if r[:ANY] == 0
            return Inf
        end

        return confidence_wilson(r[e], r[:ANY])
    end

    function ci_wald(r::SimulationResult, e::Symbol)
        if r[:ANY] == 0
            return Inf
        end

        return confidence_wald(r[e], r[:ANY])
    end

    function ci_wilson_rel(r::SimulationResult, e::Symbol)
        if r[:ANY] == 0
            return Inf
        end

        return ci_wilson(r, e) / rate(r, e)
    end

    function ci_wald_rel(r::SimulationResult, e::Symbol)
        if r[:ANY] == 0
            return Inf
        end

        return ci_wald(r, e) / rate(r, e)
    end


    import Base: getindex

    getindex(r::SimulationResult, e::Symbol) = r.n_events[e]


    import Base: +

    function +(r1::SimulationResult{evs}, r2::SimulationResult{evs}) where evs
        n_events = Dict{Symbol,Int}(:ANY => r1[:ANY] + r2[:ANY])
        for e in evs
            n_events[e] = r1[e] + r2[e]
        end
        return SimulationResult{evs}(n_events)
    end


    export tic!

    function tic!(r::SimulationResult)
        return tic!(r, :ANY)
    end

    function tic!(r::SimulationResult{evs}, e::Symbol) where evs
        @assert e in (:ANY, evs...)
        r.n_events[e] += 1
        return r
    end


    function gnuplotexport(filename, xlabel, data::Dict{T_key,SimulationResult{T_fields}}; lt=<=, xmap=identity, ymap=identity) where {T_key,T_fields}
        ylabels = "ANY"
        for e in T_fields
            e_string = repr(e)[2:end]
            ylabels = "$(ylabels) $(e_string) $(e_string)_rate $(e_string)_ci_wilson $(e_string)_ci_wilson_rel"
        end

        open(filename, "w") do fp
            println(fp, "$(xlabel) $(ylabels)")
            K = collect(keys(data))
            sort!(K; lt=lt)
            for k in K
                v = ymap(data[k])
                value = "$(v[:ANY])"
                for e in T_fields
                    value = "$(value) $(v[e]) $(rate(v, e)) $(ci_wilson(v, e)) $(ci_wilson_rel(v, e))"
                end
                println(fp, "$(xmap(k)) $(value)")
            end
        end
    end


    # some helper functions, not meant to be used outside this module

    function sim_helper_prepare_tmp_variables_01(sim_code, sim_channel, ::Type{sim_T_llr}) where {sim_T_llr<:LLRType}
        k = PolarCodes.k(sim_code)
        n = PolarCodes.n(sim_code)

        tmp_u = Vector{GF2Element}(undef, k)
        tmp_u_full = Vector{GF2Element}(undef, n)
        tmp_c = Vector{GF2Element}(undef, n)
        tmp_y = Vector{outtype(typeof(sim_channel))}(undef, n)
        tmp_llrs = Vector{sim_T_llr}(undef, n)
        tmp_ĉ = Vector{GF2Element}(undef, n)
        tmp_û = Vector{GF2Element}(undef, k)

        return (tmp_u, tmp_u_full, tmp_c, tmp_y, tmp_llrs, tmp_ĉ, tmp_û)
    end

    function sim_helper_prepare_channel_functions_01(sim_channel, ::Type{sim_T_llr}) where {sim_T_llr<:LLRType}
        ch_f = get_f(sim_channel)
        ch_llr = get_llr(sim_channel)
        ch_y2llr = y -> sim_T_llr(ch_llr(y))

        return (ch_f, ch_y2llr)
    end


    # simulation experiment: simple experiment with PM error, list error and ML-LB

    export sim_parallel_experiments_fer_simple01, sim_local_experiments_fer_simple01

    function sim_local_experiments_fer_simple01(sim_code, sim_decoder, sim_channel; N_iterations=100)
        sim_T_llr = llrtype(sim_decoder)
        sim_T_pm = pmtype(sim_decoder)

        (tmp_u, tmp_u_full, tmp_c, tmp_y, tmp_llrs, tmp_ĉ, tmp_û) = sim_helper_prepare_tmp_variables_01(sim_code, sim_channel, sim_T_llr)
        (ch_f, ch_y2llr) = sim_helper_prepare_channel_functions_01(sim_channel, sim_T_llr)

        result = SimulationResult(:PMError, :ListError, :MLError)

        for i_iteration in 1:N_iterations

            # count a new :ANY
            tic!(result)

            # generate dummy data and polar-encode
            rand!(sim_code, tmp_u)
            polarencode!(sim_code, tmp_c, tmp_u)

            # simulate channel
            map!(ch_f, tmp_y, tmp_c)
            map!(ch_y2llr, tmp_llrs, tmp_y)

            # decode (SCL)
            decode!(sim_decoder, tmp_llrs)
            pm_opt = extract_best_pm!(sim_decoder, tmp_ĉ)

            # deinterleave = polar-decode
            deinterleave!(sim_code, tmp_û, tmp_ĉ)

            # count frame errors
            if tmp_u != tmp_û

                # PM error: satisfied by "if" already!
                tic!(result, :PMError)

                # preparatory steps for list-FER and ML-LB
                (list_cs, list_pms) = extract_list(sim_decoder)
                interleave!(sim_code, tmp_u_full, tmp_u)
                pm_ml = reconstruct_pm!(sim_decoder, tmp_llrs, tmp_u_full)

                # list error: is any list item the true data?
                if any(l -> (list_cs[:,l] == tmp_u_full), 1:length(list_pms))
                    # yes -> no list error

                    # ML error
                    if pm_opt < pm_ml
                        # the PM of the decoder output is better than
                        # the PM of the true transmitted codeword -> ML error
                        tic!(result, :MLError)
                    elseif pm_opt == pm_ml
                        # the list does contain the true data, don't add it
                        # if rand() > 1. / (count(list_pms .== pm_ml) + 0)
                        if rand() > 1. / (count(x -> x == pm_ml, list_pms) + 0)
                            tic!(result, :MLError)
                        end
                    end
                else
                    # no -> list error
                    tic!(result, :ListError)

                    # ML error
                    if pm_opt < pm_ml
                        # the PM of the decoder output is better than
                        # the PM of the true transmitted codeword -> ML error
                        tic!(result, :MLError)
                    elseif pm_opt == pm_ml
                        # the list does not contain the true data, add it
                        # if rand() > 1. / (count(list_pms .== pm_ml) + 1)
                        if rand() > 1. / (count(x -> x == pm_ml, list_pms) + 1)
                            tic!(result, :MLError)
                        end
                    end
                end

            end

        end

        return result
    end

    function sim_parallel_experiments_fer_simple01(sim_code, sim_decoder, sim_channel; desired_ci_rel=0.1, desired_ci_event=:PMError, N_jobsperbatch=96, N_framesperjob=200)
        result = SimulationResult(:PMError, :ListError, :MLError)

        @timev while ci_wilson_rel(result, desired_ci_event) > desired_ci_rel
            result = result + reduce(+, pmap(i -> sim_local_experiments_fer_simple01(sim_code, sim_decoder, sim_channel; N_iterations=N_framesperjob), 1:N_jobsperbatch))
        end

        return result
    end


    # simulation experiment: experiment with ML-in-list

    export sim_parallel_experiments_fer_mlinlist01, sim_local_experiments_fer_mlinlist01

    function sim_local_experiments_fer_mlinlist01(sim_code, sim_decoder, sim_channel; N_iterations=100)
        sim_T_llr = llrtype(sim_decoder)
        sim_T_pm = pmtype(sim_decoder)

        (tmp_u, tmp_u_full, tmp_c, tmp_y, tmp_llrs, tmp_ĉ, tmp_û) = sim_helper_prepare_tmp_variables_01(sim_code, sim_channel, sim_T_llr)
        (ch_f, ch_y2llr) = sim_helper_prepare_channel_functions_01(sim_channel, sim_T_llr)

        # channel function custom to ML-in-list decoding
        ch_llr = get_llr(sim_channel)

        # # variables to set up ML-in-list decoding variables
        # const n = MySCLDecoder.n(sim_decoder)
        # const L = MySCLDecoder.L(sim_decoder)

        # # custom variables for ML-in-list decoding
        # const tmp_cs = Matrix{GF2Element}(undef, n, L)
        # const tmp_pms = Vector{sim_T_pm}(undef, L)
        # const tmp_float_llrs = Vector{Float64}(undef, n)
        # const tmp_likelihoods = Vector{Float64}(undef, L)

        # variables to set up ML-in-list decoding variables
        n = MySCLDecoder.n(sim_decoder)
        L = MySCLDecoder.L(sim_decoder)

        # custom variables for ML-in-list decoding
        tmp_cs = Matrix{GF2Element}(undef, n, L)
        tmp_pms = Vector{sim_T_pm}(undef, L)
        tmp_float_llrs = Vector{Float64}(undef, n)
        tmp_likelihoods = Vector{Float64}(undef, L)

        result = SimulationResult(:MLInListError)

        for i_iteration in 1:N_iterations

            # count a new :ANY
            tic!(result)

            # generate dummy data and polar-encode
            rand!(sim_code, tmp_u)
            polarencode!(sim_code, tmp_c, tmp_u)

            # simulate channel
            map!(ch_f, tmp_y, tmp_c)
            map!(ch_y2llr, tmp_llrs, tmp_y)

            # decode (SCL)
            L_used = decode!(sim_decoder, tmp_llrs)
            extract_list_C!(sim_decoder, tmp_cs, tmp_pms)
            @assert L_used == L

            # calculate likelihoods of list entries
            map!(ch_llr, tmp_float_llrs, tmp_y)
            for l in 1:L
                tmp_likelihoods[l] = 0.0
                for i in 1:n
                    tmp_likelihoods[l] += ((-1)^toint(tmp_cs[i,l]))*tmp_float_llrs[i]
                end
            end
            max_val, max_idx = findmax(tmp_likelihoods)
            copyto!(tmp_ĉ, tmp_cs[:,max_idx])

            # polar-decode (deinterleave is not enough, because we are working on codewords: extract_list_C!)
            polardecode!(sim_code, tmp_û, tmp_ĉ)

            # count ML-in-list frame errors
            if any(tmp_u .!= tmp_û)
                tic!(result, :MLInListError)
            end

        end

        return result
    end

    function sim_parallel_experiments_fer_mlinlist01(sim_code, sim_decoder, sim_channel; desired_ci_rel=0.1, desired_ci_event=:MLInListError, N_jobsperbatch=96, N_framesperjob=200)
        result = SimulationResult(:MLInListError)

        @timev while ci_wilson_rel(result, desired_ci_event) > desired_ci_rel
            result = result + reduce(+, pmap(i -> sim_local_experiments_fer_mlinlist01(sim_code, sim_decoder, sim_channel; N_iterations=N_framesperjob), 1:N_jobsperbatch))
        end

        return result
    end

end
