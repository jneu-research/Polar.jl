__precompile__()

module Distributions

    # infrastructure to make sure we're always operating on proper distributions with unique variables

    export isdistribution
    isdistribution(d::Dict{T,Float64}) where {T} = sum(values(d)) ≈ 1.0 && all(x -> 0 <= x, values(d))
    hasuniqueelements(c) = length(unique(c)) == length(c)


    # the distribution type

    export Distribution

    struct Distribution{col_types,N}
        col_names :: NTuple{N,String}
        data :: Dict{col_types,Float64}

        function Distribution(col_names::NTuple{N,String}, data::Dict{col_types,Float64}) where {col_types<:Tuple,N}
            @assert nfields(col_types) == N
            @assert isdistribution(data)
            @assert hasuniqueelements(col_names)

            return new{col_types,N}(col_names, data)
        end
    end

    function _prefixvars(d::Distribution, prefix::String)
        col_names = Tuple([ "$(prefix)_$(var)" for var in variables(d) ])
        return rename(d, col_names)
    end

    function _raw_transform(f::F, d_in::Dict{T_k,T_v}) where {F,T_k,T_v}
        result_vec = map(
            (k, v) -> (f(k), v),
            keys(d_in),
            values(d_in)
        )
        return _raw_vec2dict(result_vec)
    end

    function _raw_transform(f::F, d_in1::Dict{T_k1,T_v1}, d_in2::Dict{T_k2,T_v2}) where {F,T_k1,T_v1,T_k2,T_v2}
        result_vec = map(
            (k, v) -> (f(k...), prod(v)),
            Iterators.product(keys(d_in1), keys(d_in2)),
            Iterators.product(values(d_in1), values(d_in2))
        )
        return _raw_vec2dict(result_vec)
    end

    function _raw_transform(f::F, d_in1::Dict{T_k1,T_v1}, d_in2::Dict{T_k2,T_v2}, d_in3::Dict{T_k3,T_v3}) where {F,T_k1,T_v1,T_k2,T_v2,T_k3,T_v3}
        result_vec = map(
            (k, v) -> (f(k...), prod(v)),
            Iterators.product(keys(d_in1), keys(d_in2), keys(d_in3)),
            Iterators.product(values(d_in1), values(d_in2), values(d_in3))
        )
        return _raw_vec2dict(result_vec)
    end

    function _raw_vec2dict(result_vec::Array{Tuple{T_k,T_v},N}) where {T_k,T_v,N}
        result_dict = Dict{T_k,T_v}()
        for (k, v) in result_vec
            result_dict[k] = get(result_dict, k, 0.0) + v
        end
        return result_dict
    end


    export support, variables, rename

    support(d::Distribution) = keys(d.data)
    variables(d::Distribution) = d.col_names
    rename(d::Distribution{T,N}, col_names::NTuple{N,String}) where {T,N} = Distribution(col_names, d.data)


    export ×
    import Base: getindex

    getindex(d::Distribution{T,N}, val::T) where {T,N} = get(d.data, val, 0.0)

    function ×(d_in1::Distribution{col_types1,N1}, d_in2::Distribution{col_types2,N2}) where {col_types1,col_types2,N1,N2}
        col_names_out = (variables(d_in1)..., variables(d_in2)...)
        f = (t1, t2) -> (t1..., t2...)
        return Distribution(col_names_out, _raw_transform(f, d_in1.data, d_in2.data))
    end


    export transform

    function transform(f::F, col_names_out::NTuple{N,String}, d_in1::Distribution) where {F,N}
        d_in1_ = _prefixvars(d_in1, "1")
        return Distribution(col_names_out, _raw_transform(f, d_in1_.data))
    end

    function transform(f::F, col_names_out::NTuple{N,String}, d_in1::Distribution, d_in2::Distribution) where {F,N}
        d_in1_ = _prefixvars(d_in1, "1")
        d_in2_ = _prefixvars(d_in2, "2")
        return Distribution(col_names_out, _raw_transform(f, d_in1_.data, d_in2_.data))
    end

    function transform(f::F, col_names_out::NTuple{N,String}, d_in1::Distribution, d_in2::Distribution, d_in3::Distribution) where {F,N}
        d_in1_ = _prefixvars(d_in1, "1")
        d_in2_ = _prefixvars(d_in2, "2")
        d_in3_ = _prefixvars(d_in3, "3")
        return Distribution(col_names_out, _raw_transform(f, d_in1_.data, d_in2_.data, d_in3_.data))
    end


    export marginalize   # careful, not type-stable (cannot be stabilized!?)

    function marginalize(d_in::Distribution{T,N}, cols::NTuple{M,String}) where {T,N,M}
        @assert all(c -> c ∈ variables(d_in), cols)
        col_names_out = collect(Iterators.filter(v -> v ∉ cols, variables(d_in)))
        idxs_keep = findin(variables(d_in), col_names_out) :: Vector{Int}
        f = t -> t[idxs_keep]
        return transform(f, Tuple(col_names_out), d_in)
    end


    export expectedvalue

    function expectedvalue(d::Distribution{T,N}) where {T,N}
        add_tpl = (x, y) -> broadcast(+, x, y)
        return reduce(add_tpl, [ val.*prob for (val, prob) in d.data ])
    end


    export conditionalize

    function conditionalize(d_in::Distribution{T,N}, cols::NTuple{M,String}, vals) where {T,N,M}
        @assert all(c -> c ∈ variables(d_in), cols)
        @assert length(cols) == length(vals)
        col_names_out = collect(Iterators.filter(v -> v ∉ cols, variables(d_in)))
        idxs_keep = findin(variables(d_in), col_names_out) :: Vector{Int}
        idxs_remove = findin(variables(d_in), cols) :: Vector{Int}
        f = t -> (t[idxs_keep], t[idxs_remove] == vals)
        return Distribution(Tuple(col_names_out), _raw_conditionalize_rescaleremap(_raw_conditionalize_filter(_raw_transform(f, d_in.data))))
    end

    function _raw_conditionalize_filter(data::Dict{Tuple{T,Bool},Float64}) where {T}
        return filter((k, v) -> k[2], data)
    end

    function _raw_conditionalize_rescaleremap(data::Dict{Tuple{T,Bool},Float64}) where {T}
        C = sum([ v for (k, v) in data ]) :: Float64
        result_vec = map((k, v) -> (k[1], v/C :: Float64), keys(data), values(data))
        return _raw_vec2dict(result_vec)
    end

end
