__precompile__()

module Stacks

    export Stack
    export FixedCapacityStack
    export isfull

    import Base: isempty, empty!, length
    import Base: push!, pop!

    abstract type Stack end

    mutable struct FixedCapacityStack{T,N} <: Stack
        elements::Vector{T}
        index::UInt

        function FixedCapacityStack{T,N}() where {T,N}
            @assert isa(N, Int)
            s = new()
            s.elements = Vector{T}(undef, N)
            s.index = 0
            return s
        end
    end

    @inline isempty(s::FixedCapacityStack) = (s.index == 0)
    @inline isfull(s::FixedCapacityStack{T,N}) where {T,N} = (s.index == N)
    @inline length(s::FixedCapacityStack) = s.index
    # @inline endof(s::FixedCapacityStack) = length(s)

    @inline function empty!(s::FixedCapacityStack)
        s.index = 0
        return
    end

    @inline function push!(s::FixedCapacityStack{T,N}, e::T) where {T,N}
        @boundscheck if isfull(s)
            throw(BoundsError())
        end
        s.index += 1
        @inbounds s.elements[s.index] = e
        return
    end

    @inline function pop!(s::FixedCapacityStack)
        @boundscheck if isempty(s)
            throw(BoundsError())
        end
        @inbounds e = s.elements[s.index]
        s.index -= 1
        return e
    end

end
