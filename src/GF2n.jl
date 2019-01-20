__precompile__()

module GF2n

    export GF2Element, GF2nElement, FiniteFieldElement
    export GF2_0, GF2_1
    export toint, elements
    export dim, generator

    import Base: +, -, *, /, ⊻
    import Base: zero, one, iszero
    import Base: rand, length
    import Base: show


    # abstract type for elements of fields with characteristic two: GF2nElement

    abstract type FiniteFieldElement <: Number end
    abstract type GF2nElement <: FiniteFieldElement end

    @inline zero(x::T) where {T<:GF2nElement} = T(0)
    @inline zero(::Type{T}) where {T<:GF2nElement} = T(0)
    @inline one(x::T) where {T<:GF2nElement} = T(1)
    @inline one(::Type{T}) where {T<:GF2nElement} = T(1)
    @inline iszero(x::T) where {T<:GF2nElement} = (x.val == 0)

    rand(::Type{T}) where {T<:GF2nElement} = T(rand(Int))
    rand(::Type{T}, dims::Integer...) where {T<:GF2nElement} = T.(rand(Int, dims...))
    # countnz(x::Array{T,N}) where {T<:GF2nElement,N} = count(!iszero, x)
    length(::Type{T}) where {T<:GF2nElement} = 2^dim(T)

    toint(x::T) where {T<:GF2nElement} = Int(x.val)
    elements(::Type{T}) where {T<:GF2nElement} = [ T(i) for i in 0:((2^dim(T))-1) ]

    +(x::T, y::T) where {T<:GF2nElement} = T(x.val ⊻ y.val)
    -(x::T, y::T) where {T<:GF2nElement} = T(x.val ⊻ y.val)
    -(x::T) where {T<:GF2nElement} = x
    function /(x::T, y::T) where {T<:GF2nElement}
        if iszero(y)
            throw(DivideError())
        end
        return gfdiv_unsafe(x, y)
    end


    # type for elements of two-element finite field: GF2Element

    struct GF2Element <: GF2nElement
        val::UInt8

        function GF2Element(x::T) where {T<:Integer}
            return new(x & 0x01)
        end
    end

    dim(::Type{GF2Element}) = 1
    generator(::Type{GF2Element}) = GF2Element(0x01)

    *(x::GF2Element, y::GF2Element) = GF2Element(x.val & y.val)
    function gfdiv_unsafe(x::GF2Element, y::GF2Element)
        return x
    end

    function show(io::IO, x::GF2Element)
        print(io, (x.val == 0 ? "0" : "1"), "₂")
    end


    # some helpers

    const GF2_0 = GF2Element(0)
    const GF2_1 = GF2Element(1)
    ⊻(x::GF2Element, y::GF2Element) = x + y

end
