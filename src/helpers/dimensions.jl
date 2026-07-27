struct Dimension
    size::Tuple{Vararg{Union{Int, Float64}}}

    function Dimension(values::Tuple)
        size = map(values) do value
            if value isa Integer && !(value isa Bool) && value > 0
                Int(value)
            elseif value isa AbstractFloat && isinf(value) && value > 0
                Inf
            else
                throw(ArgumentError("each dimension must be a positive integer or Inf; got $value",))
            end
        end
        new(Tuple(size))
    end
end

function Dimension(val :: Real)
    Dimension((val,))
end

function Dimension(X::Dimension...)
    sizes = [x.size for x in X]
    Dimension(Tuple(Iterators.flatten(sizes)))
end


Base.iterate(dimension::Dimension, state...) =
    iterate(dimension.size, state...)

Base.length(dimension::Dimension) = length(dimension.size)
Base.eltype(::Type{Dimension}) = Union{Int, Float64}

Base.getindex(dimension::Dimension, index...) =
    getindex(dimension.size, index...)

Base.IteratorSize(::Type{Dimension}) = Base.HasLength()