struct QDDimension
    dim::Union{Float64, Integer}

    function QDDimension(x::Union{Float64,Integer})
        valid = (x isa Integer && x > 0) || (x isa Float64 && isequal(x, Inf))
        valid || throw(ArgumentError("dimension must be a positive integer or Inf"))
        return new(x)
    end
end
