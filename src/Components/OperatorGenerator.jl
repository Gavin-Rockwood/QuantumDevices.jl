"""A numerical operator factory with a positive default state count."""
struct OperatorGenerator
    f::Any
    default_dimension::Int
    function OperatorGenerator(f, dimension::Integer)
        dimension > 0 || throw(ArgumentError("default dimension must be positive"))
        return new(f, Int(dimension))
    end
end
