# Domain types and evaluation live together; explicit overloads are colocated.
function _require_finite(name, value)
    value isa Number && isfinite(value) || throw(ArgumentError("$name must be a finite number"))
    return value
end
function _require_finite_real(name, value)
    value isa Real && isfinite(value) || throw(ArgumentError("$name must be a finite real number"))
    return value
end

include("pulses/envelopes.jl")
include("pulses/carriers.jl")
include("pulses/pulse.jl")
include("pulses/overloads.jl")
