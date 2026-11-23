using Interpolations
using SpecialFunctions: erf
function _require_finite(name, value)
    concrete_parameter_value(value) || throw(ArgumentError("$name must be a finite scalar"))
    return value
end
function _require_finite_real(name, value)
    value isa Real || throw(ArgumentError("$name must be real"))
    return _require_finite(name, value)
end
include("envelopes/envelopes.jl")
include("carriers.jl")
include("recipes.jl")
include("pulse.jl")
include("PulseSequence.jl")
include("overloads.jl")
include("parameters.jl")
include("envelopes/compilation.jl")
include("compilation.jl")
include("runtime_parameters.jl")
include("timing.jl")
include("display.jl")
