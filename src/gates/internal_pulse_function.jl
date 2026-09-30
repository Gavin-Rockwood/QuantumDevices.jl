const INTERNAL_PULSE_GENERATORS = Dict{Symbol,Function}()

function _register_internal_pulse!(name::Symbol, generator::Function)
    haskey(INTERNAL_PULSE_GENERATORS, name) &&
        throw(ArgumentError("Internal pulse is already registered: $name"))
    INTERNAL_PULSE_GENERATORS[name] = generator
    return generator
end

function _internal_pulse_generator(name::Symbol)
    generator = get(INTERNAL_PULSE_GENERATORS, name, nothing)
    generator === nothing && throw(ArgumentError("Unknown internal pulse: $name"))
    return generator
end

"""
    InternalPulseFunction(name::Symbol, parameters::NamedTuple)

Parameterized pulse referencing an existing registered `f(p, t)` generator.
Prefer public pulse constructors, which supply the expected fields. Unknown names
and stored `duration` fields throw `ArgumentError`. Its name and parameters are
stored in JSON without callable artifacts.
"""
struct InternalPulseFunction{P<:NamedTuple} <: AbstractParameterizedPulse
    name::Symbol
    parameters::P

    function InternalPulseFunction(name::Symbol, parameters::P) where {P<:NamedTuple}
        _check_pulse_parameters(parameters)
        _internal_pulse_generator(name)
        new{P}(name, parameters)
    end
end

pulse_function(pulse::InternalPulseFunction) = _internal_pulse_generator(pulse.name)

function _pulse_window(start, stop, duration)
    duration isa Real && isfinite(duration) && duration >= 0 ||
        throw(ArgumentError("Gate duration must be finite and nonnegative"))
    _require_finite_real("start", start)
    stop = stop === nothing ? duration : stop
    _require_finite_real("stop", stop)
    0 <= start <= stop <= duration ||
        throw(ArgumentError("Pulse window must satisfy 0 <= start <= stop <= gate duration"))
    return start, stop
end
