function Base.getproperty(gate::DeviceGate, key::Symbol)
    if key === :duration
        explicit = getfield(gate, :duration_override)
        explicit !== nothing && return explicit
        return maximum(p.delay + p.duration for p in values(gate.parameters) if p isa AbstractPulse)
    end
    return getfield(gate, key)
end
Base.propertynames(::DeviceGate, private::Bool=false) = private ?
    (:parameters, :hamiltonian, :duration, :duration_override) : (:parameters, :hamiltonian, :duration)
function _replace_property(gate::DeviceGate, key::Symbol, value)
    key === :parameters && return DeviceGate(value, gate.hamiltonian, gate.duration_override)
    key === :hamiltonian && return DeviceGate(gate.parameters, value, gate.duration_override)
    key === :duration && return DeviceGate(gate.parameters, gate.hamiltonian, value)
    throw(ArgumentError("DeviceGate.$key is not a writable public field"))
end

_check_gate_endpoints(value, idle, duration) = value == idle
_check_gate_endpoints(pulse::AbstractPulse, idle, duration) = pulse(0.0) == idle && pulse(duration) == idle

"""
    numerical(model::DeviceModel, gate::DeviceGate)

Combine the projected idle Hamiltonian and gate drive with gate parameter overrides.
Overrides of idle parameters must equal their idle values at both endpoints.
"""
function numerical(model::DeviceModel, gate::DeviceGate)
    combined = Dict{Symbol,Any}(pairs(model.parameters))
    for (key, value) in pairs(gate.parameters)
        if haskey(combined, key) && !_check_gate_endpoints(value, combined[key], gate.duration)
            throw(ArgumentError("Endpoints of parameter $key do not agree with the idle parameter"))
        end
        combined[key] = value isa AbstractPulse ?
            (t -> _require_finite("Pulse value", value(t))) : value
    end
    required = union(parameter_keys(model.hamiltonian), parameter_keys(gate.hamiltonian))
    missing_ = setdiff(required, Set(keys(combined)))
    isempty(missing_) || throw(ArgumentError("Missing gate parameters: $missing_"))
    dims = [get(model.truncation_dimensions, c, c.dimension) for c in model.components]
    drive_terms = _projected_terms(gate.hamiltonian, model.components, dims)
    return _evaluate_terms(vcat(model.compiled_terms, drive_terms), combined)
end


function parameters(gate::DeviceGate)
    result = Dict{String,Tuple}()
    for (name, value) in pairs(gate.parameters)
        key = _parameter_name(name)
        path = hasproperty(gate, name) ? "parameters/$key" : key
        if value isa AbstractPulse
            for (nested_name, (nested_path, leaf)) in parameters(value)
                result["$key/$nested_name"] = ("$path/$nested_path", leaf)
            end
        else
            result[key] = (path, value)
        end
    end
    return result
end

SciMLBase.solve(problem::CalibrationProblem, args...; kwargs...) =
    SciMLBase.solve(problem.problem, args...; kwargs...)

function pulse_tstops(gate::DeviceGate)
    stops = Real[]
    for control in values(gate.parameters)
        control isa AbstractPulse && append!(stops, pulse_tstops(control))
    end
    filter!(t -> 0 < t < gate.duration, stops)
    sort!(unique!(stops))
    return isempty(stops) ? Float64[] : collect(promote(stops...))
end

function get_unitary(model::DeviceModel, gate::DeviceGate; kwargs...)
    options = merge((; tstops=pulse_tstops(gate)), (; kwargs...))
    return get_unitary(numerical(model, gate), gate.duration; options...)
end
