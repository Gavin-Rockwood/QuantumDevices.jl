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
    numerical(model::DeviceModel, gate::DeviceGate; scalar=1, dense=false)

Combine the projected idle Hamiltonian and gate drive with gate parameter overrides.
Overrides of idle parameters must equal their idle values at both endpoints.
`dense=true` materializes all operator matrices as dense matrices at construction.
`scalar` is folded into each operator matrix and constant coefficient before
constructing the Hamiltonian. Use `scalar=2pi` when passing it to native `sesolve`.
"""
function numerical(model::DeviceModel, gate::DeviceGate; scalar::Number=1, dense::Bool=false)
    combined = _gate_parameters(model, gate)
    dims = [get(model.truncation_dimensions, c, c.dimension) for c in model.components]
    drive_terms = _projected_terms(gate.hamiltonian, model.components, dims)
    return _evaluate_terms(vcat(model.compiled_terms, drive_terms), combined; scalar, dense)
end

function _gate_parameters(model, gate)
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
    return combined
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
    return _get_unitary_angular(numerical(model, gate; scalar=2pi), gate.duration; options...)
end

"""
    numerical(model::DeviceModel, gate::DeviceGate, states;
              output_states=states, frame=:lab, energy_shift=:auto, dense=false, kwargs...)

Return a dimensionless `QuantumObject` operator in the supplied ordered state basis.
Its dimensions are the number of selected states; `.data` gives the gate matrix.
`states` and `output_states` are vectors of orthonormal QuantumToolbox kets.
Evolve the input columns through the full retained model with `2π * H`, then
project onto `output_states`. Columns are never renormalized, retaining leakage.
This does not project the Hamiltonian before evolution.

`frame=:lab` preserves lab-frame phases; `frame=:interaction` removes idle
evolution before the final projection. Native solver keywords and pulse event
times are handled as in [`get_gate_matrix`](@ref). Automatic energy centering
uses the mean input-state idle energy and restores global phase.
Set `energy_shift=0` to disable it, or supply an offset in frequency units.
`dense=true` uses dense Hamiltonian matrices during evolution; the returned
projected operator already has dense matrix storage.

The two-argument `numerical(model, gate)` returns the Hamiltonian in frequency
units; this three-argument form returns its evolved action on the selected states.
"""
function numerical(model::DeviceModel, gate::DeviceGate, states;
    output_states=states, frame=:lab, energy_shift=:auto, dense::Bool=false, kwargs...)
    initial, output = _gate_bases(states, output_states)
    projection = _prepare_gate_projection(output, frame, model.H)
    options = merge((; tstops=pulse_tstops(gate)), (; kwargs...))
    shift = _energy_shift(model.H, initial, energy_shift)
    evaluate = _prepared_gate_numerical(model, gate; scalar=2pi, energy_shift=shift, dense)
    columns = _evolve_gate_columns(evaluate(gate), initial, gate.duration;
        phase_shift=2pi*shift, options...)
    return QuantumObject(_project_gate_columns(columns, projection, gate.duration); type=Operator())
end

get_gate_matrix(model::DeviceModel, gate::DeviceGate, states; kwargs...) =
    numerical(model, gate, states; kwargs...)

function Base.show(io::IO, gate::DeviceGate)
    controls = count(value -> value isa AbstractPulse, values(gate.parameters))
    print(io, "DeviceGate(duration=", gate.duration, ", controls=", controls,
        ", scalars=", length(gate.parameters) - controls, ")")
end

function Base.show(io::IO, ::MIME"text/plain", gate::DeviceGate)
    get(io, :compact, false) && return show(io, gate)
    show(io, gate)
    print(io, "\n  Duration: ", gate.duration,
        gate.duration_override === nothing ? " (inferred)" : " (explicit)")
    print(io, "\n  Hamiltonian: ")
    limited = IOContext(io, :compact => true, :limit => true)
    show(limited, gate.hamiltonian)
    isempty(gate.parameters) && print(io, "\n  Parameters: none")
    for (name, value) in pairs(gate.parameters)
        print(io, "\n  ", name, ": ")
        if value isa Pulse
            print(io, "Pulse\n")
            _show_pulse_details(io, value; end_time=gate.duration, indent="    ")
        else
            show(limited, value)
        end
    end
end
