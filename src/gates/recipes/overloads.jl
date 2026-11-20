function pulse_tstops(gate::GateRecipe)
    return sort!(unique(vcat(Real[0,gate.duration],[t for c in gate.controls for t in pulse_tstops(c.pulse)])))
end
pulse_tstops(gate::PreparedGate) = pulse_tstops(gate.recipe)

"""Explicitly require controlled values at `0` and duration to match reference values."""
function check_endpoints(prepared::PreparedGate; atol=0, rtol=0)
    gate = prepared.recipe
    for c in gate.controls, t in (0,gate.duration)
        isapprox(_resolve_pulse(c.pulse,prepared.model.parameters)(t),prepared.parameters[c.target];atol,rtol) ||
            throw(ArgumentError("control $(c.name) does not match its reference value at time $t"))
    end
    return true
end

function check_endpoints(model::NumericalModel,gate::GateRecipe;atol=0,rtol=0)
    _,values = _reuse_gate_data(model,gate)
    for c in gate.controls, t in (0,gate.duration)
        haskey(values,c.target) || throw(ArgumentError("missing reference for control $(c.name)"))
        isapprox(_resolve_pulse(c.pulse,model.parameters)(t),values[c.target];atol,rtol) ||
            throw(ArgumentError("control $(c.name) does not match its reference value at time $t"))
    end
    return true
end
