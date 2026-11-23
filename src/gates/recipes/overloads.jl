function pulse_tstops(gate::GateRecipe)
    return sort!(unique(vcat(Real[0,gate.duration],[t for c in gate.controls for t in pulse_tstops(c.pulse)])))
end
pulse_tstops(gate::PreparedGate) = pulse_tstops(gate.recipe)

"""Explicitly require controlled values at `0` and duration to match reference values."""
function _check_control_endpoints(gate,values;atol,rtol)
    bindings = _compile_gate_bindings(gate,values,())
    for (target,waveform) in bindings, t in (0,gate.duration)
        haskey(values,target) || throw(ArgumentError("missing physical reference for control target $target"))
        isapprox(waveform(nothing,t),values[target];atol,rtol) ||
            throw(ArgumentError("control target $target does not match its reference value at time $t"))
    end
    return true
end
function check_endpoints(prepared::PreparedGate;atol=0,rtol=0)
    return _check_control_endpoints(prepared.recipe,prepared.parameters;atol,rtol)
end
function check_endpoints(model::NumericalModel,gate::GateRecipe;atol=0,rtol=0)
    _,values = _reuse_gate_data(model,gate)
    return _check_control_endpoints(gate,values;atol,rtol)
end
