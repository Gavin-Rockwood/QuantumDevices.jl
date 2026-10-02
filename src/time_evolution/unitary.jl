"""
    get_unitary(H, duration; kwargs...)
    get_unitary(model::DeviceModel, gate::DeviceGate; kwargs...)

Evolve the identity over `[0, duration]` and return the final retained-space
propagator. `H` is a static quantum operator or `QobjEvo` in frequency units
(cycles per unit time), e.g. GHz for time in ns. Evolution explicitly uses
`2π * H` in QuantumToolbox `sesolve`. Duration must be finite and positive.

Solver kwargs pass through, with `progress_bar=false` by default. Unsuccessful
retcodes raise an error. The model/gate overload supplies `pulse_tstops(gate)`
by default so adaptive stepping resolves delayed controls. No rotating frame or rotating-wave approximation is applied.
"""
function get_unitary(H::AbstractQuantumObject, duration::Real; kwargs...)
    isfinite(duration) && duration > 0 ||
        throw(ArgumentError("Unitary evolution requires a finite positive duration"))
    initial = qeye_like(H isa QobjEvo ? H(0.0) : H)
    options = merge((; progress_bar=false), (; kwargs...))
    solution = sesolve(2pi * H, initial, [zero(duration), duration]; options...)
    SciMLBase.successful_retcode(solution.retcode) ||
        error("Unitary evolution failed with return code $(solution.retcode)")
    return solution.states[end]
end

"""
    gate_unitary(model, gate; kwargs...)

Compatibility wrapper for [`get_unitary`](@ref). Hamiltonians use frequency units;
the `2π` conversion is applied by `get_unitary` at the solver call.
"""
gate_unitary(args...; kwargs...) = get_unitary(args...; kwargs...)
