"""
    get_unitary(H, duration; kwargs...)
    get_unitary(model::DeviceModel, gate::DeviceGate; kwargs...)

Evolve the identity over `[0, duration]` and return the final retained-space
propagator. `H` is a static quantum operator or `QobjEvo` in frequency units
(cycles per unit time), e.g. GHz for time in ns. Evolution explicitly uses
`2π * H` in QuantumToolbox `sesolve`. The model/gate overload folds `2π`
into matrices using `numerical(model, gate; scalar=2pi)`.
Duration must be finite and positive.

`dense=true` (default) densifies the Hamiltonian and initial identity before
integration. For `QobjEvo`, matrix terms are converted while time-dependent
coefficients are preserved. `dense=false` preserves the supplied storage.
Solver kwargs pass through, with `progress_bar=false` by default. Unsuccessful
retcodes raise an error. The model/gate overload supplies `pulse_tstops(gate)`
by default so adaptive stepping resolves delayed controls. No rotating frame or rotating-wave approximation is applied.
"""
get_unitary(H::AbstractQuantumObject, duration::Real; kwargs...) =
    _get_unitary_angular(2pi * H, duration; kwargs...)

function _get_unitary_angular(H::AbstractQuantumObject, duration::Real; dense::Bool=true, kwargs...)
    isfinite(duration) && duration > 0 ||
        throw(ArgumentError("Unitary evolution requires a finite positive duration"))
    H = _numerical_storage(H, dense)
    initial = qeye_like(H isa QobjEvo ? H(0.0) : H)
    initial = _numerical_storage(initial, dense)
    options = merge((; progress_bar=false), (; kwargs...))
    solution = sesolve(H, initial, [zero(duration), duration]; options...)
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

# Adapt changes the operator-tree matrix storage without evaluating/freezing its
# time-dependent coefficients. Conversion happens once, before solver setup.
struct _DenseHamiltonianStorage end
Adapt.adapt_storage(::_DenseHamiltonianStorage, A::AbstractArray) = Array(A)

function _numerical_storage(H::QobjEvo, dense::Bool)
    dense || return H
    return QobjEvo(Adapt.adapt(_DenseHamiltonianStorage(), H.data), H.type, H.dimensions)
end
