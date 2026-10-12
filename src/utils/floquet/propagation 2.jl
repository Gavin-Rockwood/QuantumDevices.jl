"""
    propagate_floquet_modes(modes_t0, H, t, T; t0=0, e_quasi=nothing,
                            propagator_kwargs=Dict{Symbol,Any}(), kwargs...)
    propagate_floquet_modes(basis::FloquetBasis, t)

Propagate initial modes over `mod(t-t0, T)`. Supplying quasienergies in frequency
units removes the dynamical phase and returns periodic Floquet modes. Without
`e_quasi`, retain the upstream behavior of returning evolved kets within a period.
Use the basis overload to obtain phase-corrected modes automatically.
"""
function propagate_floquet_modes(modes_t0, H::AbstractQuantumObject, t::Real,
    T::Real; t0::Real=0, e_quasi=nothing,
    propagator_kwargs=Dict{Symbol,Any}(), kwargs...)
    _floquet_period(T)
    isfinite(t) && isfinite(t0) || throw(ArgumentError("Floquet times must be finite"))
    e_quasi === nothing || length(e_quasi) == length(modes_t0) ||
        throw(DimensionMismatch("Quasienergy and mode counts differ"))
    elapsed = mod(t - t0, T)
    iszero(elapsed) && return modes_t0
    options = merge((; propagator_kwargs...), (; kwargs...))
    U = _floquet_propagator(H, elapsed, t0; options...)
    return [e_quasi === nothing ? U * state :
        exp(2pi * im * e_quasi[i] * elapsed) * (U * state)
        for (i, state) in enumerate(modes_t0)]
end

propagate_floquet_modes(basis::FloquetBasis, t::Real) = basis.modes(t)
