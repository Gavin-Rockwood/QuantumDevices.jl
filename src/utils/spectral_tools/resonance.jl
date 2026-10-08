# Adapted from SuperconductingCircuits.jl/src/Dynamics/Drives/ResonanceFinder.jl.
"""
    ResonanceResult

Floquet resonance estimate with `frequency`, `minimum_gap`, and approximate
two-state `drive_time=1/(2*minimum_gap)`. Frequencies and gaps use cycles per unit
time; drive time uses the reciprocal time unit. `state_keys`, sorted `frequencies`,
an N×2 matrix of tracked `quasienergies`, circular `gaps`, the native `tracking`
result, and an [`AvoidedCrossingFit`](@ref) in `fit` provide diagnostics.
"""
struct ResonanceResult{K,T,F}
    frequency::Float64
    minimum_gap::Float64
    drive_time::Float64
    state_keys::K
    frequencies::Vector{Float64}
    quasienergies::Matrix{Float64}
    gaps::Vector{Float64}
    tracking::T
    fit::F
end

function _resonance_references(references, state_keys)
    references isa Union{AbstractDict,AbstractVector} ||
        throw(ArgumentError("Reference states must be a dictionary or vector"))
    available = references isa AbstractDict ? collect(keys(references)) : collect(eachindex(references))
    selected = state_keys === nothing ? available : collect(state_keys)
    length(selected) == 2 && allunique(selected) ||
        throw(ArgumentError("Select exactly two distinct reference states using state_keys"))
    all(k -> k in available, selected) || throw(ArgumentError("Unknown reference state key"))
    states = [references[k] for k in selected]
    all(s -> s isa QuantumObject && isket(s), states) ||
        throw(ArgumentError("Reference states must be quantum kets"))
    # Track with integer keys so dictionary iteration never changes selected order.
    return selected, Dict(i => states[i] for i in 1:2)
end

_circular_quasienergy_gap(e1, e2, frequency) =
    abs(mod(e1 - e2 + frequency / 2, frequency) - frequency / 2)

"""
    find_resonance(H_func, freqs, reference_states; state_keys=nothing, t0=0,
        propagator_kwargs=(;), tracking_method=QuantumDevices.expansive_tracking,
        fit_kwargs=(;), use_logging=false)
    find_resonance(H0, drive_op, amplitude, freqs, reference_states; kwargs...)

Find one isolated, bracketed Floquet avoided crossing. `H_func(frequency)`
returns a Hamiltonian in frequency units with period `1/frequency`. The second
form uses `H0 + amplitude*sin(2π*frequency*t)*drive_op`.
At least five distinct, finite, positive frequencies are required; they are
sorted ascending before tracking. Reference states are a dictionary or vector;
use both entries by default, or select two labels/indices with `state_keys`.
`t0` is the absolute mode sampling/reference time.

Fit the quasienergy separation modulo frequency using [`fit_avoided_crossing`](@ref).
Returns a [`ResonanceResult`](@ref). The transfer-time estimate assumes an isolated
two-state resonance; inspect residuals and tracking confidence before using it.
Unbracketed/failed fits and gaps ≤ `sqrt(eps(Float64))*maximum(sampled_gaps)`
raise errors. Evolution failures propagate. Solver options pass through
`propagator_kwargs`; no plots are displayed automatically.
"""
function find_resonance(H_func, freqs, reference_states; state_keys=nothing,
    t0::Real=0, propagator_kwargs=(;), tracking_method=expansive_tracking,
    fit_kwargs=(;), use_logging=false)
    frequencies, _ = _spectral_parameters(freqs)
    all(>(0), frequencies) || throw(ArgumentError("Drive frequencies must be positive"))
    isfinite(t0) || throw(ArgumentError("Reference time must be finite"))
    selected, references = _resonance_references(reference_states, state_keys)
    sweep = floquet_sweep(H_func, frequencies, 1 ./ frequencies;
        t0, states_to_track=references, propagator_kwargs, use_logging,
        tracking_kwargs=(; method=tracking_method))
    tracking = sweep["Tracking"]
    # Auxiliary energies follow the same iteration order as the tracking input.
    tracked_keys = collect(keys(references))
    columns = [findfirst(==(i), tracked_keys) for i in 1:2]
    histories = tracking.others["Quasienergies"]
    quasienergies = [histories[j][columns[i]] for j in eachindex(frequencies), i in 1:2]
    # Restore user labels and diagnostic column order together.
    labeled_states = Dict(selected[i] => tracking.states[i] for i in 1:2)
    tracking = TrackingResult(labeled_states,
        Dict("Quasienergies" => [collect(row) for row in eachrow(quasienergies)]),
        [confidence[columns] for confidence in tracking.confidence])
    gaps = [_circular_quasienergy_gap(quasienergies[j, 1], quasienergies[j, 2], frequencies[j])
        for j in eachindex(frequencies)]
    fit = fit_avoided_crossing(frequencies, gaps; fit_kwargs)
    fit.minimum_gap > sqrt(eps(Float64)) * maximum(gaps) ||
        error("Fitted resonance gap is unresolved")
    drive_time = inv(2 * fit.minimum_gap)
    isfinite(drive_time) || error("Estimated drive time is not finite")
    result = ResonanceResult(fit.center, fit.minimum_gap, drive_time, selected,
        frequencies, quasienergies, gaps, tracking, fit)
    use_logging && @info "Floquet resonance" frequency=result.frequency minimum_gap=result.minimum_gap drive_time=result.drive_time
    return result
end

function find_resonance(H0::QuantumObject, drive_op::QuantumObject, amplitude::Real,
    freqs, reference_states; kwargs...)
    isoper(H0) && isoper(drive_op) || throw(ArgumentError("Hamiltonian and drive must be operators"))
    isfinite(amplitude) || throw(ArgumentError("Drive amplitude must be finite"))
    H_func = frequency -> H0 + QobjEvo(drive_op, (p, t) -> amplitude * sin(2pi * frequency * t))
    return find_resonance(H_func, freqs, reference_states; kwargs...)
end

function Base.show(io::IO, ::MIME"text/plain", result::ResonanceResult)
    print(io, "ResonanceResult\n  Frequency: ", result.frequency,
        "\n  Minimum gap: ", result.minimum_gap,
        "\n  Approximate two-state drive time: ", result.drive_time,
        "\n  RMS fit residual: ", result.fit.rms_residual)
end
