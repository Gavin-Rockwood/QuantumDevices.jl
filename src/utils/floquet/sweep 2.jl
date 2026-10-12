"""
    floquet_sweep(H_func, sampling_points, T; sampling_times=[], use_logging=true,
                  states_to_track=nothing, propagator_kwargs=Dict{Symbol,Any}(),
                  t0=0, tracking_kwargs=NamedTuple(), kwargs...)

Compute Floquet modes and quasienergies for `H_func(point)`. `T` is one period
or a collection with one period per point. Sampling times default to `t0`.
Repeated `(point, period)` pairs reuse their basis. Returns a dictionary with
`"F_Modes"` and `"F_Energies"`; optional initial states (a dictionary or collection)
add `"Tracking"`, a `TrackingResult` with auxiliary `"Quasienergies"` histories.
Solver options pass to `get_floquet_basis`; `tracking_kwargs` pass to `track_states`.
"""
function floquet_sweep(H_func, sampling_points::AbstractArray,
    T::Union{Real,AbstractArray}; sampling_times=[], use_logging=true,
    states_to_track=nothing, propagator_kwargs=Dict{Symbol,Any}(),
    t0::Real=0, tracking_kwargs=NamedTuple(), kwargs...)
    points = vec(collect(sampling_points))
    periods = T isa Real ? fill(T, length(points)) : vec(collect(T))
    length(periods) == length(points) || throw(DimensionMismatch("Expected one period per sampling point"))
    foreach(_floquet_period, periods)
    isfinite(t0) || throw(ArgumentError("Reference time must be finite"))
    times = isempty(sampling_times) ? fill(t0, length(points)) : vec(collect(sampling_times))
    length(times) == length(points) || throw(DimensionMismatch("Expected one sampling time per point"))
    all(isfinite, times) || throw(ArgumentError("Sampling times must be finite"))
    use_logging && @info "Beginning Floquet sweep" steps=length(points)
    bases = Dict{Any,Any}()
    modes, energies = [], []
    for (point, period, time) in zip(points, periods, times)
        basis = get!(bases, (point, period)) do
            get_floquet_basis(H_func(point), period; t0, propagator_kwargs, kwargs...)
        end
        push!(modes, basis.modes(time))
        push!(energies, copy(basis.e_quasi))
    end
    result = Dict{String,Any}("F_Modes" => modes, "F_Energies" => energies)
    if states_to_track !== nothing && !isempty(states_to_track)
        result["Tracking"] = track_states(modes, states_to_track;
            other_sorts=Dict("Quasienergies" => energies), tracking_kwargs...)
    end
    use_logging && @info "Finished Floquet sweep"
    return result
end
