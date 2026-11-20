# Stored configuration contains only user-facing options; sampled trajectories
# are transient build data. Keep all dressing defaults in this function.
function normalize_dressing_options(options)
    options isa NamedTuple || throw(ArgumentError("dressing_options must be a NamedTuple"))
    defaults = (method=:quick, nsteps=20, thread_count=1, trajectory=cubic_dressing_trajectory,
                degeneracy_atol=1e-9, degeneracy_rtol=1e-10, overlap_threshold=0.5)
    unknown = setdiff(keys(options), keys(defaults))
    isempty(unknown) || throw(ArgumentError("unknown dressing_options keys: $(join(unknown, ", "))"))
    config = merge(defaults, options)
    dressing_thread_count(config.thread_count)
    tracking_method(config.method)
    config.nsteps isa Integer && !(config.nsteps isa Bool) && config.nsteps >= 2 ||
        throw(ArgumentError("nsteps must be an integer >= 2"))
    tracking_tolerances(config.degeneracy_atol, config.degeneracy_rtol, config.overlap_threshold)
    applicable(config.trajectory, 0.0) || throw(ArgumentError("trajectory must be callable with a real sample"))
    return config
end

cubic_dressing_trajectory(x) = x^3

function dressing_samples(options)
    values = map(options.trajectory, range(0, 1; length=options.nsteps))
    all(x -> x isa Real && isfinite(x) && 0 <= x <= 1, values) &&
    first(values) == 0 && last(values) == 1 && issorted(values) || throw(
        ArgumentError("trajectory samples must be finite, monotone in [0,1], with exact endpoints 0 and 1"))
    return Float64.(values)
end

function dressing_thread_count(thread_count)
    thread_count isa Int && !(thread_count isa Bool) && thread_count > 0 ||
        throw(ArgumentError("thread_count must be a positive Int"))
    return thread_count
end

# Builder-only execution setting; never merge it into captured recipe options.
function validate_dressed_threads(dressed_threads)
    isnothing(dressed_threads) && return nothing
    dressed_threads isa Int && !(dressed_threads isa Bool) && dressed_threads > 0 ||
        throw(ArgumentError("dressed_threads must be a positive Int or nothing"))
    return dressed_threads
end

dressing_sample_workers(options, dressed_threads) =
    isnothing(dressed_threads) ? options.thread_count : dressed_threads

tracking_sample_workers(options, dressed_threads, tracking_threads) =
    ismissing(tracking_threads) ? dressing_sample_workers(options, dressed_threads) : tracking_threads
