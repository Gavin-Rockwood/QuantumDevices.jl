align_tracking_phase(ket, amplitude) =
    iszero(amplitude) ? ket : ket * (conj(amplitude) / abs(amplitude))

function near_degenerate_flags(energies, atol, rtol)
    flags = falses(length(energies))
    order = sortperm(energies)
    for j = 2:length(order)
        a, b = order[j - 1], order[j]
        if abs(energies[a] - energies[b]) <= atol + rtol * max(abs(energies[a]), abs(energies[b]))
            flags[a] = true
            flags[b] = true
        end
    end
    return flags
end

saved_array(array, saved_count) = saved_count == 1 ? vec(array) : array
"""
    track_states(state_history; seed=nothing, seed_transport=nothing, keys=nothing,
                 other=(;), save_steps=:all,
                 transports=nothing, method=:quick, degeneracy_atol=1e-9,
                 degeneracy_rtol=1e-10, overlap_threshold=0.5, tracking_threads=1)

Follow normalized kets with one-to-one assignment and phase alignment. The default
`method=:quick` takes each reference's best still-available candidate in reference
order, breaking ties by candidate order. `method=:global` maximizes total squared
overlap across all references. Inputs share coordinates unless explicit
QuantumObject transports map between consecutive bases; projected reference norm loss is retained. Auxiliary
`other` histories follow the same assignments; optional finite real `energies`
provide degeneracy diagnostics. Keys default to `(0,), (1,), ...`.

An optional normalized `seed` supplies a separate step-zero reference collection.
Only its `m` trajectories are followed through each candidate sample (which may
have more states). Match seed to the first sample, then each matched, phase-aligned
collection to the next. Keys match the seed count. `seed_transport` maps seed
coordinates into the first sample; otherwise their bases must agree. It requires
`seed` and is separate from the consecutive-sample `transports`.
The seed is outside the saved history: `save_steps=[1]` denotes the first supplied
sample, `total_steps` excludes the seed, and auxiliary histories need no seed entry.
First-sample overlap accuracy includes seed matching, including projected norm loss.

Process every sample, retaining only `save_steps`. Multiple saved samples have
key rows and sample columns; one saved sample produces vectors. Accuracy is saved
at those samples; final diagnostics always summarize the entire processed path.
These are continuity diagnostics, not probabilities of correct labeling.
`tracking_threads` bounds default-pool overlap workers. Parallel execution caches
bounded windows of dense adjacent overlaps; matching and phase alignment remain
ordered. One worker retains the direct, allocation-light quick matcher without a
seed. Seeded transitions remain sequential and compute rectangular m-by-n overlaps;
workers share candidate-column batches within each transition, with task-owned
packing scratch and deterministic assignment/phase alignment.
"""
function track_states(
    state_history;
    seed = nothing,
    seed_transport = nothing,
    keys = nothing,
    other = (;),
    save_steps = :all,
    transports = nothing,
    method = :quick,
    degeneracy_atol = 1e-9,
    degeneracy_rtol = 1e-10,
    overlap_threshold = 0.5,
    tracking_threads = 1,
)
    validate_tracking_threads(tracking_threads)
    tracking_method(method)
    tracking_tolerances(degeneracy_atol, degeneracy_rtol, overlap_threshold)
    seeded = !isnothing(seed)
    !seeded && !isnothing(seed_transport) &&
        throw(ArgumentError("seed_transport requires a seed"))
    references = seeded ? collect(seed) : nothing
    seeded && validate_kets(references)
    history, auxiliary = prepare_tracking_inputs(state_history, other, transports;
        reference_count=seeded ? length(references) : nothing)
    n = seeded ? length(references) : length(first(history))
    steps = length(history)
    if seeded
        if isnothing(seed_transport)
            first(references).dimensions == first(first(history)).dimensions ||
                throw(DimensionMismatch("different seed and first-sample bases require seed_transport"))
        else
            validate_tracking_transport(seed_transport, references, first(history))
        end
    end
    match = method === :quick ? quick_overlap_assignment : global_overlap_assignment
    names = isnothing(keys) ? Tuple((i,) for i = 0:(n - 1)) : Tuple(keys)
    length(names) == n && allunique(names) || throw(
        ArgumentError(
            "initial keys must be unique and match the tracked reference count",
        ),
    )
    saved = saved_sample_indices(save_steps, steps)
    saved_count = length(saved)
    energyhistory = haskey(auxiliary, :energies) ? auxiliary.energies : nothing
    state_type = promote_type(map(eltype, history)...)
    if seeded
        phase_type = tracking_overlap_eltype(references, first(history), seed_transport)
        if !isnothing(transports)
            for transport in transports
                phase_type = promote_type(phase_type, eltype(transport.data))
            end
        end
        # A complex seed or transport may phase-align otherwise real candidates
        # into complex kets; retain that storage type in saved trajectory arrays.
        state_type = promote_type(state_type, eltype(references),
            typeof(first(first(history)) * one(phase_type)))
    end
    saved_states = Matrix{state_type}(undef, n, saved_count)
    saved_other = map(
        values -> Matrix{promote_type(map(eltype, values)...)}(undef, n, saved_count),
        auxiliary,
    )
    saved_indices = Matrix{Int}(undef, n, saved_count)
    accuracy = (
        overlaps = ones(n, saved_count),
        cumulative_log_overlaps = zeros(n, saved_count),
        minimum_overlaps = ones(n, saved_count),
        low_overlap = falses(n, saved_count),
        degenerate = falses(n, saved_count),
        ambiguous = falses(n, saved_count),
    )
    previous = seeded ? references : first(history)
    minimums = ones(n)
    logs = zeros(n)
    low = falses(n)
    degenerate = falses(n)
    assignment = collect(1:n)
    stepweights = ones(n)
    nextsaved = 1
    workers = seeded ? min(tracking_threads, maximum(length, history), Threads.nthreads(:default)) :
        min(tracking_threads, max(steps - 1, 0), Threads.nthreads(:default))
    seed_overlaps = nothing
    seed_scratch = seeded ? Any[nothing for _ in 1:workers] : nothing
    cached = !seeded && workers > 1 ? Any[nothing for _ in 1:(2 * workers)] : nothing
    scratch = !seeded && workers > 1 ? Any[nothing for _ in 1:workers] : nothing
    cache_start = 0
    cache_end = 0
    phases = !seeded && workers > 1 ? fill(one(complex(dot(first(previous), first(previous)))), n) : nothing
    # Match, reorder, update full-path metrics, then retain requested samples.
    for t = 1:steps
        if seeded || t > 1
            if seeded
                transport = t == 1 ? seed_transport : isnothing(transports) ? nothing : transports[t - 1]
                seed_overlaps = seeded_tracking_overlaps!(seed_overlaps, seed_scratch,
                    previous, history[t], transport; workers)
                assignment, amplitudes = match(seed_overlaps)
            elseif workers <= 1
                references = isnothing(transports) ? previous :
                             [transports[t - 1] * ket for ket in previous]
                assignment, amplitudes = match(references, history[t])
            else
                if t > cache_end
                    cache_start = t
                    cache_end = min(steps, t + 2 * workers - 1)
                    tracking_overlap_batch(history, transports, cache_start, cache_end;
                        matrices = cached, scratch)
                end
                raw = cached[t - cache_start + 1]
                overlaps = [conj(phases[i]) * raw[row, j]
                            for (i, row) in enumerate(assignment), j in axes(raw, 2)]
                assignment, amplitudes = match(overlaps)
                phases .= [iszero(a) ? one(a) : conj(a) / abs(a) for a in amplitudes]
            end
            previous = [
                align_tracking_phase(history[t][j], amplitude) for
                (j, amplitude) in zip(assignment, amplitudes)
            ]
            stepweights = abs2.(amplitudes)
            logs .+= log.(stepweights)
            minimums .= min.(minimums, stepweights)
            low .|= stepweights .< overlap_threshold
        end
        if !isnothing(energyhistory)
            flags = near_degenerate_flags(energyhistory[t], degeneracy_atol, degeneracy_rtol)
            degenerate .|= flags[assignment]
        end
        if nextsaved <= saved_count && t == saved[nextsaved]
            column = nextsaved
            saved_states[:, column] = previous
            saved_indices[:, column] = assignment .- 1
            for name in propertynames(saved_other)
                getproperty(saved_other, name)[:, column] = getproperty(auxiliary, name)[t][assignment]
            end
            accuracy.overlaps[:, column] = stepweights
            accuracy.cumulative_log_overlaps[:, column] = logs
            accuracy.minimum_overlaps[:, column] = minimums
            accuracy.low_overlap[:, column] = low
            accuracy.degenerate[:, column] = degenerate
            accuracy.ambiguous[:, column] = low .| degenerate
            nextsaved += 1
        end
    end
    diagnostics = [LabelDiagnostic(minimums[i], logs[i], low[i], degenerate[i]) for i = 1:n]
    return TrackingResult(
        names,
        saved,
        steps,
        saved_array(saved_states, saved_count),
        map(a -> saved_array(a, saved_count), saved_other),
        saved_array(saved_indices, saved_count),
        map(a -> saved_array(a, saved_count), accuracy),
        diagnostics,
    )
end

# Two raw result slots and one packing workspace per worker are reused across
# windows. Each stage assigns a workspace to one task, independent of thread IDs;
# the stage barrier makes reuse safe even when Julia migrates tasks.
function tracking_overlap_batch(history, transports, first_step, last_step;
    matrices = Any[nothing for _ in first_step:last_step],
    scratch = Any[nothing for _ in 1:cld(last_step - first_step + 1, 2)])
    for parity in (0, 1)
        @sync for t in first_step:last_step
            (t - 2) % 2 == parity || continue
            slot = t - first_step + 1
            workspace = cld(slot, 2)
            Threads.@spawn :default begin
                transport = isnothing(transports) ? nothing : transports[t - 1]
                amplitudes, packed = tracking_pair_overlaps!(matrices[slot], scratch[workspace],
                    history[t - 1], history[t], transport)
                all(a -> isfinite(abs2(a)), amplitudes) ||
                    throw(ArgumentError("nonfinite state overlaps"))
                matrices[slot] = amplitudes
                scratch[workspace] = packed
            end
        end
    end
    return matrices
end

# BLAS-compatible dense kets use reusable packed columns and in-place GEMM.
# Sparse and other array/element types retain generic matrix multiplication.
function tracking_pair_overlaps!(amplitudes, packed, references, candidates, transport)
    dense = all(ket -> ket.data isa StridedVector && eltype(ket.data) <: LinearAlgebra.BlasFloat,
        references) && all(ket -> ket.data isa StridedVector && eltype(ket.data) <: LinearAlgebra.BlasFloat,
        candidates)
    T = tracking_overlap_eltype(references, candidates, transport)
    if !dense || !(T <: LinearAlgebra.BlasFloat)
        left = reduce(hcat, [ket.data for ket in references])
        right = reduce(hcat, [ket.data for ket in candidates])
        !isnothing(transport) && (left = transport.data * left)
        return adjoint(left) * right, nothing
    end
    left_size = (length(first(references).data), length(references))
    right_size = (length(first(candidates).data), length(candidates))
    transported_size = (right_size[1], left_size[2])
    if isnothing(packed) || eltype(packed.left) != T ||
       size(packed.left) != left_size || size(packed.right) != right_size ||
       (!isnothing(transport) &&
        (isnothing(packed.transported) || size(packed.transported) != transported_size))
        packed = (left = Matrix{T}(undef, left_size), right = Matrix{T}(undef, right_size),
            transported = isnothing(transport) ? nothing : Matrix{T}(undef, transported_size))
    end
    for (j, ket) in enumerate(references)
        copyto!(view(packed.left, :, j), ket.data)
    end
    for (j, ket) in enumerate(candidates)
        copyto!(view(packed.right, :, j), ket.data)
    end
    left = packed.left
    if !isnothing(transport)
        mul!(packed.transported, transport.data, left)
        left = packed.transported
    end
    output_size = (length(references), length(candidates))
    if !(amplitudes isa AbstractMatrix{T}) || size(amplitudes) != output_size
        amplitudes = Matrix{T}(undef, output_size)
    end
    mul!(amplitudes, adjoint(left), packed.right)
    return amplitudes, packed
end

function tracking_overlap_eltype(references, candidates, transport)
    T = promote_type(mapreduce(ket -> eltype(ket.data), promote_type, references),
        mapreduce(ket -> eltype(ket.data), promote_type, candidates))
    return isnothing(transport) ? T : promote_type(T, eltype(transport.data))
end

# Seeded references depend on the previous assignment, so only candidate columns
# within the current transition run in parallel. Output is m×n, never n×n. Each
# task owns its packing workspace and disjoint output columns across migrations.
function seeded_tracking_overlaps!(amplitudes, scratch, references, candidates, transport; workers)
    T = tracking_overlap_eltype(references, candidates, transport)
    output_size = (length(references), length(candidates))
    if !(amplitudes isa Matrix{T}) || size(amplitudes) != output_size
        amplitudes = Matrix{T}(undef, output_size)
    end
    count = min(workers, length(candidates))
    function columns!(worker)
        firstcolumn = fld((worker - 1) * length(candidates), count) + 1
        lastcolumn = fld(worker * length(candidates), count)
        output = view(amplitudes, :, firstcolumn:lastcolumn)
        result, packed = tracking_pair_overlaps!(output, scratch[worker], references,
            view(candidates, firstcolumn:lastcolumn), transport)
        result === output || copyto!(output, result)
        scratch[worker] = packed
        return nothing
    end
    if count == 1
        columns!(1)
    else
        @sync for worker in 1:count
            Threads.@spawn :default columns!(worker)
        end
    end
    all(a -> isfinite(abs2(a)), amplitudes) ||
        throw(ArgumentError("nonfinite state overlaps"))
    return amplitudes
end
