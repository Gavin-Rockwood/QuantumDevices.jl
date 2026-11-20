include("options.jl")
include("DressedStates.jl")
include("overloads.jl")

# Prepare product labels, generate independent ramp eigensystems, then track and
# organize them in order. Model construction consumes the resulting eigensystem.
function composite_eigensystem(recipe,models,selections,initial,interaction,reference;dressed_threads=nothing,tracking_threads=missing)
    if recipe_stage_operation(recipe,[model.recipe.definition for model in models]) === :projection
        return project_single_eigensystem(recipe,only(models),only(selections))
    end
    options = recipe.dressing_options
    if !recipe.dressing
        components = Tuple(component for model in models for component in model.eigensystem.components)
        if iszero(interaction.data)
            diagonal = real.(diag(initial.data))
            permutation = sortperm(diagonal; alg=Base.Sort.MergeSort)
            states = [QuantumObject(basis(length(diagonal),i-1).data; type=Ket(),
                dims=(Tuple(initial.dims[1]),(1,))) for i in permutation]
            return untracked_state_eigensystem(diagonal[permutation], states, components)
        end
        energies, rotation = diagonalize(reference)
        return untracked_eigensystem(energies, rotation, components)
    end
    preparation = prepare_product_dressing(models,selections,length.(selections))
    histories, energyhistory = dressing_sample_histories(preparation,initial,interaction,reference,options;dressed_threads)
    tracking = track_states(
        histories;
        keys = preparation.keys,
        other = (energies = energyhistory,),
        save_steps = :end,
        method = options.method,
        tracking_threads = tracking_sample_workers(options, dressed_threads, tracking_threads),
        degeneracy_atol = options.degeneracy_atol,
        degeneracy_rtol = options.degeneracy_rtol,
        overlap_threshold = options.overlap_threshold,
    )
    organized = energy_ordered_tracking(inherited_diagnostics(tracking, preparation.inherited))
    return DressedStates(organized, preparation.components)
end

function leaf_eigensystem(owner, energies, rotation)
    flags = near_degenerate_flags(energies, 1e-9, 1e-10)
    diagnostics = [LabelDiagnostic(1.0, 0.0, false, flag) for flag in flags]
    count = length(energies)
    accuracy = (
        overlaps = ones(count),
        cumulative_log_overlaps = zeros(count),
        minimum_overlaps = ones(count),
        low_overlap = falses(count),
        degenerate = flags,
        ambiguous = copy(flags),
    )
    tracking = TrackingResult(
        Tuple((i,) for i = 0:(count - 1)),
        (1,),
        1,
        eigensystem_kets(rotation, length(energies)),
        (energies = energies,),
        collect(0:(count - 1)),
        accuracy,
        diagnostics,
    )
    return DressedStates(tracking, (owner,))
end

function untracked_eigensystem(energies, rotation, components)
    return untracked_state_eigensystem(energies, eigensystem_kets(rotation, length(energies)), components)
end

function untracked_state_eigensystem(energies, states, components)
    tracking = TrackingResult(
        (),
        (1,),
        1,
        states,
        (energies = collect(energies),),
        collect(0:(length(energies) - 1)),
        nothing,
        nothing;
        tracked = false,
    )
    return DressedStates(tracking, Tuple(components))
end

function eigensystem_kets(rotation, count)
    count isa Integer && !(count isa Bool) && 0 <= count <= size(rotation, 2) ||
        throw(ArgumentError("count must be an integer between zero and the rotation column count"))
    dimensions = (Tuple(rotation.dims[1]), (1,))
    # Column indexing copies storage, so saved kets remain independent of the
    # rotation and each other, as with the former matrix-vector products.
    return [QuantumObject(rotation.data[:, index]; type = Ket(), dims = dimensions)
            for index in 1:count]
end

# Last tensor factor varies fastest, matching QuantumToolbox.tensor/basis.
function product_tracking_keys(models, selections)
    keys = [()]
    diagnostics = [LabelDiagnostic(1.0, 0.0, false, false)]
    components = Tuple(component for model in models for component in model.eigensystem.components)
    for (model, selected) in zip(models, selections)
        tracking = model.eigensystem.tracking
        tracking.tracked || throw(
            ArgumentError(
                "a labeled parent requires upstream labeling; enable dressing on the upstream recipe or set dressing=false on the parent",
            ),
        )
        nextkeys = typeof((first(keys)..., first(tracking.keys)...))[]
        nextdiagnostics = LabelDiagnostic[]
        for (key, diagnostic) in zip(keys, diagnostics), column in selected
            childkey = tracking.keys[column]
            childdiagnostic = tracking.diagnostics[column]
            push!(nextkeys, (key..., childkey...))
            push!(
                nextdiagnostics,
                LabelDiagnostic(
                    min(diagnostic.minimum_overlap, childdiagnostic.minimum_overlap),
                    diagnostic.cumulative_log_overlap + childdiagnostic.cumulative_log_overlap,
                    diagnostic.low_overlap || childdiagnostic.low_overlap,
                    diagnostic.degenerate || childdiagnostic.degenerate,
                ),
            )
        end
        keys = nextkeys
        diagnostics = nextdiagnostics
    end
    return keys, diagnostics, components
end

# Product kets preserve physical occupations even when the initial Hamiltonian
# is degenerate; do not replace them with a newly diagonalized initial basis.
function prepare_product_dressing(models, selections, dimensions)
    keys, inherited, components = product_tracking_keys(models, selections)
    initialenergies = Float64[0.0]
    partial = []
    for (index, dimension) in enumerate(dimensions)
        if index == 1
            partial = [basis(dimension, i) for i = 0:(dimension - 1)]
        else
            partial = [
                tensor(previous, basis(dimension, i)) for previous in partial for
                i = 0:(dimension - 1)
            ]
        end
        initialenergies =
            [a + b for a in initialenergies for b in models[index].spectrum[selections[index]]]
    end
    return (; keys, inherited, components, states = partial, energies = initialenergies)
end

# Each task owns a contiguous, balanced batch. Spawn explicitly in the default
# pool, including when the caller runs on an interactive thread.
function foreach_dressing_sample(f, sample_count, thread_count)
    dressing_thread_count(thread_count)
    workers = min(thread_count, sample_count, Threads.nthreads(:default))
    if workers <= 1
        for index in 1:sample_count
            f(index)
        end
    else
        @sync for batch in 1:workers
            first_index = fld((batch - 1) * sample_count, workers) + 1
            last_index = fld(batch * sample_count, workers)
            Threads.@spawn :default begin
                for index in first_index:last_index
                    f(index)
                end
            end
        end
    end
    return nothing
end

# Each diagonalization starts from the same initial Hamiltonian. Disjoint history
# slots retain ramp order for the sequential tracking stage.
function dressing_sample_histories(preparation,initial,interaction,reference,options;dressed_threads=nothing)
    count = length(preparation.keys)
    if iszero(interaction.data)
        permutation = sortperm(preparation.energies; alg=Base.Sort.MergeSort)
        return [preparation.states, preparation.states[permutation]],
               [preparation.energies, preparation.energies[permutation]]
    end
    samples = dressing_samples(options)[2:end]
    histories = Vector{typeof(preparation.states)}(undef, length(samples) + 1)
    energyhistory = Vector{Vector{Float64}}(undef, length(samples) + 1)
    histories[1] = preparation.states
    energyhistory[1] = preparation.energies
    foreach_dressing_sample(length(samples), dressing_sample_workers(options, dressed_threads)) do index
        sampleH = index == length(samples) ? reference : initial + samples[index]*interaction
        sample_energies, sample_rotation = diagonalize(sampleH)
        histories[index + 1] = eigensystem_kets(sample_rotation, count)
        energyhistory[index + 1] = sample_energies
    end
    return histories, energyhistory
end

function inherited_diagnostics(result, inherited)
    diagnostics = [
        LabelDiagnostic(
            min(upstream.minimum_overlap, localdiagnostic.minimum_overlap),
            upstream.cumulative_log_overlap + localdiagnostic.cumulative_log_overlap,
            upstream.low_overlap || localdiagnostic.low_overlap,
            upstream.degenerate || localdiagnostic.degenerate,
        ) for (upstream, localdiagnostic) in zip(inherited, getfield(result, :_diagnostics))
    ]
    accuracy = merge(
        getfield(result, :_accuracy),
        (
            minimum_overlaps = [d.minimum_overlap for d in diagnostics],
            cumulative_log_overlaps = [d.cumulative_log_overlap for d in diagnostics],
            low_overlap = [d.low_overlap for d in diagnostics],
            degenerate = [d.degenerate for d in diagnostics],
            ambiguous = [d.ambiguous for d in diagnostics],
        ),
    )
    return TrackingResult(
        result.keys,
        result.saved_steps,
        result.total_steps,
        getfield(result, :_states),
        getfield(result, :_other),
        getfield(result, :_indices),
        accuracy,
        diagnostics,
    )
end

function energy_ordered_tracking(result)
    permutation = sortperm(getfield(result, :_indices))
    return TrackingResult(
        Tuple(result.keys[i] for i in permutation),
        result.saved_steps,
        result.total_steps,
        getfield(result, :_states)[permutation],
        map(a -> a[permutation], getfield(result, :_other)),
        getfield(result, :_indices)[permutation],
        map(a -> a[permutation], getfield(result, :_accuracy)),
        getfield(result, :_diagnostics)[permutation],
    )
end

# Selection already established the retained energy order. Inherit those exact
# indices/diagnostics and expose coordinate kets without sorting or solving again.
function project_single_eigensystem(recipe,model,selected)
    energies = collect(model.spectrum[selected])
    # Selected upstream eigenvectors define the new coordinates. No redundant
    # rotation matrix, diagonalization, or tracking ramp is constructed here.
    basisstates = [basis(length(selected), i-1) for i in 1:length(selected)]
    upstream = model.eigensystem.tracking
    components = model.eigensystem.components
    if !recipe.dressing
        return untracked_state_eigensystem(energies, basisstates, components)
    end
    upstream.tracked || throw(ArgumentError("a labeled parent requires upstream labeling; enable dressing on the upstream recipe or set dressing=false on the parent"))
    accuracy = map(a -> a[selected], getfield(upstream, :_accuracy))
    tracking = TrackingResult(
        Tuple(upstream.keys[i] for i in selected), upstream.saved_steps, upstream.total_steps,
        basisstates, (energies=energies,),
        collect(0:(length(selected)-1)), accuracy, getfield(upstream, :_diagnostics)[selected])
    return DressedStates(tracking, components)
end
