"""
    get_dressed_states(H0s, Hint; trajectory=x -> x^3, nsteps=20,
                       method=QuantumDevices.expansive_tracking)

Track bare product eigenstates while interpolating `sum(H0s) + g * Hint`.
`H0s` contains local Hamiltonians in tensor order and `Hint` is embedded in their
product space. At least two steps are required; the sampled trajectory must begin
at zero and end at one. Increase `nsteps` to test label convergence.

Returns a tracking result whose dictionaries use zero-based product-level tuples,
with final states, energies (`others`), and cumulative overlap confidence.
Level labels originate in local energy ordering, not computational array indices.
"""
function get_dressed_states(
    H0s,
    Hint;
    trajectory = x -> x^3,
    nsteps = 20,
    method = expansive_tracking,
)
    isempty(H0s) && throw(ArgumentError("At least one local Hamiltonian is required"))
    nsteps isa Integer && nsteps >= 2 || throw(ArgumentError("nsteps must be at least two"))
    preparation = _prepare_dressing(H0s)
    return _track_dressing(preparation, Hint; trajectory, nsteps, method)
end

# Read-only preparation shared by models with unchanged component Hamiltonians.
struct _DressingPreparation
    H0
    states
    energies
end

function _local_eigensystem(H)
    values, states = eigenstates(H)
    return (values, states)
end

function _prepare_dressing(H0s, eigensystems=[_local_eigensystem(H) for H in H0s])
    λ0s = [system[1] for system in eigensystems]
    ψ0s = [system[2] for system in eigensystems]

    # Labeled bare product states
    bare_states = Dict()
    bare_energies = Dict()

    level_ranges = map(ψ -> 0:length(ψ)-1, ψ0s)

    for levels in Iterators.product(level_ranges...)
        key = Tuple(levels)

        bare_states[key] = tensor(
            (ψ0s[i][levels[i] + 1] for i in eachindex(ψ0s))...
        )

        bare_energies[key] = sum(
            λ0s[i][levels[i] + 1]
            for i in eachindex(λ0s)
        )
    end

    dims = [size(H)[1] for H in H0s]
    H0 = sum([wrap(H0s[i], i, dims) for i in 1:length(dims)])

    return _DressingPreparation(H0, bare_states, bare_energies)
end

function _track_dressing(preparation::_DressingPreparation, Hint;
    trajectory=x -> x^3, nsteps=20, method=expansive_tracking)
    nsteps isa Integer && nsteps >= 2 || throw(ArgumentError("nsteps must be at least two"))
    gs = trajectory.(LinRange(0, 1, nsteps))
    first(gs) == 0 && last(gs) == 1 || throw(ArgumentError("Trajectory must start at zero and end at one"))

    H0, bare_states, bare_energies = preparation.H0, preparation.states, preparation.energies
    keys_ = collect(keys(bare_states))

    previous_states = [bare_states[key] for key in keys_]

    confidence = Dict(key => 1.0 for key in keys_)

    dressed_states = copy(bare_states)
    dressed_energies = copy(bare_energies)


    for g in gs[2:end]
        λ, ψ = eigenstates(H0 + g * Hint)

        result = track_states(
            [previous_states, ψ],
            Dict(key => previous_states[i] for (i, key) in enumerate(keys_));
            other_sorts = Dict("Energy" => [[dressed_energies[key] for key in keys_], λ]),
            method = method,
            return_steps = :end
            )

        current_states = [result.states[key][1] for key in keys_]

        current_energies = result.others["Energy"][1]

        for (i, key) in enumerate(keys_)
            dressed_states[key] = current_states[i]
            dressed_energies[key] = current_energies[i]

            confidence[key] *= result.confidence[1][i]
        end

        previous_states = current_states
    end

    return TrackingResult(dressed_states, dressed_energies, confidence)
end

function _tracking_step_count(steps, reference)
    reference === nothing || reference isa DeviceModel || throw(ArgumentError("Reference must be a DeviceModel or nothing"))
    steps === nothing && return reference === nothing ? 18 : 0
    steps isa Integer && !(steps isa Bool) && 0 <= steps < typemax(Int)-1 ||
        throw(ArgumentError("tracking_steps must be a nonnegative integer"))
    return Int(steps)
end

function _component_basis_overlap(old,new,d)
    old.dimension == new.dimension || throw(DimensionMismatch("Reference parent dimensions differ"))
    if old === new || old.operators === new.operators
        return qeye(d)
    end
    a,b = old.basis,new.basis
    a !== nothing && b !== nothing && a.id == b.id ||
        throw(ArgumentError("No common physical basis for component $(new.name)"))
    a.vectors === nothing && b.vectors === nothing && return qeye(d)
    A = a.vectors === nothing ? Matrix{ComplexF64}(I,old.dimension,d) : a.vectors[:,1:d]
    B = b.vectors === nothing ? Matrix{ComplexF64}(I,new.dimension,d) : b.vectors[:,1:d]
    return QuantumObject(B'*A)
end

function _reference_transport(reference,components,dims)
    length(reference.components)==length(components) &&
        all(a.name==b.name for (a,b) in zip(reference.components,components)) ||
        throw(ArgumentError("Reference component names and order differ"))
    old_dims = _model_dimensions(reference.components,reference.truncation_dimensions,reference.max_dimension)
    old_dims == dims || throw(DimensionMismatch("Reference retained dimensions differ"))
    return tensor((_component_basis_overlap(a,b,d) for (a,b,d) in zip(reference.components,components,dims))...)
end

# No transported reference normalization: projection loss belongs in overlap confidence.
function _match_model_reference(H,transport,reference)
    energies,states = eigenstates(H)
    labels = collect(keys(reference.states))
    length(labels)==length(states) || throw(DimensionMismatch("Reference state count differs"))
    old = [transport*reference.states[key] for key in labels]
    indices,overlaps = expansive_tracking(old,states)
    matched = Dict()
    spectrum = Dict()
    confidence = Dict()
    for (i,key) in enumerate(labels)
        state = states[indices[i]]
        phase = dot(old[i],state)
        matched[key] = iszero(abs(phase)) ? state : state*(conj(phase)/abs(phase))
        spectrum[key] = energies[indices[i]]
        confidence[key] = overlaps[i]
    end
    return TrackingResult(matched,spectrum,confidence)
end

function _interpolate_parameters(a,b,t)
    Set(keys(a))==Set(keys(b)) || throw(ArgumentError("Intermediate tracking requires identical parameter names"))
    return NamedTuple{keys(b)}(Tuple(map(keys(b)) do key
        x,y = a[key],b[key]
        x == y && return y
        x isa Number && y isa Number || throw(ArgumentError("Cannot interpolate parameter $key"))
        (1-t)*x+t*y
    end))
end

function _interpolate_component(a,b,t)
    a.type==b.type && a.hamiltonian==b.hamiltonian ||
        throw(ArgumentError("Intermediate tracking requires unchanged component structure"))
    a.parameters==b.parameters && return b
    parameters = _interpolate_parameters(a.parameters,b.parameters,t)
    if a.type == "TunableTransmon"
        parameters = NamedTuple{Tuple(k for k in keys(parameters) if k ∉ (:EJmax,:d))}(
            Tuple(parameters[k] for k in keys(parameters) if k ∉ (:EJmax,:d)))
    end
    return _rebuild_component(a,Dict(:parameters=>parameters))
end

function _track_model_reference(components,interaction,coupling,H,dims,reference,steps;
    max_dimension,threaded_components)
    transport = _reference_transport(reference,components,dims)
    previous = reference
    confidence = Dict(key=>1.0 for key in keys(reference.states))
    if steps>0
        interaction==reference.interactions || throw(ArgumentError("Intermediate tracking requires unchanged interactions"))
        for k in 1:steps
            t = k/(steps+1)
            mid_components = Vector{Component}(undef,length(components))
            make_component(i) = (mid_components[i]=_interpolate_component(reference.components[i],components[i],t))
            if threaded_components && length(components)>1 && Threads.nthreads()>1
                Threads.@threads for i in eachindex(components)
                    make_component(i)
                end
            else
                foreach(make_component,eachindex(components))
            end
            previous = make_model(mid_components,interaction,
                _interpolate_parameters(reference.coupling_parameters,coupling,t);
                max_dimension,truncation_dimensions=Dict(c=>d for (c,d) in zip(mid_components,dims)),
                reference_model=previous,tracking_steps=0,threaded_components)
            for key in keys(confidence)
                confidence[key] *= previous.confidence[key]
            end
        end
        transport = _reference_transport(previous,components,dims)
    end
    result = _match_model_reference(H,transport,previous)
    return TrackingResult(result.states,result.others,
        Dict(key=>confidence[key]*result.confidence[key] for key in keys(confidence)))
end
