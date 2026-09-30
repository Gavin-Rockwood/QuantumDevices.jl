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
    gs = trajectory.(LinRange(0, 1, nsteps))
    first(gs) == 0 && last(gs) == 1 || throw(ArgumentError("Trajectory must start at zero and end at one"))

    λ0s = []
    ψ0s = []

    for H in H0s
        λ, ψ = eigenstates(H)
        push!(λ0s, λ)
        push!(ψ0s, ψ)
    end

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
