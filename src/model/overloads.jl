function Base.getproperty(model::NumericalModel, name::Symbol)
    if name == :control_maps
        return getfield(model,:recipe).definition.control_maps
    elseif name == :parameters
        return recipe_parameter_view(getfield(model,:recipe);source=:model)
    elseif name == :variables
        return recipe_variables(getfield(model,:recipe))
    elseif name == :resolved_parameters
        return getfield(getfield(model,:recipe),:parameters)
    elseif name == :states
        return getfield(model, :eigensystem).states
    elseif name in (:energies, :spectrum)
        return getfield(model, :eigensystem).energies
    elseif name == :labeling
        tracking = getfield(model, :eigensystem).tracking
        return tracking.tracked ? tracking : nothing
    end
    return getfield(model, name)
end

Base.propertynames(model::NumericalModel, private::Bool = false) =
    (fieldnames(typeof(model))..., :parameters, :variables, :resolved_parameters, :states, :energies, :spectrum, :labeling, :control_maps)
"""Return energy-ordered eigenstates in the model construction basis."""
states(model::NumericalModel) =
    collect(model.states)

function Base.show(io::IO, model::NumericalModel)
    return print(
        io,
        "NumericalModel(owners=",
        sort!(collect(model.owners)),
        ", dimension=",
        length(model.spectrum),
        ", operators=",
        length(model.operators),
        ")",
    )
end

"""Zero-based energy index for a complete occupation tuple in model component order."""
stateindex(model::NumericalModel, key::Tuple) = tracking_row(model.eigensystem.tracking, key) - 1
stateindex(model::NumericalModel, first::Int, rest::Int...) = stateindex(model, (first, rest...))
"""Return the labeled eigenstate in model construction coordinates, as `model.states(key)`."""
state(model::NumericalModel, key::Tuple) = model.states(key)
state(model::NumericalModel, first::Int, rest::Int...) = state(model, (first, rest...))
"""Return the energy (E/h) associated with a complete occupation tuple."""
stateenergy(model::NumericalModel, key::Tuple) = model.energies(key)
stateenergy(model::NumericalModel, first::Int, rest::Int...) = stateenergy(model, (first, rest...))

function Base.show(io::IO, ::MIME"text/plain", model::NumericalModel)
    previewcount = min(length(model.spectrum), 6)
    preview =
        join(model.spectrum[1:previewcount], ", ") *
        (previewcount < length(model.spectrum) ? " … ($(length(model.spectrum)) total)" : "")
    labeling =
        isnothing(model.labeling) ? "unavailable" :
        "$(length(model.eigensystem.keys)) labels, $(count(d->d.ambiguous,model.labeling.diagnostics)) ambiguous, $(count(d->d.low_overlap,model.labeling.diagnostics)) low-overlap"
    return summary_lines(
        io,
        model,
        [
            "NumericalModel · dimension $(length(model.spectrum))",
            "  owners: " * list_summary(io, model.owners),
            "  spectrum (E/h): " * preview,
            "  operators ($(length(model.operators))): " * list_summary(io, keys(model.operators)),
            "  labeling: " * labeling,
            "  provenance: $(nameof(typeof(model.recipe.definition))) · resolved recipe with $(length(model.recipe.inputs)) immediate inputs",
        ],
    )
end

function model_identity(model::NumericalModel)
    n = length(model.spectrum)
    dimensions = Tuple(model.H.dims[1])
    return qeye(n; dims=dimensions)
end

"""Retrieve the model's physical symbolic Hamiltonian; numerical reference data remains in `.H`."""
get_hamiltonian(model::NumericalModel) = getfield(model, :hamiltonian)
