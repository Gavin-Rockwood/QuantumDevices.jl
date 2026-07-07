function Base.show(io::IO, interaction::InteractionSpec)
    print(io, "InteractionSpec(:", interaction.id, ", ", _expression_text(io, interaction.expr), ")")
end

function Base.show(io::IO, ::MIME"text/plain", interaction::InteractionSpec)
    println(io, "InteractionSpec :", interaction.id)
    println(io, "  Expression: ", _expression_text(io, interaction.expr))
    participants = _display_sorted(_participating_components(interaction.expr))
    println(io, "  Components: ", isempty(participants) ? "none" : join(participants, ", "))
    print(io, "  Metadata: ", _display_metadata(io, interaction.metadata))
end

function Base.show(io::IO, dressing::DressingSpec)
    print(
        io,
        "DressingSpec(",
        dressing.steps,
        " steps, ",
        dressing.schedule,
        ", exponent=",
        dressing.exponent,
        ")",
    )
end

function Base.show(io::IO, ::MIME"text/plain", dressing::DressingSpec)
    println(io, "DressingSpec")
    println(io, "  Steps: ", dressing.steps)
    println(io, "  Schedule: ", dressing.schedule)
    println(io, "  Exponent: ", dressing.exponent)
    println(io, "  Minimum overlap: ", dressing.minimum_overlap)
    print(io, "  Low-overlap policy: ", dressing.on_low_overlap)
end

function Base.show(io::IO, spec::ModelSpec)
    print(
        io,
        "ModelSpec(:",
        spec.id,
        "; ",
        length(spec.children),
        " ",
        _display_plural(length(spec.children), "child", "children"),
        ", ",
        length(spec.interactions),
        " ",
        _display_plural(length(spec.interactions), "interaction"),
        ", dimension=",
        prod(values(spec.dims)),
        ", selected=",
        length(spec.states_to_keep),
        ")",
    )
end

function Base.show(io::IO, ::MIME"text/plain", spec::ModelSpec)
    seen = Set{UInt}()
    _display_model_tree(io, spec, "", true, seen, true)
end

function _display_model_tree(io, spec, prefix, is_last, seen, root = false)
    marker = root ? "" : (is_last ? _display_symbol(io, "└─ ", "`- ") : _display_symbol(io, "├─ ", "|- "))
    identity = objectid(spec)
    if identity in seen
        print(io, prefix, marker, "ModelSpec :", spec.id, " ", _display_symbol(io, "↩", "<ref>"))
        return
    end
    push!(seen, identity)
    component_ids = union(Set(keys(spec.components)), Set(keys(spec.children)))
    component_count = length(component_ids)

    print(io, prefix, marker, "ModelSpec :", spec.id)
    print(
        io,
        " [dims=",
        join((string(key, "=", spec.dims[key]) for key in _display_sorted(keys(spec.dims))), ","),
        "; dimension=",
        prod(values(spec.dims)),
        "]",
    )
    print(
        io,
        " (",
        component_count,
        " ",
        _display_plural(component_count, "component"),
        ", ",
        length(spec.interactions),
        " ",
        _display_plural(length(spec.interactions), "interaction"),
        ")",
    )

    details_prefix = root ?
        string(prefix, "  ") :
        string(prefix, is_last ? "   " : _display_symbol(io, "│  ", "|  "))
    !isempty(component_ids) && print(
        io,
        "\n",
        details_prefix,
        "components: ",
        join(
            (
                string(key, "[", spec.dims[key], "]")
                for key in _display_sorted(component_ids)
            ),
            ", ",
        ),
    )
    !isempty(spec.initialization_dims) && print(
        io,
        "\n",
        details_prefix,
        "initialization: ",
        join(
            (
                string(
                    key isa Tuple ? join(key, ".") : key,
                    "=",
                    spec.initialization_dims[key],
                )
                for key in _display_sorted(keys(spec.initialization_dims))
            ),
            ", ",
        ),
    )
    !isempty(spec.states_to_keep) && print(
        io,
        "\n",
        details_prefix,
        "states: ",
        join(string.(spec.states_to_keep), ", "),
        "\n",
        details_prefix,
        "dressing: ",
        sprint(show, spec.dressing; context = IOContext(io, :compact => true)),
    )
    !isempty(spec.parameters) && print(
        io,
        "\n",
        details_prefix,
        "parameters: ",
        join(
            (
                sprint(show, spec.parameters[key].default; context = IOContext(io, :compact => true)) |>
                value -> string(key, "=", value)
                for key in _display_sorted(keys(spec.parameters))
            ),
            ", ",
        ),
    )
    !isempty(spec.interactions) && print(
        io,
        "\n",
        details_prefix,
        "interactions: ",
        join(
            (string(key, ": ", _expression_text(io, spec.interactions[key].expr)) for key in _display_sorted(keys(spec.interactions))),
            "; ",
        ),
    )

    children = _display_sorted(keys(spec.children))
    shown, omitted = _display_items(io, children)
    for (index, child_id) in enumerate(shown)
        print(io, "\n")
        child_is_last = index == length(shown) && omitted == 0
        _display_model_tree(io, spec.children[child_id], details_prefix, child_is_last, seen)
    end
    omitted > 0 && print(
        io,
        "\n",
        details_prefix,
        _display_symbol(io, "└─ … ", "`- ... "),
        omitted,
        " ",
        _display_plural(omitted, "descendant"),
    )
end

function Base.show(io::IO, model::QuantumDeviceModel)
    full_dimension = size(model.hamiltonian, 1)
    print(
        io,
        "QuantumDeviceModel(:",
        model.spec.id,
        "; basis=",
        model.basis,
        ", dimension=",
        full_dimension,
        ", ",
        length(model.operators),
        " ",
        _display_plural(length(model.operators), "operator"),
        ")",
    )
end

function Base.show(io::IO, ::MIME"text/plain", model::QuantumDeviceModel)
    println(io, "QuantumDeviceModel :", model.spec.id)
    println(io, "  Basis: ", model.basis)
    println(io, "  Hamiltonian: ", _display_shape(io, model.hamiltonian))
    println(io, "  Dressed states: ", length(model.states), " labeled states")
    println(io, "  Energies: ", length(model.energies), " labeled values")
    println(
        io,
        "  State order: ",
        isempty(model.state_order) ? "none" : join(model.state_order, ", "),
    )
    if isempty(model.spec.states_to_keep)
        println(io, "  State labels: ", length(model.state_labels), " product labels")
    else
        println(io, "  Selected labels: ", join(string.(model.state_labels), ", "))
    end
    println(
        io,
        "  Minimum tracked overlap: ",
        minimum(model.state_overlaps),
    )
    _display_named_lines(io, "  Operators:", keys(model.operators)) do output, path
        print(output, _display_path(path), " [", _display_shape(output, model.operators[path]), "]")
    end
    if isempty(model.resolved_parameters)
        println(io, "  Resolved parameters: none")
    else
        println(
            io,
            "  Resolved parameters: ",
            join(
                (
                    sprint(show, model.resolved_parameters[key]; context = IOContext(io, :compact => true)) |>
                    value -> string(key, "=", value)
                    for key in _display_sorted(keys(model.resolved_parameters))
                ),
                ", ",
            ),
        )
    end
    println(io, "  Provenance:")
    _display_model_tree(io, model.spec, "    ", true, Set{UInt}(), true)
end
