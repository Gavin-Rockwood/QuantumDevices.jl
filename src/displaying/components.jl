function _display_parameter_summary(io::IO, parameter::DeviceParameter)
    print(io, _display_path(parameter.path), " = ")
    _display_value(io, parameter.default)
    flags = String[]
    parameter.fixed && push!(flags, "fixed")
    parameter.required && push!(flags, "required")
    !isempty(flags) && print(io, " [", join(flags, ", "), "]")
    print(io, " ∈ ", _display_domain_text(parameter.domain))
end

function _display_operator_name(operator)
    return operator isa Symbol ? string(operator) : _display_path(operator)
end

function _display_spec_line(io::IO, spec)
    parameters = _spec_parameters(spec)
    operators = try
        available_operators(spec)
    catch
        Any[]
    end
    print(
        io,
        nameof(typeof(spec)),
        "(",
        length(parameters),
        " ",
        _display_plural(length(parameters), "parameter"),
        ", ",
        length(operators),
        " ",
        _display_plural(length(operators), "operator"),
        ")",
    )
end

function _display_spec_plain(io::IO, spec)
    println(io, nameof(typeof(spec)))
    parameters = _spec_parameters(spec)
    _display_named_lines(io, "Parameters:", keys(parameters)) do output, key
        _display_parameter_summary(output, parameters[key])
    end
    operators = try
        available_operators(spec)
    catch
        Any[]
    end
    _display_named_lines(io, "Operators:", operators) do output, operator
        print(output, _display_operator_name(operator))
    end
    hamiltonian = default_hamiltonian(spec)
    hamiltonian !== nothing && print(io, "Hamiltonian: ", _expression_text(io, hamiltonian))
end

for SpecType in (:QubitSpec, :TransmonSpec, :FluxTunableTransmonSpec, :ResonatorSpec, :FrozenModelSpec)
    @eval begin
        Base.show(io::IO, spec::$SpecType) = _display_spec_line(io, spec)
        Base.show(io::IO, ::MIME"text/plain", spec::$SpecType) = _display_spec_plain(io, spec)
    end
end

function Base.show(io::IO, spec::GenericSpec)
    print(io, "GenericSpec(", length(spec.spectrum), " levels)")
end

function Base.show(io::IO, ::MIME"text/plain", spec::GenericSpec)
    println(io, "GenericSpec")
    println(io, "  Dimension: ", spec.dimension.default)
    print(io, "  Spectrum: ")
    values = [parameter.default for parameter in spec.spectrum]
    shown, omitted = _display_items(io, values)
    print(io, join(shown, ", "))
    omitted > 0 && print(io, ", ", _display_symbol(io, "…", "..."), " ", omitted, " more")
end

function Base.show(io::IO, component::Component)
    print(
        io,
        "Component(:",
        component.id,
        ", ",
        nameof(typeof(component.spec)),
        "; ",
        length(component.parameters),
        " ",
        _display_plural(length(component.parameters), "parameter"),
        ", ",
        length(component.operators),
        " ",
        _display_plural(length(component.operators), "operator"),
        ")",
    )
end

function Base.show(io::IO, ::MIME"text/plain", component::Component)
    println(io, "Component :", component.id)
    println(io, "  Name: ", component.name)
    println(io, "  Spec: ", nameof(typeof(component.spec)))
    _display_named_lines(io, "  Parameters:", keys(component.parameters)) do output, key
        _display_parameter_summary(output, component.parameters[key])
    end
    _display_named_lines(io, "  Operators:", component.operators) do output, operator
        print(output, _display_operator_name(operator))
    end
    print(io, "  Hamiltonian: ")
    if component.hamiltonian === nothing
        print(io, "none")
    else
        print(io, _expression_text(io, component.hamiltonian))
    end
    print(io, "\n  Metadata: ", _display_metadata(io, component.metadata))
end
