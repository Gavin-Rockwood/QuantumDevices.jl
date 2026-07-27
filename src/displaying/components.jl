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
    return operator isa Symbol ? string(operator) : _display_operator_path(operator)
end

function _display_spec_line(io::IO, spec)
    component_parameters = parameters(spec)
    component_operators = operators(spec)
    print(
        io,
        nameof(typeof(spec)),
        "(",
        length(component_parameters),
        " ",
        _display_plural(length(component_parameters), "parameter"),
        ", ",
        length(component_operators),
        " ",
        _display_plural(length(component_operators), "operator"),
        ")",
    )
end

function _display_spec_plain(io::IO, spec)
    println(io, nameof(typeof(spec)))
    println(io, "Dimension: ", spec.dimension.size)
    component_parameters = parameters(spec)
    _display_named_lines(io, "Parameters:", keys(component_parameters)) do output, key
        _display_parameter_summary(output, component_parameters[key])
    end
    _display_named_lines(io, "Operators:", operators(spec)) do output, operator
        print(output, _display_operator_name(operator))
    end
    print(io, "Hamiltonian: ", _expression_text(io, hamiltonian(spec)))
end

for SpecType in (:QubitSpec, :TransmonSpec, :FluxTunableTransmonSpec, :ResonatorSpec)
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
    println(io, "  Dimension: ", spec.dimension.size)
    print(io, "  Spectrum: ")
    values = spec.spectrum
    shown, omitted = _display_items(io, values)
    print(io, join(shown, ", "))
    omitted > 0 && print(io, ", ", _display_symbol(io, "…", "..."), " ", omitted, " more")
    spec.source === nothing || print(io, "\n  Source: ", nameof(typeof(spec.source)))
end

function Base.show(io::IO, component::Component)
    print(
        io,
        "Component(:",
        component.name,
        ", ",
        nameof(typeof(component.spec)),
        "; ",
        length(parameters(component)),
        " ",
        _display_plural(length(parameters(component)), "parameter"),
        ", ",
        length(operators(component)),
        " ",
        _display_plural(length(operators(component)), "operator"),
        ")",
    )
end

function Base.show(io::IO, ::MIME"text/plain", component::Component)
    println(io, "Component :", component.name)
    println(io, "  Spec: ", nameof(typeof(component.spec)))
    println(io, "  Dimension: ", component.spec.dimension.size)
    component_parameters = parameters(component)
    _display_named_lines(io, "  Parameters:", keys(component_parameters)) do output, key
        _display_parameter_summary(output, component_parameters[key])
    end
    _display_named_lines(io, "  Operators:", operators(component)) do output, operator
        print(output, _display_operator_name(operator))
    end
    print(io, "  Hamiltonian: ", _expression_text(io, hamiltonian(component)))
end
