function Base.show(io::IO, dressing::DressingSpec)
    print(io, "DressingSpec(", dressing.steps, " steps, ", dressing.schedule,
        ", exponent=", dressing.exponent, ")")
end

function Base.show(io::IO, ::MIME"text/plain", dressing::DressingSpec)
    println(io, "DressingSpec")
    println(io, "  Steps: ", dressing.steps)
    println(io, "  Schedule: ", dressing.schedule)
    println(io, "  Exponent: ", dressing.exponent)
    println(io, "  Minimum overlap: ", dressing.minimum_overlap)
    print(io, "  Diagnostics: ", dressing.diagnostics)
end

function Base.show(io::IO, spec::ModelSpec)
    print(io, "ModelSpec(:", spec.name, "; ", length(spec.subsystems), " ",
        _display_plural(length(spec.subsystems), "component"), ", ",
        length(spec.interactions), " ", _display_plural(length(spec.interactions), "interaction"),
        ", dimension=", prod(spec.dimension.size), ")")
end

function Base.show(io::IO, ::MIME"text/plain", spec::ModelSpec)
    println(io, "ModelSpec :", spec.name)
    println(io, "  Components: ", join(keys(spec.subsystems), ", "))
    println(io, "  Dimensions: ", spec.dimension.size)
    println(io, "  States: ", spec.states_to_keep === nothing ? "all" : join(spec.states_to_keep, ", "))
    println(io, "  Parameters:")
    for path in sort!(collect(keys(spec.parameters)); by = path -> path.parts)
        print(io, "    ", join(path.parts, "/"), " = ")
        _display_value(io, spec.defaults[path.parts])
        println(io)
    end
    println(io, "  Interactions:")
    for interaction in spec.interactions
        println(io, "    ", _expression_text(io, interaction))
    end
    print(io, "  Dressing: ")
    show(io, spec.dressingspec)
end

function Base.show(io::IO, model::QuantumDeviceModel)
    print(io, "QuantumDeviceModel(:", model.spec.name, "; dimension=",
        size(model.hamiltonian, 1), ", ", length(model.operators), " operators)")
end

function Base.show(io::IO, ::MIME"text/plain", model::QuantumDeviceModel)
    println(io, "QuantumDeviceModel :", model.spec.name)
    free = parameters(model)
    parameterized = model.hamiltonian isa QuantumObjectEvolution
    println(io, "  Hamiltonian: ", join(size(model.hamiltonian), "×"),
        parameterized ? " (QobjEvo)" : "")
    println(io, "  Free parameters: ",
        isempty(free) ? "none" :
        join((join(path.parts, "/") for path in sort!(collect(keys(free));
            by = path -> path.parts)), ", "))
    println(io, "  Operators: ", length(model.operators))
    println(io, "  States: ", length(model.states))
    println(io, "  Energies: ", length(model.energies))
    print(io, "  Minimum dressing overlap: ", minimum(model.dressing_res.overlaps))
end
