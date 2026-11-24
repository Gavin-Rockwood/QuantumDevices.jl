Base.show(io::IO, d::QuantumDevice) = print(
    io,
    "QuantumDevice(",
    d.name,
    ", components=",
    length(d.components),
    ", interactions=",
    length(d.interactions),
    ", recipes=",
    length(d.recipes),
    ")",
)

function Base.show(io::IO, ::MIME"text/plain", device::QuantumDevice)
    lines = ["QuantumDevice $(device.name)"]
    components = sort!(collect(device.components); by = p -> string(first(p)))
    push!(lines, "  components ($(length(components))):")
    isempty(components) && push!(lines, "    none")
    for (name, c) in components
        label = isnothing(c.label) || c.label == string(name) ? "" : " · label=$(repr(c.label))"
        push!(lines, "    $name · $(c.type)$label")
    end
    push!(
        lines,
        "  registered physical interactions: " *
        list_summary(io, keys(device.interactions)),
    )
    push!(lines, "  registered recipes ($(length(device.recipes))):")
    isempty(device.recipes) && push!(lines, "    none")
    for (name, r) in sort!(collect(device.recipes); by = p -> string(first(p)))
        push!(lines, "    $name · components: " * coverage_description(r))
    end
    push!(lines, "  recipe-local interactions are stored on their recipes")
    return summary_lines(io, device, lines)
end

function Base.getproperty(device::QuantumDevice, name::Symbol)
    if name === :defaults
        return device_defaults(device)
    elseif name === :parameters
        return device_parameter_view(device)
    elseif name === :variables
        return device_variables(device)
    end
    return getfield(device,name)
end

Base.propertynames(::QuantumDevice, private::Bool=false) = (fieldnames(QuantumDevice)..., :variables, :parameters, :defaults)
