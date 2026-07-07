function Base.show(io::IO, path::ParamPath)
    print(io, _display_path(path))
end

Base.show(io::IO, ::MIME"text/plain", path::ParamPath) = show(io, path)

function Base.show(io::IO, domain::Domain)
    print(io, "Domain(", _display_domain_text(domain), ")")
end

function Base.show(io::IO, ::MIME"text/plain", domain::Domain)
    print(io, _display_domain_text(domain))
end

function Base.show(io::IO, parameter::DeviceParameter)
    print(io, "Parameter(", _display_path(parameter.path), "=", parameter.default)
    parameter.fixed && print(io, ", fixed")
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", parameter::DeviceParameter)
    println(io, "DeviceParameter ", _display_path(parameter.path))
    print(io, "  Default: ")
    _display_value(io, parameter.default)
    println(io)
    println(io, "  Domain: ", _display_domain_text(parameter.domain))
    flags = String[]
    parameter.fixed && push!(flags, "fixed")
    parameter.required && push!(flags, "required")
    println(io, "  Flags: ", isempty(flags) ? "none" : join(flags, ", "))
    !isempty(parameter.description) && println(io, "  Description: ", parameter.description)
    print(io, "  Metadata: ", _display_metadata(io, parameter.metadata))
end

function Base.show(io::IO, dimension::QDDimension)
    print(io, "QDDimension(", dimension.dim, ")")
end

Base.show(io::IO, ::MIME"text/plain", dimension::QDDimension) = show(io, dimension)
