_display_unicode(io::IO) = get(io, :unicode, true)::Bool
_display_symbol(io::IO, unicode::AbstractString, ascii::AbstractString) =
    _display_unicode(io) ? unicode : ascii

_display_plural(count::Integer, singular::AbstractString, plural = string(singular, "s")) =
    count == 1 ? singular : plural

_display_sorted(values) = sort!(collect(values); by = string)

function _display_path(path)
    parts = path isa ParamPath ? path.parts : Tuple(path)
    return join(string.(parts), "/")
end

function _display_operator_path(path)
    parts = Tuple(path)
    if length(parts) >= 4 && parts[1] == :components && parts[3] == :operators
        return join(string.((parts[2], parts[4:end]...)), ".")
    elseif length(parts) >= 2 && parts[1] == :operators
        return join(string.(parts[2:end]), ".")
    end
    return join(string.(parts), ".")
end

function _display_items(io::IO, values)
    items = _display_sorted(values)
    isempty(items) && return (items, 0)
    get(io, :limit, false) || return (items, 0)
    rows, columns = displaysize(io)
    available = max(1, min(rows ÷ 4, columns ÷ 18))
    limit = min(length(items), available)
    return (items[1:limit], length(items) - limit)
end

function _display_metadata(io::IO, metadata)
    isempty(metadata) && return "none"
    keys_shown, omitted = _display_items(io, keys(metadata))
    text = join(string.(keys_shown), ", ")
    omitted > 0 && (text *= ", " * _display_symbol(io, "…", "..."))
    return string(length(metadata), " ", _display_plural(length(metadata), "key"), ": ", text)
end

_display_shape(value) = join(size(value), _display_shape_separator())
_display_shape_separator() = "×"
_display_shape(io::IO, value) = join(size(value), _display_symbol(io, "×", "x"))

function _display_value(io::IO, value)
    show(IOContext(io, :compact => true, :limit => true), value)
end

function _display_named_lines(io::IO, label::AbstractString, values, formatter)
    items = collect(values)
    if isempty(items)
        println(io, label, " none")
        return
    end
    println(io, label)
    leading = match(r"^\s*", label).match
    item_prefix = string(leading, "  ")
    shown, omitted = _display_items(io, items)
    for value in shown
        print(io, item_prefix)
        formatter(io, value)
        println(io)
    end
    omitted > 0 && println(io, item_prefix, _display_symbol(io, "…", "..."), " ", omitted, " more")
end

_display_named_lines(formatter, io::IO, label::AbstractString, values) =
    _display_named_lines(io, label, values, formatter)

function _display_domain_text(domain::Domain)
    data = domain.data
    if data === nothing
        return domain isa Domain{Any} ? "any" : string(nameof(typeof(domain).parameters[1]))
    elseif data === :real
        return "real"
    elseif data === :positive
        return "positive real"
    elseif data === :nonnegative
        return "nonnegative real"
    elseif data isa Tuple && length(data) == 2
        return string(data[1], "..", data[2])
    elseif data isa AbstractSet
        return string("{", join(_display_sorted(data), ", "), "}")
    end
    return string(data)
end
