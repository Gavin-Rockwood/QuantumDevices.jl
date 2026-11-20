Base.show(io::IO, result::TrackingResult) = print(
    io,
    "TrackingResult(",
    length(result.keys),
    " keys, ",
    length(result.saved_steps),
    " saved / ",
    result.total_steps,
    " processed samples, ",
    result.tracked ? "tracking available" : "untracked",
    ")",
)
Base.show(io::IO, data::TrackedData) = print(
    io,
    "TrackedData(objects=",
    size(data, 1),
    ", saved steps=",
    length(data.tracking.saved_steps),
    ")",
)
function tracked_key_preview(io, values, prefix; kind = "keys")
    items = [sprint(show, value; context = io) for value in values]
    if !get(io, :limit, false)
        return "[" * join(items, ", ") * "]"
    end
    available = max(8, displaysize(io)[2] - textwidth(prefix))
    selected = String[]
    for item in items
        candidate = "[" * join(vcat(selected, item), ", ") * "]"
        if textwidth(candidate) > available
            break
        end
        push!(selected, item)
    end
    if length(selected) < length(items)
        suffix = "… ($(length(items)) $kind)"
        while !isempty(selected) &&
              textwidth("[" * join(selected, ", ") * ", " * suffix * "]") > available
            pop!(selected)
        end
        return "[" * (isempty(selected) ? "" : join(selected, ", ") * ", ") * suffix * "]"
    end
    return "[" * join(selected, ", ") * "]"
end

function tracked_key_ranges(io, keys)
    isempty(keys) && return ["    unavailable"]
    integer_tuples = all(
        key ->
            key isa Tuple &&
            length(key) == length(first(keys)) &&
            all(value -> value isa Integer && !(value isa Bool), key),
        keys,
    )
    if !integer_tuples || isempty(first(keys))
        return ["    " * tracked_key_preview(io, keys, "    ")]
    end
    lines = String[]
    positions = length(first(keys))
    shown = get(io, :limit, false) ? min(positions, max(1, displaysize(io)[1] - 7)) : positions
    for position = 1:shown
        values = sort!(unique([key[position] for key in keys]))
        prefix = "    $position: "
        description = if length(values) == 1
            string(only(values))
        elseif all(i -> values[i] == first(values) + i - 1, eachindex(values))
            string(first(values), ":", last(values))
        else
            tracked_key_preview(io, values, prefix; kind = "values")
        end
        push!(lines, prefix * description)
    end
    shown < positions && push!(lines, "    … ($positions key positions)")
    return lines
end

function Base.show(io::IO, ::MIME"text/plain", data::TrackedData)
    get(io, :compact, false) && return show(io, data)
    lines = [
        "TrackedData",
        "  Number of Tracked Objects: $(size(data,1))",
        "  Number of Saved Steps:     $(length(data.tracking.saved_steps))",
        "",
        "  Key ranges:",
    ]
    append!(lines, tracked_key_ranges(io, data.tracking.keys))
    return print(io, join(lines, "\n"))
end

function Base.show(io::IO, ::MIME"text/plain", result::TrackingResult)
    quality =
        result.tracked ?
        "$(count(d->d.ambiguous,result.diagnostics)) ambiguous · $(count(d->d.low_overlap,result.diagnostics)) low-overlap" :
        "tracking unavailable"
    return summary_lines(
        io,
        result,
        [
            "TrackingResult · $(length(result.keys)) keys",
            "  saved samples: $(result.saved_steps) · processed samples: $(result.total_steps)",
            "  " * quality,
            "  saved arrays have key rows and sample columns; a single saved sample is a vector",
        ],
    )
end
