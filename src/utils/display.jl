# Inspection summaries never traverse numerical caches or call operator factories.
function list_summary(io, values; ordered = false)
    items = string.(collect(values))
    !ordered && sort!(items)
    isempty(items) && return "none"
    budget = max(1, min(8, displaysize(io)[1]÷3))
    if get(io, :limit, false) && length(items) > budget
        return join(items[1:budget], ", ") * " … ($(length(items)) total)"
    end
    return join(items, ", ")
end

function summary_lines(io, object, lines)
    return get(io, :compact, false) ? show(io, object) : print(io, join(lines, "\n"))
end

ordered_parameters(parameters) =
    sort!(collect(pairs(parameters)); by = p -> string(first(p)))
parameter_values(parameters) =
    isempty(parameters) ? "none" :
    join((string(p, " = ", v) for (p, v) in ordered_parameters(parameters)), ", ")
