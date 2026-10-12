function Base.show(io::IO, result::TrackingResult)
    if result.confidence isa AbstractDict
        print(io, "TrackingResult(", length(result.states), " dressed states)")
    else
        print(io, "TrackingResult(", length(result.states), " states, ",
            length(result.confidence), " saved steps)")
    end
end

function Base.show(io::IO, ::MIME"text/plain", result::TrackingResult)
    show(io, result)
    limited = IOContext(io, :compact => true, :limit => true)
    print(io, "\n  State labels: ")
    show(limited, collect(keys(result.states)))
    if result.confidence isa AbstractDict
        print(io, "\n  Energies (frequency units): ")
        show(limited, result.others)
    else
        print(io, "\n  Auxiliary quantities: ")
        show(limited, collect(keys(result.others)))
    end
    print(io, "\n  Final confidence: ")
    if isempty(result.confidence)
        print(io, "unavailable")
    else
        show(limited, result.confidence isa AbstractDict ? result.confidence : last(result.confidence))
    end
end
