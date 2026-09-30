function Base.show(io::IO, result::TrackingResult)
    print(io, "TrackingResult(", length(result.states), " states, ",
        length(result.confidence), " saved steps)")
end

function Base.show(io::IO, ::MIME"text/plain", result::TrackingResult)
    show(io, result)
    limited = IOContext(io, :compact => true, :limit => true)
    print(io, "\n  State labels: ")
    show(limited, collect(keys(result.states)))
    print(io, "\n  Auxiliary quantities: ")
    show(limited, collect(keys(result.others)))
    print(io, "\n  Final confidence: ")
    if isempty(result.confidence)
        print(io, "unavailable")
    else
        show(limited, last(result.confidence))
    end
end
