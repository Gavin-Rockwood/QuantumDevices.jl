# Dressed-state forwarding
function Base.getproperty(dressed::DressedStates, name::Symbol)
    if name == :states
        return getfield(dressed, :tracking).states
    elseif name == :energies
        return getfield(dressed, :tracking).other.energies
    elseif name == :keys
        return getfield(dressed, :tracking).keys
    end
    return getfield(dressed, name)
end

Base.propertynames(::DressedStates, private::Bool = false) =
    (:tracking, :components, :states, :energies, :keys)


Base.show(io::IO, dressed::DressedStates) = print(
    io,
    "DressedStates(",
    length(dressed.energies),
    " states, components=",
    dressed.components,
    ", keys=",
    length(dressed.keys),
    ")",
)

function Base.show(io::IO, ::MIME"text/plain", dressed::DressedStates)
    return summary_lines(
        io,
        dressed,
        [
            "DressedStates · $(length(dressed.energies)) states in energy order",
            "  physical component order: $(dressed.components)",
            "  keys: $(length(dressed.keys)) · saved samples: $(dressed.tracking.saved_steps) · processed samples: $(dressed.tracking.total_steps)",
            dressed.tracking.tracked ? "  tracking diagnostics available" :
            "  tracking unavailable; normal indexing remains available",
        ],
    )
end
