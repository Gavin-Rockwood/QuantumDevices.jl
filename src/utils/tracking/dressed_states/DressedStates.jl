"""Final energy-ordered construction-basis states and energies with physical tuple order."""
struct DressedStates{R <: TrackingResult, C <: Tuple}
    tracking::R
    components::C
    function DressedStates(tracking::R, components::C) where {R <: TrackingResult, C <: Tuple}
        length(tracking.saved_steps) == 1 && haskey(getfield(tracking, :_other), :energies) ||
            throw(ArgumentError("DressedStates requires a single saved sample with energies"))
        energies = getfield(tracking, :_other).energies
        all(e -> e isa Real && isfinite(e), energies) && issorted(energies) ||
            throw(ArgumentError("dressed energies must be finite real values in ascending order"))
        all(c -> c isa Symbol, components) && allunique(components) ||
            throw(ArgumentError("components must be distinct physical owner symbols"))
        all(
            key -> key isa Tuple &&
                 length(key) == length(components) &&
                 all(i -> i isa Integer && !(i isa Bool) && i >= 0, key),
            tracking.keys,
        ) || throw(
            ArgumentError(
                "occupation keys must contain one nonnegative integer level per physical component",
            ),
        )
        return new{R, C}(tracking, components)
    end
end
