"""
    Sech(width; center=nothing)

Unit-peak `sech((t-center)/width)`. Width uses physical time units and the default
center is half the pulse duration. Truncation retains nonzero endpoint values.
"""
struct Sech{W<:Real,C} <: AbstractEnvelope
    width::W
    center::C
    function Sech(width::W, center::C) where {W<:Real,C}
        _require_positive_width("width", width)
        center === nothing || _require_finite_real("center", center)
        new{W,C}(width, center)
    end
end
Sech(width; center=nothing) = Sech(width, center)
