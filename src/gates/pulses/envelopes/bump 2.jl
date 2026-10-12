"""
    Bump(k=2; center=nothing)

Smooth compact-support, unit-peak envelope `exp(k*x^2/(x^2-1))` for
`abs(x) < 1`, zero otherwise, where `x = (t-center)/(duration/2)`.
The positive dimensionless `k` controls sharpness; larger values narrow the peak.
The center defaults to half the pulse duration and uses pulse-local time units.
Changing `Pulse.duration` automatically rescales the width. A shifted center can
truncate the bump at the pulse window boundaries.
"""
struct Bump{K<:Real,C} <: AbstractEnvelope
    k::K
    center::C
    function Bump(k::K, center::C) where {K<:Real,C}
        _require_positive_width("k", k)
        center === nothing || _require_finite_real("center", center)
        new{K,C}(k, center)
    end
end
Bump(k=2; center=nothing) = Bump(k, center)
