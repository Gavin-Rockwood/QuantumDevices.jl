"""
    Gaussian(sigma; center=nothing)

Unit-peak Gaussian with width `sigma` in physical time units. Default center is
half the pulse duration. Truncation does not subtract its nonzero endpoint values.
"""
struct Gaussian{S<:Real,C} <: AbstractEnvelope
    sigma::S
    center::C
    function Gaussian(sigma::S, center::C) where {S<:Real,C}
        _require_finite_real("sigma", sigma) > 0 || throw(ArgumentError("sigma must be positive"))
        center === nothing || _require_finite_real("center", center)
        new{S,C}(sigma, center)
    end
end
Gaussian(sigma; center=nothing) = Gaussian(sigma, center)
