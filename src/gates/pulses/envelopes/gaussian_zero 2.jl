"""
    GaussianZero(sigma)

Centered, unit-peak Gaussian with the endpoint baseline subtracted. `sigma` uses
physical time units. Both endpoints are exactly zero; unlike `Gaussian`, the
center is always half the pulse duration.
"""
struct GaussianZero{S<:Real} <: AbstractEnvelope
    sigma::S
    function GaussianZero(sigma::S) where {S<:Real}
        _require_positive_width("sigma", sigma)
        new{S}(sigma)
    end
end
