"""
    GaussianSquare(sigma; ramp_time=2sigma)

Unit-height flattop with symmetric, endpoint-normalized Gaussian edges.
`sigma` and `ramp_time` use physical time units; both ramps must fit in the pulse
window. Positive ramp times give exact zero endpoints. `ramp_time=0` gives an
abrupt unit-height square. The plateau and each ramp join have value one.
"""
struct GaussianSquare{S<:Real,R<:Real} <: AbstractEnvelope
    sigma::S
    ramp_time::R
    function GaussianSquare(sigma::S, ramp_time::R) where {S<:Real,R<:Real}
        _require_positive_width("sigma", sigma)
        _require_finite_real("ramp_time", ramp_time) >= 0 || throw(ArgumentError("ramp_time must be nonnegative"))
        new{S,R}(sigma, ramp_time)
    end
end
GaussianSquare(sigma; ramp_time=2sigma) = GaussianSquare(sigma, ramp_time)
