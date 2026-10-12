"""
    DRAG(sigma; beta=0)

Complex baseband envelope `G(t) + im*beta*G′(t)`, where `G` is `GaussianZero(sigma)`
and the derivative is with respect to physical time. `beta` has units of time and
may have either sign. Normalization applies to the in-phase Gaussian, not the
complex magnitude. The derivative quadrature can have nonzero endpoints.
Use `Pulse(DRAG(...); carrier=IQCarrier(...), ...)` for a real laboratory-frame drive.
"""
struct DRAG{S<:Real,B<:Real} <: AbstractEnvelope
    sigma::S
    beta::B
    function DRAG(sigma::S, beta::B) where {S<:Real,B<:Real}
        _require_positive_width("sigma", sigma)
        _require_finite_real("beta", beta)
        new{S,B}(sigma, beta)
    end
end
DRAG(sigma; beta=0) = DRAG(sigma, beta)
