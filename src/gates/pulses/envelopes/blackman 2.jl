"""
    Blackman()

Symmetric unit-peak Blackman window:
`0.42 - 0.5cospi(2t/duration) + 0.08cospi(4t/duration)`.
Both endpoints are exactly zero.
"""
struct Blackman <: AbstractEnvelope end
