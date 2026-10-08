"""
    Cosine()

Single unit-peak cosine lobe centered on half the pulse duration:
`cospi((t-duration/2)/duration)`, equivalently `sinpi(t/duration)`.
Endpoints are exactly zero. This is not the squared lobe of `SineSquared`.
"""
struct Cosine <: AbstractEnvelope end
