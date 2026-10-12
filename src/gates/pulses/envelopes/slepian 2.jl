"""
    Slepian(time_bandwidth=2.5; samples=129)

Leading discrete prolate spheroidal sequence (DPSS), normalized to unit peak and
linearly interpolated over the pulse duration. `time_bandwidth` is the standard
dimensionless DPSS time-half-bandwidth product NW, with `0 < NW < samples/2`.
`samples` must be an odd integer ≥3. Endpoints are not shifted to zero.
The Float64 taper is cached at construction and rebuilt by `setpath` updates.
Only `time_bandwidth` is exposed for calibration; `samples` is configuration.
"""
struct Slepian{W<:Real} <: AbstractEnvelope
    time_bandwidth::W
    samples::Int
    weights::Vector{Float64}
    function Slepian(nw::W, samples::Integer) where {W<:Real}
        _require_finite_real("time_bandwidth", nw)
        n = Int(samples)
        n >= 3 && isodd(n) || throw(ArgumentError("samples must be an odd integer ≥3"))
        0 < nw < n/2 || throw(ArgumentError("time_bandwidth must lie strictly between 0 and samples/2"))
        w = Float64(nw)/n
        diagonal = [((n-1-2j)/2)^2 * cospi(2w) for j in 0:n-1]
        offdiagonal = [j*(n-j)/2 for j in 1:n-1]
        weights = vec(eigen(SymTridiagonal(diagonal, offdiagonal), n:n).vectors)
        weights[n÷2+1] < 0 && (weights .*= -1)
        weights ./= maximum(weights)
        new{W}(nw, n, weights)
    end
end
Slepian(nw=2.5; samples=129) = Slepian(nw, samples)
