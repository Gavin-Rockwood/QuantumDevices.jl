"""
    Slepian(time_bandwidth=2.5;duration,samples=129)

Peak-normalized DPSS window. Odd sample count >=3; time_bandwidth is dimensionless. Weights are cached in the normalized shape.
Returns an Envelope with a physical reference duration and normalized shape payload.
"""
struct Slepian{W<:Real,T<:Tuple}
    time_bandwidth::W
    samples::Int
    weights::T
    function Slepian(::Val{:raw},nw::W, samples::Integer) where {W<:Real}
        _require_finite_real("time_bandwidth", nw)
        n = Int(samples)
        n >= 3 && isodd(n) || throw(ArgumentError("samples must be an odd integer ≥3"))
        0 < nw < n/2 || throw(ArgumentError("time_bandwidth must lie strictly between 0 and samples/2"))
        w = Float64(nw)/n
        diagonal = [((n-1-2j)/2)^2 * cospi(2w) for j = 0:(n-1)]
        offdiagonal = [j*(n-j)/2 for j = 1:(n-1)]
        weights = vec(eigen(SymTridiagonal(diagonal, offdiagonal), n:n).vectors)
        weights[n÷2+1] < 0 && (weights .*= -1)
        weights ./= maximum(weights)
        new{W,typeof(Tuple(weights))}(nw, n, Tuple(weights))
    end
end
function Slepian(nw=2.5;duration=nothing,samples=129,rise=_default_rise(duration),fall=_default_fall(duration))
    return _wrap_shape(Slepian(Val(:raw),nw,samples),duration,rise,fall)
end
_physical_parameters(s::Slepian,d)=(time_bandwidth=s.time_bandwidth,samples=s.samples)
_shape_from_physical(::Slepian,p,d)=Slepian(Val(:raw),p.time_bandwidth,p.samples)
