# Bounded v1 physical-time implementation, never used by new constructors.
module _LegacyPhysical
using LinearAlgebra, Interpolations
using SpecialFunctions: erf
using ..QuantumDevices: _require_finite, _require_finite_real, _require_positive_width
abstract type AbstractEnvelope end
function envelope_value end
validate_envelope(e::AbstractEnvelope,d)=(_require_finite("Envelope value",envelope_value(e,0.,d));_require_finite("Envelope value",envelope_value(e,d,d));nothing)
function envelope_tstops end
_validated_envelope_tstops(e,d)=collect(envelope_tstops(e,d))
"""`Constant()` is the unit envelope throughout the pulse window."""
struct Constant <: AbstractEnvelope end

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
Gaussian(sigma; center = nothing) = Gaussian(sigma, center)

"""`SineSquared()` is one unit-height lobe, `sinpi(t / duration)^2`."""
struct SineSquared <: AbstractEnvelope end

"""
    Envelope(f, parameters=(;))

Custom dimensionless envelope with explicit named parameters. `f(p, t, duration)`
receives local pulse time and duration. Parameters are exposed recursively to
`getpath`, `setpath`, and calibration. Callable persistence uses JLD2 and requires
its definition to be available when loading; built-in shapes have no such requirement.
"""
struct Envelope{F,P<:NamedTuple} <: AbstractEnvelope
    callable::F
    parameters::P
    function Envelope(f::F, p::P) where {F,P<:NamedTuple}
        applicable(f, p, 0.0, 1.0) || throw(ArgumentError("Envelope must accept f(p, t, duration)"))
        new{F,P}(f, p)
    end
end
Envelope(f) = Envelope(f, (;))

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
Sech(width; center = nothing) = Sech(width, center)

"""
    Cosine()

Single unit-peak cosine lobe centered on half the pulse duration:
`cospi((t-duration/2)/duration)`, equivalently `sinpi(t/duration)`.
Endpoints are exactly zero. This is not the squared lobe of `SineSquared`.
"""
struct Cosine <: AbstractEnvelope end

"""
    Blackman()

Symmetric unit-peak Blackman window:
`0.42 - 0.5cospi(2t/duration) + 0.08cospi(4t/duration)`.
Both endpoints are exactly zero.
"""
struct Blackman <: AbstractEnvelope end

"""
    Bump(width; k=2, center=nothing)

Compact-support unit-peak envelope with physical support half-width `width`.
The dimensionless positive `k` controls sharpness. Center defaults to duration/2;
changing duration never stretches the support width.
"""
struct Bump{W<:Real,K<:Real,C} <: AbstractEnvelope
    width::W
    k::K
    center::C
    function Bump(width::W, k::K, center::C) where {W<:Real,K<:Real,C}
        _require_positive_width("width", width)
        _require_positive_width("k", k)
        center===nothing || _require_finite_real("center", center)
        new{W,K,C}(width, k, center)
    end
end
Bump(width; k = 2, center = nothing) = Bump(width, k, center)

"""
    DRAG(sigma; beta=0)

Complex baseband envelope `G(t) + im*beta*G′(t)`, where `G` is `GaussianZero(sigma)`
and the derivative is with respect to physical time. `beta` has units of time and
may have either sign. Normalization applies to the in-phase Gaussian, not the
complex magnitude. The derivative quadrature can have nonzero endpoints.
Use `Pulse(DRAG(...), frequency; ...)` for a real laboratory-frame drive.
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
DRAG(sigma; beta = 0) = DRAG(sigma, beta)

"""
    Slepian(time_bandwidth=2.5; samples=129)

Leading discrete prolate spheroidal sequence (DPSS), normalized to unit peak and
linearly interpolated over the pulse duration. `time_bandwidth` is the standard
dimensionless DPSS time-half-bandwidth product NW, with `0 < NW < samples/2`.
`samples` must be an odd integer ≥3. Endpoints are not shifted to zero.
The Float64 taper is cached at construction and rebuilt by `setpath` updates.
`time_bandwidth` is a numeric calibration candidate; `samples` is discoverable configuration.
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
        diagonal = [((n-1-2j)/2)^2 * cospi(2w) for j = 0:(n-1)]
        offdiagonal = [j*(n-j)/2 for j = 1:(n-1)]
        weights = vec(eigen(SymTridiagonal(diagonal, offdiagonal), n:n).vectors)
        weights[n÷2+1] < 0 && (weights .*= -1)
        weights ./= maximum(weights)
        new{W}(nw, n, weights)
    end
end
Slepian(nw = 2.5; samples = 129) = Slepian(nw, samples)

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
        _require_finite_real("ramp_time", ramp_time) >= 0 ||
            throw(ArgumentError("ramp_time must be nonnegative"))
        new{S,R}(sigma, ramp_time)
    end
end
GaussianSquare(sigma; ramp_time = 2sigma) = GaussianSquare(sigma, ramp_time)

"""
    ErfSquare(sigma; ramp_time=2sigma)

Unit-height flattop with symmetric, endpoint-normalized error-function edges.
`sigma` and `ramp_time` use physical time units; both ramps must fit in the pulse
window. Positive ramp times give exact zero endpoints. `ramp_time=0` gives an
abrupt unit-height square. The plateau and each ramp join have value one.
"""
struct ErfSquare{S<:Real,R<:Real} <: AbstractEnvelope
    sigma::S
    ramp_time::R
    function ErfSquare(sigma::S, ramp_time::R) where {S<:Real,R<:Real}
        _require_positive_width("sigma", sigma)
        _require_finite_real("ramp_time", ramp_time) >= 0 ||
            throw(ArgumentError("ramp_time must be nonnegative"))
        new{S,R}(sigma, ramp_time)
    end
end
ErfSquare(sigma; ramp_time = 2sigma) = ErfSquare(sigma, ramp_time)

"""
    InterpolatedEnvelope(times, values; interpolation=:linear)

Piecewise linear physical-time samples, starting at zero. Nodes must be finite,
strictly increasing, and cover exactly the enclosing duration. No extrapolation
or duration rescaling is performed. Real and complex samples are supported.
"""
struct InterpolatedEnvelope{T,V,I} <: AbstractEnvelope
    times::T
    values::V
    interpolation::Symbol
    interpolant::I
    function InterpolatedEnvelope(times, values; interpolation = :linear)
        interpolation === :linear || throw(ArgumentError("Only :linear interpolation is supported"))
        times isa AbstractVector && values isa AbstractVector || throw(ArgumentError("Nodes must be vectors"))
        length(times) == length(values) >= 2 ||
            throw(ArgumentError("Need equally sized node vectors with at least two nodes"))
        all(t -> t isa Real && isfinite(t), times) || throw(ArgumentError("Times must be finite real nodes"))
        first(times) == 0 && all(>(0), diff(times)) ||
            throw(ArgumentError("Times must start at zero and strictly increase"))
        foreach(v -> _require_finite("sample", v), values)
        ts, vs = collect(times), collect(values)
        itp = interpolate((ts,), vs, Gridded(Linear()))
        new{typeof(ts),typeof(vs),typeof(itp)}(ts, vs, interpolation, itp)
    end
end
validate_envelope(e::InterpolatedEnvelope, d) =
    (last(e.times) == d || throw(ArgumentError("Interpolation domain must equal duration")); nothing)
envelope_value(e::InterpolatedEnvelope, t, d) = e.interpolant(t)
envelope_tstops(e::InterpolatedEnvelope, d) = e.times[2:(end-1)]

struct Segment{E<:AbstractEnvelope,D<:Real}
    shape::E
    duration::D
    function Segment(shape::E; duration) where {E<:AbstractEnvelope}
        _require_positive_width("segment duration", duration)
        validate_envelope(shape, duration)
        new{E,typeof(duration)}(shape, duration)
    end
end
struct EnvelopeSequence{S<:Tuple} <: AbstractEnvelope
    segments::S
    function EnvelopeSequence(segments::Segment...)
        isempty(segments) && throw(ArgumentError("EnvelopeSequence needs segments"))
        isfinite(sum(s.duration for s in segments)) || throw(ArgumentError("Duration overflow"))
        new{typeof(segments)}(segments)
    end
end
_envelope_duration(e::EnvelopeSequence) = sum(s.duration for s in e.segments)
_envelope_duration(::AbstractEnvelope) = nothing
validate_envelope(e::EnvelopeSequence, d) = (
    d == _envelope_duration(e) ||
        throw(ArgumentError("Sequence duration must equal sum of segment durations"));
    nothing
)
function envelope_value(e::EnvelopeSequence, t, d)
    start = zero(t)
    for (i, s) in enumerate(e.segments)
        stop = start+s.duration
        if t < stop || i == length(e.segments)
            return s.shape(t == stop ? s.duration : t-start, s.duration)
        end
        start = stop
    end
end
function envelope_tstops(e::EnvelopeSequence, d)
    result = Real[]
    start = zero(d)
    for s in e.segments
        append!(result, start .+ _validated_envelope_tstops(s.shape, s.duration))
        start += s.duration
        push!(result, start)
    end
    return sort!(unique(filter(t -> 0<t<d, result)))
end

const _RAMP_NAMES = (:sine_squared, :linear, :smoothstep)
struct Ramp{S} <: AbstractEnvelope
    source::S
    rising::Bool
    function Ramp(source, rising::Bool)
        source isa AbstractEnvelope || source in _RAMP_NAMES || throw(ArgumentError("Unknown ramp"))
        new{typeof(source)}(source, rising)
    end
end
function envelope_value(r::Ramp, t, d)
    x = r.rising ? t/d : 1-t/d
    r.source === :linear && return x
    r.source === :smoothstep && return x^2*(3-2x)
    r.source === :sine_squared && return sinpi(x/2)^2
    # Envelope sources receive physical ramp time, without normalized-width sampling.
    return r.source(r.rising ? t : d-t, d)
end
validate_envelope(r::Ramp, d) = r.source isa AbstractEnvelope ? validate_envelope(r.source, d) : nothing
envelope_tstops(r::Ramp, d) =
    r.source isa AbstractEnvelope ?
    sort(r.rising ? _validated_envelope_tstops(r.source, d) : d .- _validated_envelope_tstops(r.source, d)) :
    Real[]
"""
    RampedFlattop(ramp_time; ramp=:sine_squared, rise_time=ramp_time,
                  fall_time=ramp_time, ramp_up=ramp, ramp_down=ramp)

Physical-time rise, constant plateau, and fall composition. Arbitrary envelope
sources are evaluated over each ramp's physical duration; they are not normalized.
"""
struct RampedFlattop{R,F,U,D} <: AbstractEnvelope
    rise_time::R
    fall_time::F
    ramp_up::U
    ramp_down::D
    function RampedFlattop(rise, fall, up, down)
        _require_finite_real("rise_time", rise)>=0 || throw(ArgumentError("Negative rise"))
        _require_finite_real("fall_time", fall)>=0 || throw(ArgumentError("Negative fall"))
        Ramp(up, true)
        Ramp(down, false)
        new{typeof(rise),typeof(fall),typeof(up),typeof(down)}(rise, fall, up, down)
    end
end
RampedFlattop(r; ramp = :sine_squared, rise_time = r, fall_time = r, ramp_up = ramp, ramp_down = ramp) =
    RampedFlattop(rise_time, fall_time, ramp_up, ramp_down)
function _flattop_sequence(e::RampedFlattop, d)
    e.rise_time+e.fall_time<=d || throw(ArgumentError("Ramps must fit duration"))
    segments = Segment[]
    e.rise_time>0 && push!(segments, Segment(Ramp(e.ramp_up, true); duration = e.rise_time))
    plateau=d-e.rise_time-e.fall_time
    plateau>0 && push!(segments, Segment(Constant(); duration = plateau))
    e.fall_time>0 && push!(segments, Segment(Ramp(e.ramp_down, false); duration = e.fall_time))
    return EnvelopeSequence(segments...)
end
validate_envelope(e::RampedFlattop, d) = (_flattop_sequence(e, d); nothing)
# Evaluate the three segments directly, keeping sequence boundary and physical-time rules.
function envelope_value(e::RampedFlattop, t, d)
    e.rise_time+e.fall_time<=d || throw(ArgumentError("Ramps must fit duration"))
    _require_positive_width("duration", d)
    plateau = d-e.rise_time-e.fall_time
    start = zero(t)
    if e.rise_time > 0
        stop = start+e.rise_time
        if t < stop || (plateau <= 0 && e.fall_time == 0)
            return Ramp(e.ramp_up, true)(t == stop ? e.rise_time : t-start, e.rise_time)
        end
        start = stop
    end
    if plateau > 0
        stop = start+plateau
        if t < stop || e.fall_time == 0
            return one(t == stop ? plateau : t-start)
        end
        start = stop
    end
    stop = start+e.fall_time
    return Ramp(e.ramp_down, false)(t == stop ? e.fall_time : t-start, e.fall_time)
end
envelope_tstops(e::RampedFlattop, d) = envelope_tstops(_flattop_sequence(e, d), d)

envelope_value(::Constant, t, duration) = one(t)
function envelope_value(shape::Gaussian, t, duration)
    center = shape.center === nothing ? duration / 2 : shape.center
    return exp(-(t - center)^2 / (2shape.sigma^2))
end
function validate_envelope(shape::Gaussian, duration)
    center = shape.center === nothing ? duration / 2 : shape.center
    0 <= center <= duration || throw(ArgumentError("Gaussian center must lie inside the pulse window"))
    return nothing
end
envelope_value(::SineSquared, t, duration) = t == 0 || t == duration ? zero(t) : sinpi(t / duration)^2
envelope_value(shape::Envelope, t, duration) = shape.callable(shape.parameters, t, duration)
(shape::AbstractEnvelope)(t, duration) = _require_finite("Envelope value", envelope_value(shape, t, duration))

envelope_tstops(::AbstractEnvelope, duration) = Float64[]
envelope_tstops(shape::Union{GaussianSquare,ErfSquare}, duration) =
    iszero(shape.ramp_time) ? Float64[] : sort!(unique([shape.ramp_time, duration-shape.ramp_time]))
envelope_tstops(shape::Slepian, duration) = [duration*j/(shape.samples-1) for j = 1:(shape.samples-2)]

# expm1 avoids cancellation when a Gaussian is much wider than the pulse.
function _zero_gaussian(sigma, t, duration)
    (t == 0 || t == duration) && return zero(float(t))
    x = (t-duration/2)/sigma
    h = duration/(2sigma)
    !isfinite(h^2) && return exp(-x^2/2)
    denominator = -expm1(-h^2/2)
    if iszero(denominator)
        # Limit as sigma/duration → infinity.
        u = t/duration
        return 4u*(1-u)
    end
    return exp(-x^2/2) * (-expm1((x^2-h^2)/2)) / denominator
end
function _zero_gaussian_derivative(sigma, t, duration)
    h = duration/(2sigma)
    denominator = -expm1(-h^2/2)
    iszero(denominator) && return 4*(1-2t/duration)/duration
    x = (t-duration/2)/sigma
    !isfinite(x) && return zero(float(t))
    return -x*exp(-x^2/2)/(sigma*denominator)
end
envelope_value(shape::GaussianZero, t, duration) = _zero_gaussian(shape.sigma, t, duration)
envelope_value(shape::DRAG, t, duration) = complex(
    _zero_gaussian(shape.sigma, t, duration),
    shape.beta*_zero_gaussian_derivative(shape.sigma, t, duration),
)
envelope_value(::Cosine, t, duration) = t == 0 || t == duration ? zero(float(t)) : sinpi(t/duration)
envelope_value(::Blackman, t, duration) =
    t == 0 || t == duration ? zero(float(t)) : 0.42-0.5cospi(2t/duration)+0.08cospi(4t/duration)
function envelope_value(shape::Sech, t, duration)
    center = shape.center === nothing ? duration/2 : shape.center
    return inv(cosh((t-center)/shape.width))
end
function validate_envelope(shape::Sech, duration)
    center = shape.center === nothing ? duration/2 : shape.center
    0 <= center <= duration || throw(ArgumentError("Sech center must lie inside the pulse window"))
    return nothing
end
function envelope_value(shape::Bump, t, duration)
    center = shape.center === nothing ? duration/2 : shape.center
    x = (t-center)/shape.width
    abs(x) >= 1 && return zero(float(x))
    return exp(shape.k * (x^2/(x^2-1)))
end
function validate_envelope(shape::Bump, duration)
    center = shape.center === nothing ? duration/2 : shape.center
    0 <= center <= duration || throw(ArgumentError("Bump center must lie inside the pulse window"))
    return nothing
end
function validate_envelope(shape::Union{GaussianSquare,ErfSquare}, duration)
    2shape.ramp_time <= duration || throw(ArgumentError("Both ramps must fit inside the pulse duration"))
    return nothing
end
function envelope_value(shape::GaussianSquare, t, duration)
    r = shape.ramp_time
    iszero(r) && return one(float(t))
    r <= t <= duration-r && return one(float(t))
    distance = min(t, duration-t)
    distance >= r && return one(float(t))
    # Left half of a centered, baseline-subtracted Gaussian lobe.
    return _zero_gaussian(shape.sigma, distance, 2r)
end
function envelope_value(shape::ErfSquare, t, duration)
    r = shape.ramp_time
    iszero(r) && return one(float(t))
    r <= t <= duration-r && return one(float(t))
    distance = min(t, duration-t)
    distance <= 0 && return zero(float(t))
    distance >= r && return one(float(t))
    h = r/(2sqrt(2)*shape.sigma)
    denominator = erf(h)
    iszero(denominator) && return distance/r
    return (erf((distance-r/2)/(sqrt(2)*shape.sigma))+denominator)/(2denominator)
end
function envelope_value(shape::Slepian, t, duration)
    t <= 0 && return first(shape.weights)
    t >= duration && return last(shape.weights)
    x = (shape.samples-1)*(t/duration)
    index = min(floor(Int, x)+1, shape.samples-1)
    fraction = x-(index-1)
    return muladd(fraction, shape.weights[index+1]-shape.weights[index], shape.weights[index])
end

end
# Only these historical constructors can be selected by a persisted record.
const _LEGACY_SHAPE_TYPES=(Constant=_LegacyPhysical.Constant,Gaussian=_LegacyPhysical.Gaussian,
    GaussianZero=_LegacyPhysical.GaussianZero,SineSquared=_LegacyPhysical.SineSquared,
    Cosine=_LegacyPhysical.Cosine,Blackman=_LegacyPhysical.Blackman,Sech=_LegacyPhysical.Sech,
    Bump=_LegacyPhysical.Bump,DRAG=_LegacyPhysical.DRAG,Slepian=_LegacyPhysical.Slepian,
    GaussianSquare=_LegacyPhysical.GaussianSquare,ErfSquare=_LegacyPhysical.ErfSquare,
    RampedFlattop=_LegacyPhysical.RampedFlattop,Ramp=_LegacyPhysical.Ramp,
    InterpolatedEnvelope=_LegacyPhysical.InterpolatedEnvelope,Envelope=_LegacyPhysical.Envelope,
    EnvelopeSequence=_LegacyPhysical.EnvelopeSequence)
_legacy_field_names(T)=filter(k->!(k in (:interpolant,:weights)),fieldnames(T))
struct LegacyEnvelope{F,D,S} <: AbstractEnvelope
    implementation::Symbol
    fields::F
    duration::D
    legacy_shape::S
    LegacyEnvelope(::Val{:validated},kind,f::F,d::D,s::S) where {F,D,S}=new{F,D,S}(kind,f,d,s)
end
_legacy_input(e::LegacyEnvelope)=getfield(e,:legacy_shape)
function _legacy_input(e::EnvelopeSequence)
    all(s->s.scale==1,e.segments) || throw(ArgumentError("Historical ramp sequences do not support scaled pieces"))
    return _LegacyPhysical.EnvelopeSequence((_legacy_input(s) for s in e.segments)...)
end
_legacy_input(s::Segment)=_LegacyPhysical.Segment(_legacy_input(s.shape);duration=s.duration)
_legacy_input(x::Tuple)=map(_legacy_input,x)
_legacy_input(x::AbstractEnvelope)=throw(ArgumentError("A legacy physical ramp requires a legacy physical source"))
_legacy_input(x)=x
_legacy_field(e::_LegacyPhysical.AbstractEnvelope)=_legacy_from_source(e,nothing)
_legacy_field(s::_LegacyPhysical.Segment)=Segment(_legacy_from_source(s.shape,s.duration);duration=s.duration)
_legacy_field(x::Tuple)=map(_legacy_field,x)
_legacy_field(x)=x
function _legacy_from_source(source,duration)
    kind=nameof(typeof(source))
    haskey(_LEGACY_SHAPE_TYPES,kind) && source isa _LEGACY_SHAPE_TYPES[kind] ||
        throw(ArgumentError("Unsupported historical shape $kind"))
    names=_legacy_field_names(typeof(source))
    fields=NamedTuple{Tuple(names)}(Tuple(_legacy_field(getproperty(source,k)) for k in names))
    return LegacyEnvelope(Val(:validated),kind,fields,duration,source)
end
function LegacyEnvelope(kind::Symbol,fields::NamedTuple,duration)
    haskey(_LEGACY_SHAPE_TYPES,kind) || throw(ArgumentError("Unsupported historical shape $kind"))
    T=_LEGACY_SHAPE_TYPES[kind]
    names=_legacy_field_names(T)
    Set(keys(fields))==Set(names) || throw(ArgumentError("Historical shape field schema mismatch for $kind"))
    f=NamedTuple{Tuple(names)}(Tuple(_legacy_input(getproperty(fields,k)) for k in names))
    source=if kind===:InterpolatedEnvelope
        T(f.times,f.values;interpolation=f.interpolation)
    elseif kind===:Slepian
        T(f.time_bandwidth;samples=f.samples)
    elseif kind===:EnvelopeSequence
        T(f.segments...)
    else
        T(values(f)...)
    end
    e=_legacy_from_source(source,duration)
    if duration!==nothing
        _require_positive_width("duration",duration)
        validate_envelope(e,duration)
    end
    return e
end
_legacy_value(e::LegacyEnvelope)=getfield(e,:legacy_shape)
envelope_value(e::LegacyEnvelope,t,d)=_LegacyPhysical.envelope_value(_legacy_value(e),t,d)
validate_envelope(e::LegacyEnvelope,d)=_LegacyPhysical.validate_envelope(_legacy_value(e),d)
envelope_tstops(e::LegacyEnvelope,d)=_LegacyPhysical.envelope_tstops(_legacy_value(e),d)
_envelope_duration(e::LegacyEnvelope)=e.duration
function _stretch_envelope(e::LegacyEnvelope,d)
    _require_positive_width("duration",d)
    validate_envelope(e,d)
    return LegacyEnvelope(Val(:validated),e.implementation,e.fields,d,_legacy_value(e))
end
function Base.getproperty(e::LegacyEnvelope,k::Symbol)
    k in (:implementation,:fields,:duration,:legacy_shape) && return getfield(e,k)
    k in (:rise,:fall) && return nothing
    return getproperty(getfield(e,:fields),k)
end
Base.propertynames(e::LegacyEnvelope,private::Bool=false)=(:implementation,:fields,:duration,:rise,:fall,keys(e.fields)...,(private ? (:legacy_shape,) : ())...)

# An explicitly duration-bound custom subtype retains its documented physical-time contract.
struct PhysicalEnvelope{E,D} <: AbstractEnvelope
    source::E
    duration::D
end
envelope_value(e::PhysicalEnvelope,t,d)=envelope_value(e.source,t,d)
validate_envelope(e::PhysicalEnvelope,d)=validate_envelope(e.source,d)
envelope_tstops(e::PhysicalEnvelope,d)=envelope_tstops(e.source,d)
_envelope_duration(e::PhysicalEnvelope)=e.duration
_stretch_envelope(e::PhysicalEnvelope,d)=PhysicalEnvelope(e.source,d)
_stretch_envelope(e::AbstractEnvelope,d)=PhysicalEnvelope(e,d)
