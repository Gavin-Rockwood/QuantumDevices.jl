_shape_value(::Constant,u)=one(u)
_shape_value(s::Gaussian,u)=exp(-(u-0.5)^2/(2s.sigma^2))
_shape_value(::SineSquared,u)=u==0 || u==1 ? zero(u) : sinpi(u)^2
(e::AbstractEnvelope)(t,d)=(_require_finite_real("time",t);_require_finite("Envelope value",envelope_value(e,t,d)))
(e::AbstractEnvelope)(t)=e(t,_envelope_duration(e))
envelope_tstops(::AbstractEnvelope,d)=()
_shape_tstops(s::Union{GaussianSquare,ErfSquare})=iszero(s.ramp_time) ? () : (s.ramp_time,1-s.ramp_time)
_shape_tstops(s::Slepian)=Tuple(j/(s.samples-1) for j=1:s.samples-2)
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
_shape_value(shape::GaussianZero, t) = _zero_gaussian(shape.sigma, t, 1)
_shape_value(shape::DRAG, t) = complex(
    _zero_gaussian(shape.sigma, t, 1),
    shape.beta*_zero_gaussian_derivative(shape.sigma, t, 1),
)
_shape_value(::Cosine, t) = t == 0 || t == 1 ? zero(float(t)) : sinpi(t/1)
_shape_value(::Blackman, t) =
    t == 0 || t == 1 ? zero(float(t)) : 0.42-0.5cospi(2t/1)+0.08cospi(4t/1)
function _shape_value(shape::Sech, t)
    center = shape.center === nothing ? 1/2 : shape.center
    return inv(cosh((t-center)/shape.width))
end
function _shape_value(shape::Bump, t)
    center = shape.center === nothing ? 1/2 : shape.center
    x = (t-center)/shape.width
    abs(x) >= 1 && return zero(float(x))
    return exp(shape.k * (x^2/(x^2-1)))
end
function _shape_value(shape::GaussianSquare, t)
    r = shape.ramp_time
    iszero(r) && return one(float(t))
    r <= t <= 1-r && return one(float(t))
    distance = min(t, 1-t)
    distance >= r && return one(float(t))
    # Left half of a centered, baseline-subtracted Gaussian lobe.
    return _zero_gaussian(shape.sigma, distance, 2r)
end
function _shape_value(shape::ErfSquare, t)
    r = shape.ramp_time
    iszero(r) && return one(float(t))
    r <= t <= 1-r && return one(float(t))
    distance = min(t, 1-t)
    distance <= 0 && return zero(float(t))
    distance >= r && return one(float(t))
    h = r/(2sqrt(2)*shape.sigma)
    denominator = erf(h)
    iszero(denominator) && return distance/r
    return (erf((distance-r/2)/(sqrt(2)*shape.sigma))+denominator)/(2denominator)
end
function _shape_value(shape::Slepian, t)
    t <= 0 && return first(shape.weights)
    t >= 1 && return last(shape.weights)
    x = (shape.samples-1)*(t/1)
    index = min(floor(Int, x)+1, shape.samples-1)
    fraction = x-(index-1)
    return muladd(fraction, shape.weights[index+1]-shape.weights[index], shape.weights[index])
end

# Validate normalized records explicitly, including raw data restored from IO.
_validate_shape(s)=nothing
function _validate_shape(s::Union{Gaussian,GaussianZero,DRAG})
    _require_positive_width("normalized sigma",s.sigma)
    s isa DRAG && _require_finite_real("normalized beta",s.beta)
    return nothing
end
function _validate_shape(s::Union{Sech,Bump})
    _require_positive_width("normalized width",s.width)
    _require_finite_real("normalized center",s.center)
    0<=s.center<=1 || throw(ArgumentError("center must lie inside duration"))
    s isa Bump && _require_positive_width("k",s.k)
    return nothing
end
function _validate_shape(s::Union{GaussianSquare,ErfSquare})
    _require_positive_width("normalized sigma",s.sigma)
    _require_finite_real("normalized ramp_time",s.ramp_time)
    0<=s.ramp_time<=.5 || throw(ArgumentError("Both ramps must fit duration"))
    return nothing
end
function _validate_shape(s::InterpolatedEnvelope)
    length(s.times)==length(s.values)>=2 || throw(ArgumentError("Need matching samples"))
    all(t->t isa Real && isfinite(t),s.times) && all(isfinite,s.values) || throw(ArgumentError("Nonfinite samples"))
    first(s.times)==0 && last(s.times)==1 && all(>(0),diff(collect(s.times))) || throw(ArgumentError("Invalid normalized sample nodes"))
    s.interpolation===:linear || throw(ArgumentError("Only linear interpolation is supported"))
    return nothing
end
function _validate_shape(s::RampedFlattop)
    _require_finite_real("rise_time",s.rise_time)>=0 || throw(ArgumentError("Negative rise"))
    _require_finite_real("fall_time",s.fall_time)>=0 || throw(ArgumentError("Negative fall"))
    _checked_time_add(s.rise_time,s.fall_time)<=1 || throw(ArgumentError("Ramps must fit duration"))
    for source in (s.ramp_up,s.ramp_down)
        if source isa Symbol
            _ramp_source(source)
        elseif source isa NamedTuple && keys(source)==(:shape,:rise,:fall)
            _validate_shape(source.shape)
            _normalized_interval(source.rise)
            _normalized_interval(source.fall)
        else
            throw(ArgumentError("Invalid normalized ramp source"))
        end
    end
    return nothing
end
