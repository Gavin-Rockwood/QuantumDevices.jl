"""
    AvoidedCrossingFit

Fit of `sqrt(minimum_gap^2 + slope^2*(parameter-center)^2)`. Fields `center`,
`minimum_gap`, and `slope` use the original data units. `parameters` and `gaps`
are sorted samples; `fitted_gaps`, `residuals` (fit minus data), and `rms_residual`
describe fit quality. `fit` is the native LsqFit result in normalized coordinates:
`x=(parameter-x_offset)/x_scale`, `y=gap/gap_scale`. These three normalization
fields allow interpretation of its parameters and Jacobian. Residuals are
diagnostics, not a guarantee that an isolated two-state model is appropriate.
"""
struct AvoidedCrossingFit{F}
    center::Float64
    minimum_gap::Float64
    slope::Float64
    parameters::Vector{Float64}
    gaps::Vector{Float64}
    fitted_gaps::Vector{Float64}
    residuals::Vector{Float64}
    rms_residual::Float64
    fit::F
    x_offset::Float64
    x_scale::Float64
    gap_scale::Float64
end

function _spectral_parameters(parameters)
    parameters isa AbstractVector{<:Real} ||
        throw(ArgumentError("Parameters must be a real-valued vector"))
    length(parameters) >= 5 || throw(ArgumentError("At least five samples are required"))
    x = Float64.(parameters)
    all(isfinite, x) || throw(ArgumentError("Parameters must be finite"))
    allunique(x) || throw(ArgumentError("Parameters must be distinct"))
    order = sortperm(x)
    return x[order], order
end

_avoided_crossing(x, p) = sqrt.(p[2]^2 .+ p[3]^2 .* (x .- p[1]).^2)

"""
    fit_avoided_crossing(parameters, gaps; fit_kwargs=(;))

Fit an isolated avoided crossing to
`sqrt(minimum_gap^2 + slope^2*(parameter-center)^2)` and return an
[`AvoidedCrossingFit`](@ref). Inputs require at least five distinct finite real
parameters and matching finite, nonnegative gaps. Samples are sorted ascending.
The sampled minimum must be interior and the data must be nonflat.

Parameter and gap scales are normalized internally. `fit_kwargs` forwards
LsqFit solver options (e.g. `maxIter`); `lower`, `upper`, and `inplace` are reserved
to preserve the nonnegative-gap/slope model. Nonconvergence, nonfinite results,
an unbracketed fitted center, or a rank-deficient Jacobian raise errors.
"""
function fit_avoided_crossing(parameters, gaps; fit_kwargs=(;))
    x, order = _spectral_parameters(parameters)
    gaps isa AbstractVector{<:Real} || throw(ArgumentError("Gaps must be a real-valued vector"))
    length(gaps) == length(x) || throw(DimensionMismatch("Parameter and gap counts differ"))
    y = Float64.(gaps[order])
    all(v -> isfinite(v) && v >= 0, y) || throw(ArgumentError("Gaps must be finite and nonnegative"))
    ymax, ymin = maximum(y), minimum(y)
    ymax - ymin > sqrt(eps(Float64)) * ymax || throw(ArgumentError("Gap data are flat or unresolved"))
    i = argmin(y)
    1 < i < length(x) || throw(ArgumentError("Sampled minimum must be bracketed by the parameter range"))
    any(k -> k in (:lower, :upper, :inplace), keys(fit_kwargs)) &&
        throw(ArgumentError("Fit bounds and inplace mode are reserved"))
    offset = x[i]
    scale = last(x) - first(x)
    isfinite(scale) && scale > 0 || throw(ArgumentError("Parameter range cannot be normalized"))
    xn, yn = (x .- offset) ./ scale, y ./ ymax
    slope = maximum(sqrt(max(yn[j]^2 - yn[i]^2, 0)) / abs(xn[j])
        for j in (1, length(x)))
    p0 = [0.0, max(yn[i], sqrt(eps(Float64))), slope]
    fit = LsqFit.curve_fit(_avoided_crossing, xn, yn, p0;
        lower=[-Inf, 0.0, 0.0], upper=fill(Inf, 3), fit_kwargs...)
    fit.converged || error("Avoided-crossing fit did not converge")
    all(isfinite, fit.param) && all(isfinite, fit.jacobian) ||
        error("Avoided-crossing fit produced nonfinite results")
    rank(fit.jacobian) == 3 || error("Avoided-crossing fit has a rank-deficient Jacobian")
    center = offset + scale * fit.param[1]
    first(x) < center < last(x) || error("Fitted minimum lies outside the sampled interval")
    gap, fitted_slope = ymax * fit.param[2], ymax / scale * fit.param[3]
    all(isfinite, (center, gap, fitted_slope)) || error("Fit cannot be represented in the original units")
    fitted = ymax .* _avoided_crossing(xn, fit.param)
    residuals = fitted .- y
    rms = norm(residuals) / sqrt(length(y))
    return AvoidedCrossingFit(center, gap, fitted_slope, x, y, fitted, residuals,
        rms, fit, offset, scale, ymax)
end

function Base.show(io::IO, ::MIME"text/plain", result::AvoidedCrossingFit)
    print(io, "AvoidedCrossingFit\n  Center: ", result.center,
        "\n  Minimum gap: ", result.minimum_gap, "\n  Slope: ", result.slope,
        "\n  RMS residual: ", result.rms_residual)
end
