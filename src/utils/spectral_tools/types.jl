"""Portable numerical fit of sqrt(gap²+slope²*(x-center)²), with residual and rank diagnostics."""
struct AvoidedCrossingFit{N,D}
    center::Float64
    minimum_gap::Float64
    slope::Float64
    parameters::Vector{Float64}
    gaps::Vector{Float64}
    fitted_gaps::Vector{Float64}
    residuals::Vector{Float64}
    rms_residual::Float64
    normalization::N
    normalized_parameters::Vector{Float64}
    diagnostics::D
    function AvoidedCrossingFit(center,gap,slope,x,y,fitted,residuals,rms,normalization,normalized,diagnostics)
        length(x)==length(y)==length(fitted)==length(residuals)>=5 || throw(DimensionMismatch("fit sample counts differ"))
        all(isfinite,vcat(x,y,fitted,residuals,normalized,[center,gap,slope,rms])) || throw(ArgumentError("fit data must be finite"))
        length(normalized)==3 && normalization.x_scale>0 && normalization.gap_scale>0 || throw(ArgumentError("invalid normalized fit parameters"))
        first(x)<center<last(x) || throw(ArgumentError("fitted center must be bracketed"))
        issorted(x) && allunique(x) && gap>=0 && slope>=0 && rms>=0 || throw(ArgumentError("invalid avoided-crossing fit"))
        new{typeof(normalization),typeof(diagnostics)}(center,gap,slope,Float64.(x),Float64.(y),Float64.(fitted),Float64.(residuals),rms,normalization,Float64.(normalized),diagnostics)
    end
end
"""An isolated-crossing estimate owning its FloquetSweep and portable fit, not a drive closure."""
struct ResonanceResult{S,F,M,D}
    sweep::S
    fit::F
    reference_frequency::Union{Nothing,Float64}
    photon_order::Int
    metadata::M
    diagnostics::D
    function ResonanceResult(sweep,fit,reference_frequency,photon_order,metadata,diagnostics)
        length(sweep.tracking.keys)==2 && length(sweep.points)==length(fit.parameters) || throw(DimensionMismatch("resonance requires two complete tracked branches"))
        sweep.tracking.saved_steps==Tuple(eachindex(sweep.points)) || throw(ArgumentError("resonance requires all samples saved"))
        Float64.(sweep.points)==fit.parameters || throw(ArgumentError("fit and sweep frequency order differ"))
        photon_order isa Integer && !(photon_order isa Bool) && photon_order>0 || throw(ArgumentError("photon_order must be a positive integer"))
        isnothing(reference_frequency) || (reference_frequency isa Real && isfinite(reference_frequency) && reference_frequency>0) || throw(ArgumentError("reference frequency must be finite and positive"))
        fit.minimum_gap>0 || throw(ArgumentError("resonance gap must be resolved"))
        new{typeof(sweep),typeof(fit),typeof(metadata),typeof(diagnostics)}(sweep,fit,isnothing(reference_frequency) ? nothing : Float64(reference_frequency),Int(photon_order),metadata,diagnostics)
    end
end
function Base.getproperty(result::ResonanceResult,name::Symbol)
    name===:frequency && return getfield(result,:fit).center
    name===:minimum_gap && return getfield(result,:fit).minimum_gap
    name===:drive_time && return inv(2getfield(result,:fit).minimum_gap)
    name===:gaps && return getfield(result,:fit).gaps
    tracking=getfield(result,:sweep).tracking
    name===:tracking && return tracking
    name===:state_keys && return tracking.keys
    name===:quasienergies && return tracking.other.quasienergies
    name===:frequencies && return getfield(result,:sweep).points
    return getfield(result,name)
end
Base.propertynames(r::ResonanceResult,private::Bool=false)=(fieldnames(typeof(r))...,:frequency,:minimum_gap,:drive_time,:gaps,:tracking,:state_keys,:quasienergies,:frequencies)
