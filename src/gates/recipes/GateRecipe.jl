"""
    GateRecipe(model_recipe; controls=[], terms=[], duration=nothing, metadata=(;))

A physical gate definition with numerical pulses and optional model-reference offsets.
An omitted duration follows the latest control endpoint; an explicit duration
remains fixed, including when it initially equals that endpoint. Optimization
variables and pulse links belong to CalibrationProblem.
"""
struct GateRecipe{R,C,T,M}
    model_recipe::R
    controls::C
    terms::T
    duration::Real
    metadata::M
    duration_policy::Symbol
end

function GateRecipe(model_recipe::Union{Recipe,Symbol}; controls=Control[], terms=Interaction[],
    duration=nothing, metadata=(;))
    return _physical_gate(model_recipe,collect(controls),collect(Interaction,terms),duration,
        isnothing(duration) ? :auto : :explicit,deepcopy(metadata))
end

# Legacy positional numerical definitions have explicitly supplied durations.
GateRecipe(model_recipe, controls, terms, duration::Real, metadata) =
    GateRecipe(model_recipe; controls, terms, duration, metadata)
