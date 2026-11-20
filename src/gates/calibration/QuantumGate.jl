"""A working gate recipe and at most one independently captured calibration result."""
mutable struct QuantumGate
    name::Symbol
    recipe::GateRecipe
    calibration::Union{Nothing,CalibrationResult}
    metadata::Any
    artifacts::Dict{Symbol,Any}
end

QuantumGate(name::Symbol,recipe::GateRecipe,calibration,metadata) =
    QuantumGate(name,recipe,calibration,metadata,Dict{Symbol,Any}())

# Own fields take precedence; forwarded reads always reflect the current working recipe.
function Base.getproperty(gate::QuantumGate, name::Symbol)
    name in fieldnames(QuantumGate) && return getfield(gate,name)
    return getproperty(getfield(gate,:recipe),name)
end
Base.propertynames(gate::QuantumGate, private::Bool=false) =
    Tuple(unique((fieldnames(QuantumGate)...,propertynames(getfield(gate,:recipe),private)...)))

# Working recipes remain editable, separate from any calibration snapshot.
QuantumGate(name, recipe::GateRecipe; metadata=(;)) =
    QuantumGate(Symbol(name),recipe.model_recipe isa Symbol ? deepcopy(recipe) : draft(recipe),
                nothing,deepcopy(metadata))

"""Adopt a result explicitly without validation or propagation; failed solves require opt-in."""
function QuantumGate(name, result::CalibrationResult; metadata=(;), allow_unsuccessful=false)
    SciMLBase.successful_retcode(result.solution.retcode) || allow_unsuccessful ||
        throw(ArgumentError("calibration solve was unsuccessful ($(result.solution.retcode)); use allow_unsuccessful=true to adopt its candidate explicitly"))
    return QuantumGate(Symbol(name),draft(result.recipe),result,deepcopy(metadata))
end

"""Build the current working recipe, independently of the calibrated snapshot."""
build(device::QuantumDevice, gate::QuantumGate; parameters=Dict{Num,Number}(), dressed_threads=nothing, tracking_threads=missing) =
    build(device,gate.recipe;parameters,dressed_threads,tracking_threads)

"""Evolve the current working recipe; use evolve(gate.calibration) for captured replay."""
evolve(device::QuantumDevice,gate::QuantumGate;parameters=Dict{Num,Number}(),dressed_threads=nothing,tracking_threads=missing,kwargs...) =
    evolve(build(device,gate;parameters,dressed_threads,tracking_threads);kwargs...)

"""Compile the current working recipe against a compatible captured physical model."""
make_hamiltonian(model::NumericalModel,gate::QuantumGate;kwargs...) = make_hamiltonian(model,gate.recipe;kwargs...)
check_endpoints(model::NumericalModel,gate::QuantumGate;kwargs...) = check_endpoints(model,gate.recipe;kwargs...)
