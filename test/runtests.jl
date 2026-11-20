using Test
using QuantumDevices
import QuantumToolbox as QT
using LinearAlgebra
using LaTeXStrings
import Symbolics

# Reconstruct a transient reference map from the actual stored eigenstates.
function input_basis_map(model)
    dimension = length(model.states)
    return reduce(+, (ket * QT.basis(dimension, index - 1)' for (index, ket) in enumerate(model.states)))
end

include("v7/model_fixtures.jl")

@testset "V7 package" begin
    include("v7/pulses.jl")
    include("v7/symbolics.jl")
    include("v7/symbolic_parameters.jl")
    include("v7/interaction_defaults.jl")
    include("v7/components.jl")
    include("v7/operator_registry.jl")
    include("v7/transmons.jl")
    include("v7/models.jl")
    include("v7/validation.jl")
    include("v7/validation_boundaries.jl")
    include("v7/parameter_reconciliation.jl")
    include("v7/composites.jl")
    include("v7/reductions.jl")
    include("v7/construction_basis.jl")
    include("v7/model_rebuild.jl")
    include("v7/projection_cleanup.jl")
    include("v7/recipe_interactions.jl")
    include("v7/recipes_rebuild.jl")
    include("v7/numerical.jl")
    include("v7/numerical_scalars.jl")
    include("v7/reference_vacuum.jl")
    include("v7/vacuum_evolution.jl")
    include("v7/gate_recipes.jl")
    include("v7/model_gate_reuse.jl")
    include("v7/gate_evolution.jl")
    include("v7/calibration.jl")
    include("v7/calibration_sciml.jl")
    include("v7/calibration_results.jl")
    include("v7/calibration_parameterization.jl")
    include("v7/displays.jl")
    include("v7/gate_displays.jl")
    include("v7/tracking.jl")
    include("v7/seeded_tracking.jl")
    include("v7/snapshots.jl")
    include("v7/device_io.jl")
    include("v7/device_io_values.jl")
    include("v7/device_io_sources.jl")
    include("v7/device_io_safety.jl")
    include("v7/floquet.jl")
    include("v7/spectral_tools.jl")
end
