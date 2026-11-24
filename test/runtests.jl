# Pkg.test() runs all coverage; core omits broader optimizer backend integrations.
length(ARGS) <= 1 || throw(ArgumentError("expected at most one test mode: core or full"))
const TEST_MODE = isempty(ARGS) ? "full" : only(ARGS)
TEST_MODE in ("core", "full") || throw(ArgumentError("test mode must be core or full (or omitted)"))

using Test
using QuantumDevices
import QuantumToolbox as QT
using LinearAlgebra
using LaTeXStrings
import Symbolics

# Shared helpers require the imports above and perform no work until called.
include("fixtures/models.jl")
include("fixtures/generator_records.jl")

@testset "QuantumDevices" begin
    # Pulse definitions and physical parameter paths.
    include("gates/pulses/envelopes.jl")
    include("gates/pulses/analytic_shapes.jl")
    include("gates/pulses/sequences.jl")
    include("gates/pulses/interpolation_flattops.jl")
    include("gates/pulses/paths.jl")
    include("gates/pulses/runtime.jl")

    include("symbolics/expressions.jl")
    include("symbolics/parameters.jl")
    include("components/components.jl")
    include("components/corrections.jl")
    include("components/transmons.jl")
    include("device/operator_registry.jl")
    include("device/interactions.jl")
    include("device/parameter_inventory.jl")
    include("device/recipe_parameter_views.jl")
    include("control_maps.jl")

    # Captured definitions and numerical construction.
    include("recipes/declarations.jl")
    include("recipes/interactions.jl")
    include("recipes/snapshots.jl")
    include("recipes/validation.jl")
    include("model/construction.jl")
    include("model/validation.jl")
    include("model/composites.jl")
    include("model/reductions.jl")
    include("model/construction_basis.jl")
    include("model/rebuild.jl")
    include("model/projection.jl")
    include("model/numerical.jl")
    include("model/scalars.jl")
    include("model/reference_vacuum.jl")

    include("gates/recipes/recipes.jl")
    include("gates/recipes/parameter_reconciliation.jl")
    include("gates/recipes/model_reuse.jl")
    include("gates/recipes/runtime.jl")
    include("gates/evolution/evolution.jl")
    include("gates/evolution/vacuum.jl")
    include("gates/calibration/calibration.jl")
    include("gates/calibration/sciml.jl")
    if TEST_MODE == "full"
        include("gates/calibration/backend_integration.jl")
    end
    include("gates/calibration/runtime.jl")
    include("gates/calibration/results.jl")
    include("gates/calibration/parameterization.jl")

    include("display/models.jl")
    include("display/gates.jl")
    include("utils/tracking/tracking.jl")
    include("utils/tracking/seeded_tracking.jl")
    include("utils/tracking/parallel.jl")
    include("utils/tracking/overlap_buffers.jl")
    include("utils/tracking/dressed_states/dressed_states.jl")
    include("utils/tracking/dressed_states/options.jl")
    include("utils/tracking/dressed_states/parallel.jl")

    include("io/pulses.jl")
    include("io/device.jl")
    include("io/generators.jl")
    include("io/parameter_migration.jl")
    include("io/recipe_input_migration.jl")
    include("io/values.jl")
    include("io/sources.jl")
    include("io/safety.jl")
    include("utils/floquet.jl")
    include("utils/spectral_tools.jl")
end
