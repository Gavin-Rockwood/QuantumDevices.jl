"""
Quantum-device Hamiltonians and energy parameters use frequency units `E/h`
(cycles per unit time). For time in ns, use GHz. QuantumToolbox time evolution
requires the caller to supply the corresponding `2π` factor.
"""
module QuantumDevices
    using SciMLBase: solve
    export solve
    using LinearAlgebra
    using QuantumToolbox
    using LaTeXStrings
    import Latexify
    import Symbolics
    using Symbolics: Num

    # Shared formatting and generic utilities load before core feature types.
    include("utils/utils.jl")
    export LabelDiagnostic, TrackedData, TrackingResult, track_states
    export DressedStates

    include("symbolics/symbolics.jl")
    export component_parameters, component_operators, OperatorRef, Sym, op, val, evaluate, latex

    include("Components/Components.jl")
    export OperatorGenerator, Component, make_qubit, make_resonator
    export make_transmon, make_flux_tunable_transmon

    include("interactions/interactions.jl")
    export Interaction

    include("control_maps/control_maps.jl")
    export AbstractControlMap, AffineControlMap, TabulatedControlMap, map_values, check_control_map
    export add_control_map!, replace_control_map!, flux_map_from_frequency, control_map_diagnostics

    include("recipes/recipes.jl")
    export Recipe, Leaf, Combine, Input, ResolvedRecipe, draft

    include("device/device.jl")
    export OperatorRegistry, get_operators, get_components
    export QuantumDevice, add_parameters!, add_component!, add_interaction!, add_recipe!, replace_recipe!

    # Graph inspection needs both recipe definitions and the QuantumDevice type.
    include("recipes/graph.jl")
    export graph

    include("gates/pulses/pulses.jl")
    export AbstractPulse, AbstractPulseRecipe, PulseRecipe, PulseSequenceRecipe, Pulse, PulseSequence, pulse_duration, AbstractEnvelope, Envelope, InterpolatedEnvelope
    export EnvelopeSequence, Segment, Constant, Gaussian, GaussianZero, GaussianSquare
    export SineSquared, Cosine, Blackman, Sech, Bump, DRAG, Slepian, ErfSquare, RampedFlattop
    export pulse_tstops, envelope_tstops, parameters, calibration_parameters, getpath, setpath, haspath
    include("io/pulses.jl")
    export save_pulse, load_pulse

    include("model/model.jl")
    export NumericalModel, build, make_hamiltonian, get_hamiltonian, states
    export stateindex, state, stateenergy

    include("control_maps/frequency_conversion.jl")

    include("gates/recipes/recipes.jl")
    export ControlInput, Control, GateRecipe, PreparedGate, check_endpoints
    include("gates/evolution/evolution.jl")
    export GateEvolution, evolve, unitary, idle_phases, gate_matrix
    include("gates/calibration/calibration.jl")
    export QuantumGate, CalibrationProblem, CalibrationResult, calibrated_recipe
    export validate_calibration, gate_infidelity
    include("utils/floquet/floquet.jl")
    include("utils/spectral_tools/spectral_tools.jl")
    export FloquetBasis, FloquetSweep, get_floquet_basis, propagate_floquet_modes, propagate_floquet_states, floquet_sweep
    export AvoidedCrossingFit, ResonanceResult, fit_avoided_crossing, find_resonance, stark_shift, plot_resonance
    include("gates/repository.jl")
    export add_gate!, replace_gate!, move_gate!, get_gate, attach_artifact!
    include("io/device/device.jl")
    export save_device, load_device, register_io_codec!
end
