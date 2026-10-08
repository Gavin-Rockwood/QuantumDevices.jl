"""
    QuantumDevices

Quantum-device Hamiltonians and energy parameters use frequency units, `E/h`
(cycles per unit time), not angular frequency. For time in ns, use GHz.
`numerical` preserves these units. Time evolution uses `2π * H`; `get_unitary`
applies that factor at its solver call. Direct QuantumToolbox `sesolve`/`mesolve`
calls require the caller to supply the same factor explicitly. Dimensionless
parameters (for example flux and offset charge) and phases in radians retain
their own conventions.
"""
module QuantumDevices
    using LinearAlgebra
    using SparseArrays
    import JSON3
    import JLD2
    import SciMLBase
    import CairoMakie
    import LsqFit
    import UnicodePlots
    using SpecialFunctions: erf
    using SciMLBase: MatrixOperator, solve
    using QuantumToolbox

    # Base Overloads
    import Base.truncate

    include("utils/utils.jl")
    export getpath, haspath, setpath, setpath!, track_states, get_dressed_states
    export state_amplitudes, plot_trajectories, plot_trajectories!
    export FloquetBasis, get_floquet_basis, propagate_floquet_modes, floquet_sweep
    export AvoidedCrossingFit, ResonanceResult, fit_avoided_crossing, find_resonance, plot_resonance
    include("symbolics/symbolics.jl")
    export Sym, op, param, val, call, numerical

    include("Components/Components.jl")
    export Component, make_qubit, make_resonator, make_transmon, make_tunable_transmon

    include("model/model.jl")
    export DeviceModel, make_model
    include("gates/gates.jl")
    export AbstractPulse, Pulse, DeviceGate, pulse_tstops
    export AbstractEnvelope, Constant, Gaussian, SineSquared, RampedFlattop, Envelope, envelope_value, validate_envelope
    export GaussianZero, GaussianSquare, Sech, Cosine, Blackman, Bump, ErfSquare, Slepian, DRAG, envelope_tstops
    export AbstractCarrier, SineCarrier, IQCarrier, Carrier, carrier_value
    export parameters, calibration_values, calibration_problem, calibrated_gate, calibrate
    export CalibrationProblem, solve
    export AbstractCalibrationSetup, SciMLCalibrationSetup, get_unitary, gate_unitary, unitary_fidelity, unitary_infidelity, gate_infidelity

    include("time_evolution/time_evolution.jl")
    export get_gate_matrix
    include("io/io.jl")
    export save, load




end
