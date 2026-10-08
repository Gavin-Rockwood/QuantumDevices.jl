using Test
using LinearAlgebra
using QuantumToolbox
using QuantumDevices
import QuantumDevices as QD

@testset "QuantumDevices V5" begin
    include("paths.jl")
    include("symbolics.jl")
    include("operators.jl")
    include("components.jl")
    include("tracking.jl")
    include("floquet.jl")
    include("spectral_tools.jl")
    include("models.jl")
    include("pulses.jl")
    include("gates.jl")
    include("fidelities.jl")
    include("time_evolution.jl")
    include("trajectories.jl")
    include("gate_matrices.jl")
    include("envelopes.jl")
    include("bump.jl")
    include("io_bundle.jl")
    include("calibration.jl")
end
