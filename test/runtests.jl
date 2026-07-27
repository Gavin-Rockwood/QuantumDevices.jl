import QuantumDevices as QD
using Test
import LinearAlgebra: norm
import QuantumToolbox as qt
import UnicodePlots as up
@testset "QuantumDevices.jl" begin
    # Write your tests here.
    startup_tests = joinpath(@__DIR__, "StartupTesting", "StartupTesting.jl")
    isfile(startup_tests) && include(startup_tests)
    include("operator_expressions_test.jl")
    include("component_hamiltonian_test.jl")
    include("model_tree_test.jl")
    include("component_model_test.jl")
    include("display_test.jl")
    include("generic_test.jl")
    include("transmon_test.jl")
    include("flux_tunable_transmon_test.jl")
    include("resonator_test.jl")
    include("gate_spec_test.jl")
    include("quantum_device_test.jl")
    include("state_tracking_test.jl")
end
