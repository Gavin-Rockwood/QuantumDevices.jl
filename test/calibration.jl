import SciMLBase
import OptimizationOptimJL

struct TestCalibrationPulse <: AbstractPulse
    value::Float64
end

QD.pulse_value(pulse::TestCalibrationPulse, t, duration) = pulse.value
QD.parameters(pulse::TestCalibrationPulse) =
    Dict{String,Tuple}("value" => ("value", pulse.value))

@testset "Gate calibration" begin
    q = make_qubit("calibration", 0.0)
    model = make_model([q], 0 * op(:calibration_x), (;))
    gate = DeviceGate(
        (drive = constant_pulse(0.4),),
        param(:drive) * op(:calibration_x),
        1.0,
    )
    target = -im * sigmax()
    amplitude(gate) = getpath(gate, "drive/amplitude")
    offset = "parameters/drive/parameters/offset"
    start = "parameters/drive/parameters/start"
    values = parameters(gate)
    @test values["drive/amplitude"] == ("drive/amplitude", 0.4)
    @test values["drive/offset"] == ("drive/offset", 0)
    @test values["drive/start"] == ("drive/start", 0)
    @test getpath(gate, values["drive/amplitude"][1]) == 0.4
    @test !haskey(values, "drive")
    setup = calibration_problem(model, gate, ["drive/amplitude"], target)

    @test setup isa SciMLCalibrationSetup
    @test setup.problem isa SciMLBase.OptimizationProblem
    @test setup.problem.u0 == [0.4]

    solution, calibrated = calibrate(
        setup,
        OptimizationOptimJL.NelderMead();
        maxiters = 80,
    )
    @test solution.objective < 1e-7
    @test calibrated.parameters.drive.parameters.amplitude ≈ pi / 2 atol = 1e-3
    @test gate.parameters.drive.parameters.amplitude == 0.4

    direct_solution = SciMLBase.solve(
        setup.problem,
        OptimizationOptimJL.NelderMead();
        maxiters = 80,
    )
    direct_gate = calibrated_gate(setup, direct_solution.u)
    @test direct_solution.objective < 1e-7
    @test direct_gate.parameters.drive.parameters.amplitude ≈ pi / 2 atol = 1e-3

    scalar_gate = DeviceGate((drive = 0.4,), param(:drive) * op(:calibration_x), 1.0)
    scalar_drive = "parameters/drive"
    scalar_objective = (_, candidate, target_value) ->
        (parameters(candidate)["drive"][2] - target_value)^2
    @test parameters(scalar_gate) == Dict("drive" => ("drive", 0.4))
    scalar_setup = calibration_problem(
        model,
        scalar_gate,
        ["drive"],
        1.25;
        objective = scalar_objective,
        lb = [0.0],
        ub = [2.0],
    )
    @test scalar_setup.problem.lb == [0.0]
    @test scalar_setup.problem.ub == [2.0]
    @test scalar_setup.problem.f([1.25], SciMLBase.NullParameters()) == 0

    custom_gate = DeviceGate(
        (drive = TestCalibrationPulse(0.2),),
        param(:drive) * op(:calibration_x),
        1.0,
    )
    custom_value = "parameters/drive/value"
    custom_objective = (_, candidate, target_value) ->
        (parameters(candidate)["drive/value"][2] - target_value)^2
    @test parameters(custom_gate)["drive/value"] == ("drive/value", 0.2)
    custom_setup = calibration_problem(
        model,
        custom_gate,
        ["drive/value"],
        0.8;
        objective = custom_objective,
    )
    rebuilt_custom = calibrated_gate(custom_setup, [0.8])
    @test rebuilt_custom.parameters.drive isa TestCalibrationPulse
    @test rebuilt_custom.parameters.drive.value == 0.8
    @test custom_gate.parameters.drive.value == 0.2

    changed_type = setpath(gate, "drive/amplitude", 1)
    @test amplitude(changed_type) === 1
    @test amplitude(gate) === 0.4
    @test_throws ArgumentError setpath(gate, "duration", -1.0)
    @test_throws ArgumentError setpath(
        gate,
        "parameters/drive/parameters/stop",
        2.0,
    )
    @test_throws ArgumentError setpath(gate.parameters.drive, "name", :unknown)

    generic_gate = DeviceGate(
        (drive = GenericPulseFunction((p, t) -> p.amplitude, (amplitude = 0.4,)),),
        param(:drive) * op(:calibration_x), 1.0)
    @test parameters(generic_gate)["drive/amplitude"] == ("drive/amplitude", 0.4)
    generic_setup = calibration_problem(model, generic_gate, ["drive/amplitude"], target)
    generic_solution, generic_calibrated = calibrate(generic_setup, OptimizationOptimJL.NelderMead(); maxiters = 80)
    @test generic_solution.objective < 1e-7
    @test amplitude(generic_calibrated) ≈ pi / 2 atol = 1e-3
    @test amplitude(generic_gate) == 0.4
    @test pulse_function(generic_calibrated.parameters.drive) === pulse_function(generic_gate.parameters.drive)


    @test unitary_infidelity(target, target) ≈ 0 atol = 1e-14
    @test_throws ArgumentError calibration_problem(model, gate, [], target)
    @test_throws ArgumentError calibration_problem(model, gate, ["missing"], target)
    @test_throws ArgumentError calibration_problem(
        model,
        gate,
        ["drive"],
        target,
    )
    @test_throws ArgumentError calibration_problem(model, gate, ["drive/amplitude", "drive/amplitude"], target)
    @test_throws DimensionMismatch calibrated_gate(setup, [0.1, 0.2])
    @test_throws ArgumentError calibration_problem(model, gate, ["drive/stop"], target)
    @test_throws ArgumentError calibration_problem(model, gate, "drive/amplitude", target)
    @test_throws ArgumentError calibration_problem(model, gate, [:amplitude], target)
    @test parameters(gate)["drive/stop"] == ("drive/stop", nothing)
    named_gate = DeviceGate((q1_nu = 0.2, drive = constant_pulse(0.4)), op(:calibration_x), 1.0)
    named_setup = calibration_problem(model, named_gate, ["q1_nu", "drive/amplitude"], target)
    @test named_setup.problem.u0 == [0.2, 0.4]
    named_updated = calibrated_gate(named_setup, [0.7, 0.8])
    @test named_updated.parameters.q1_nu == 0.7
    @test named_updated.parameters.drive.parameters.amplitude == 0.8
    @test named_gate.parameters.q1_nu == 0.2
    collision_gate = DeviceGate((duration = 0.2,), 0, 1.0)
    @test parameters(collision_gate)["duration"] == ("parameters/duration", 0.2)
    collision_setup = calibration_problem(model, collision_gate, ["duration"], target)
    collision_updated = calibrated_gate(collision_setup, [0.3])
    @test collision_updated.duration == 1.0
    @test collision_updated.parameters.duration == 0.3

end
