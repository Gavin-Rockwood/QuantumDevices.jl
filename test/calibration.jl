import SciMLBase
import OptimizationOptimJL

struct TestCalibrationEnvelope <: AbstractEnvelope
    value::Float64
end

QD.envelope_value(pulse::TestCalibrationEnvelope, t, duration) = pulse.value
QD.parameters(pulse::TestCalibrationEnvelope) =
    Dict{String,Tuple}("value" => ("value", pulse.value))

@testset "Gate calibration" begin
    q = make_qubit("calibration", 0.0)
    model = make_model([q], 0 * op(:calibration_x), (;))
    gate = DeviceGate(
        (drive = Pulse(Constant(); amplitude=0.1, duration=1.0),),
        param(:drive) * op(:calibration_x),
        1.0,
    )
    target = -im * sigmax()
    amplitude(gate) = getpath(gate, "drive/amplitude")
    offset = "parameters/drive/offset"
    start = "parameters/drive/delay"
    values = parameters(gate)
    @test values["drive/amplitude"] == ("drive/amplitude", 0.1)
    @test values["drive/offset"] == ("drive/offset", 0)
    @test values["drive/delay"] == ("drive/delay", 0)
    @test getpath(gate, values["drive/amplitude"][1]) == 0.1
    @test !haskey(values, "drive")
    setup = CalibrationProblem(model, gate, ["drive/amplitude"], target)

    @test setup isa CalibrationProblem
    @test setup isa SciMLCalibrationSetup
    @test solve === SciMLBase.solve
    @test calibration_problem(model, gate, ["drive/amplitude"], target) isa CalibrationProblem
    @test setup.problem isa SciMLBase.OptimizationProblem
    @test setup.problem.u0 == [0.1]

    solution = solve(
        setup,
        OptimizationOptimJL.NelderMead();
        maxiters = 80,
    )
    calibrated = calibrated_gate(setup, solution.u)
    @test solution isa SciMLBase.AbstractOptimizationSolution
    @test solution.objective < 1e-7
    @test calibrated.parameters.drive.amplitude ≈ 0.25 atol = 1e-3
    @test gate.parameters.drive.amplitude == 0.1

    dense_setup = CalibrationProblem(model, gate, ["drive/amplitude"], target;
        dense=true, lb=[0.0], ub=[1.0])
    @test dense_setup.problem.lb == [0.0]
    @test dense_setup.problem.ub == [1.0]
    # Use an unbounded problem with NelderMead, retaining the native solve interface.
    dense_setup = CalibrationProblem(model, gate, ["drive/amplitude"], target; dense=true)
    dense_solution = solve(dense_setup, OptimizationOptimJL.NelderMead(); maxiters=80)
    @test dense_solution isa SciMLBase.AbstractOptimizationSolution
    @test dense_solution.objective < 1e-7
    @test calibrated_gate(dense_setup, dense_solution.u).parameters.drive.amplitude ≈ 0.25 atol=1e-3

    direct_solution = SciMLBase.solve(
        setup.problem,
        OptimizationOptimJL.NelderMead();
        maxiters = 80,
    )
    direct_gate = calibrated_gate(setup, direct_solution.u)
    @test direct_solution.objective < 1e-7
    @test direct_gate.parameters.drive.amplitude ≈ 0.25 atol = 1e-3
    @test typeof(solution) == typeof(direct_solution)

    scalar_gate = DeviceGate((drive = 0.1,), param(:drive) * op(:calibration_x), 1.0)
    scalar_drive = "parameters/drive"
    scalar_objective = (_, candidate, target_value) ->
        (parameters(candidate)["drive"][2] - target_value)^2
    @test parameters(scalar_gate) == Dict("drive" => ("drive", 0.1))
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
    @test_throws ArgumentError CalibrationProblem(model, scalar_gate, ["drive"], 1.25;
        objective=scalar_objective, dense=true)
    @test scalar_setup.problem.f([1.25], SciMLBase.NullParameters()) == 0

    # Solver options and callbacks must behave exactly as on the wrapped problem.
    passthrough = CalibrationProblem(model, scalar_gate, ["drive"], 1.25;
        objective=scalar_objective)
    calls = Ref(0)
    callback = (args...) -> (calls[] += 1; false)
    limited = solve(passthrough, OptimizationOptimJL.NelderMead(); maxiters=4, callback)
    direct_limited = SciMLBase.solve(passthrough.problem,
        OptimizationOptimJL.NelderMead(); maxiters=4)
    @test calls[] > 0
    @test limited.u == direct_limited.u
    @test limited.objective == direct_limited.objective
    @test limited.retcode == direct_limited.retcode

    custom_gate = DeviceGate(
        (drive = Pulse(TestCalibrationEnvelope(0.2); duration=1.0),),
        param(:drive) * op(:calibration_x),
        1.0,
    )
    custom_value = "parameters/drive/envelope/value"
    custom_objective = (_, candidate, target_value) ->
        (parameters(candidate)["drive/envelope/value"][2] - target_value)^2
    @test parameters(custom_gate)["drive/envelope/value"] == ("drive/envelope/value", 0.2)
    custom_setup = calibration_problem(
        model,
        custom_gate,
        ["drive/envelope/value"],
        0.8;
        objective = custom_objective,
    )
    rebuilt_custom = calibrated_gate(custom_setup, [0.8])
    @test rebuilt_custom.parameters.drive.envelope isa TestCalibrationEnvelope
    @test rebuilt_custom.parameters.drive.envelope.value == 0.8
    @test custom_gate.parameters.drive.envelope.value == 0.2

    changed_type = setpath(gate, "drive/amplitude", 1)
    @test amplitude(changed_type) === 1
    @test amplitude(gate) === 0.1
    @test_throws ArgumentError setpath(gate, "duration", -1.0)
    @test_throws ArgumentError setpath(
        gate,
        "parameters/drive/duration",
        2.0,
    )
    @test_throws ArgumentError setpath(gate.parameters.drive, "carrier", :unknown)

    generic_gate = DeviceGate(
        (drive = Pulse(Envelope((p, t, duration) -> p.scale, (; scale=1.0)); amplitude=0.1, duration=1.0),),
        param(:drive) * op(:calibration_x), 1.0)
    @test parameters(generic_gate)["drive/amplitude"] == ("drive/amplitude", 0.1)
    generic_setup = calibration_problem(model, generic_gate, ["drive/amplitude"], target)
    generic_solution, generic_calibrated = calibrate(generic_setup, OptimizationOptimJL.NelderMead(); maxiters = 80)
    @test generic_solution.objective < 1e-7
    @test amplitude(generic_calibrated) ≈ 0.25 atol = 1e-3
    @test amplitude(generic_gate) == 0.1
    @test generic_calibrated.parameters.drive.envelope.callable === generic_gate.parameters.drive.envelope.callable


    # Nested calibration paths rebuild shapes/carriers and retain inferred timing.
    nested_drive = Pulse(RampedFlattop(0.2; ramp_up=Gaussian(0.2));
        duration=1.0, delay=0.1, amplitude=0.1, carrier=SineCarrier(0.5))
    nested_gate = DeviceGate((drive=nested_drive,), param(:drive)*op(:calibration_x))
    names = ["drive/amplitude", "drive/carrier/frequency", "drive/envelope/ramp_up/sigma", "drive/duration", "drive/delay"]
    wanted = [0.2, 0.7, 0.15, 1.2, 0.3]
    loss = (_,g,v) -> sum((getpath(g,k)-x)^2 for (k,x) in zip(names,v))
    nested_setup = CalibrationProblem(model, nested_gate, names, wanted; objective=loss)
    @test nested_setup.problem.u0 == [0.1, 0.5, 0.2, 1.0, 0.1]
    rebuilt_nested = calibrated_gate(nested_setup, wanted)
    @test rebuilt_nested.duration ≈ 1.5
    @test rebuilt_nested.duration_override === nothing
    @test nested_gate.duration ≈ 1.1
    @test nested_setup.problem.f(wanted, SciMLBase.NullParameters()) == 0
    nested_solution = solve(nested_setup, OptimizationOptimJL.NelderMead(); maxiters=1000, abstol=1e-10)
    @test nested_solution.objective < 1e-6
    @test loss(model, calibrated_gate(nested_setup, nested_solution.u), wanted) < 1e-6

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
    @test_throws ArgumentError calibration_problem(model, gate, ["drive/carrier"], target)
    @test_throws ArgumentError calibration_problem(model, gate, "drive/amplitude", target)
    @test_throws ArgumentError calibration_problem(model, gate, [:amplitude], target)
    @test parameters(gate)["drive/carrier"] == ("drive/carrier", nothing)
    named_gate = DeviceGate((q1_nu = 0.2, drive = Pulse(Constant(); amplitude=0.1, duration=1.0)), op(:calibration_x), 1.0)
    named_setup = calibration_problem(model, named_gate, ["q1_nu", "drive/amplitude"], target)
    @test named_setup.problem.u0 == [0.2, 0.1]
    named_updated = calibrated_gate(named_setup, [0.7, 0.8])
    @test named_updated.parameters.q1_nu == 0.7
    @test named_updated.parameters.drive.amplitude == 0.8
    @test named_gate.parameters.q1_nu == 0.2
    collision_gate = DeviceGate((duration = 0.2,), 0, 1.0)
    @test parameters(collision_gate)["duration"] == ("parameters/duration", 0.2)
    collision_setup = calibration_problem(model, collision_gate, ["duration"], target)
    collision_updated = calibrated_gate(collision_setup, [0.3])
    @test collision_updated.duration == 1.0
    @test collision_updated.parameters.duration == 0.3

end
