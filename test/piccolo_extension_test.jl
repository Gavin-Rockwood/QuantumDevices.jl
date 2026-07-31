using Test
import Piccolo
import QuantumDevices as QD

@testset "Piccolo gate optimization extension" begin
    qubit = QD.Component(QD.QubitSpec(0.0; fixed_frequency = true), :q)
    drive =
        QD.DeviceParameter(:drive; domain = QD.realdomain(), fixed = false, default = 0.0)
    spec = QD.ModelSpec(
        :piccolo_qubit,
        [qubit];
        interactions = (QD.param(drive) * QD.op(:q, :x),),
        parameters = (drive,),
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    path = QD.ParamPath(:drive)
    initial = zeros(1, 8)
    initial[1, 2:7] .= 0.4
    objective = QD.UnitaryObjective(ComplexF64[0 1; 1 0]; states = [(0,), (1,)])
    problem = QD.PiccoloOptimizationProblem(
        :piccolo_x,
        spec,
        Dict(path => (-1.0, 1.0)),
        objective;
        duration = 1.0,
        knots = 8,
        initial_controls = initial,
        template_kwargs = Dict(:ddu_bound => 2.0),
        solve_kwargs = Dict(:max_iter => 2, :verbose => false, :print_level => 0),
    )
    result = QD.optimize_gate(problem)
    @test result.gate isa QD.GateSpec
    @test result.settings[:backend] == :Piccolo
    @test isfinite(result.fidelity)
    @test result.gate.controls[path](0.0) ≈ 0.0 atol = 1e-10
    @test result.gate.controls[path](1.0) ≈ 0.0 atol = 1e-10

    verified = QD._evaluate_gate_objective(result.gate, objective)
    @test result.fidelity ≈ verified.fidelity atol = 1e-10
    @test result.leakage ≈ verified.leakage atol = 1e-10
    @test result.backend_result isa
          Piccolo.Control.QuantumControlProblems.QuantumControlProblem
end

@testset "streamlined Piccolo optimization" begin
    qubit = QD.Component(QD.QubitSpec(0.0; fixed_frequency = true), :q)
    drive =
        QD.DeviceParameter(:drive; domain = QD.realdomain(), fixed = false, default = 0.0)
    detuning = QD.DeviceParameter(
        :detuning;
        domain = QD.realdomain(),
        fixed = false,
        default = 0.0,
    )
    spec = QD.ModelSpec(
        :streamlined_piccolo,
        [qubit];
        interactions = (
            QD.param(drive) * QD.op(:q, :x),
            QD.param(detuning) * QD.op(:q, :z),
        ),
        parameters = (drive, detuning),
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    drive_path = QD.ParamPath(:drive)
    detuning_path = QD.ParamPath(:detuning)
    gate = QD.GateSpec(
        :streamlined_piccolo_x,
        spec;
        duration = 1.0,
        controls = Dict(drive_path => (time -> 0.4 * sinpi(time)^2)),
    )
    target = ComplexF64[0 1; 1 0]
    result = QD.optimize(
        gate,
        target,
        drive_path => (-1.0, 1.0);
        backend = :piccolo,
        knots = 8,
        solve_kwargs = Dict(:max_iter => 1, :verbose => false, :print_level => 0),
        template_kwargs = Dict(:ddu_bound => 2.0),
    )
    @test result.gate isa QD.GateSpec
    @test result.gate.control_recipes === nothing
    @test result.gate.controls[drive_path](0.0) == 0.0
    @test result.gate.controls[drive_path](1.0) == 0.0
    @test result.settings[:backend] == :Piccolo
    @test isfinite(result.fidelity)

    minimum_time = QD.optimize(
        gate,
        target,
        (drive_path => (-1.0, 1.0), :duration => (0.8, 1.2));
        backend = :piccolo,
        knots = 8,
        final_fidelity = 0.1,
        solve_kwargs = Dict(:max_iter => 1, :verbose => false, :print_level => 0),
        template_kwargs = Dict(:ddu_bound => 2.0),
    )
    @test minimum_time.settings[:template] == :MinimumTimeProblem
    @test 0.8 <= minimum_time.gate.duration <= 1.2
    @test minimum_time.minimizer[:duration] == minimum_time.gate.duration
    @test minimum_time.gate.controls[drive_path](0.0) == 0.0
    @test minimum_time.gate.controls[drive_path](minimum_time.gate.duration) == 0.0

    two_controls = QD.GateSpec(
        :two_controls,
        spec;
        duration = 1.0,
        controls = Dict(
            drive_path => (time -> 0.4 * sinpi(time)^2),
            detuning_path => (time -> 0.1 * sinpi(time)^2),
        ),
    )
    @test_throws Exception QD.optimize(
        two_controls,
        target,
        drive_path => (-1.0, 1.0);
        backend = :piccolo,
        knots = 5,
    )

    recipe_gate = QD.GateSpec(
        :recipe_piccolo,
        spec;
        duration = 1.0,
        parameters = (amplitude = 0.4,),
        controls = Dict(
            drive_path => ((p, time) -> p.amplitude * sinpi(time / p.duration)^2),
        ),
    )
    @test_throws Exception QD.optimize(
        recipe_gate,
        target,
        :amplitude => (0.0, 1.0);
        backend = :piccolo,
        knots = 5,
    )
end

@testset "Piccolo carrier-envelope optimization" begin
    qubit = QD.Component(QD.QubitSpec(0.2; fixed_frequency = true), :q)
    drive_i = QD.DeviceParameter(
        :drive_i;
        domain = QD.realdomain(),
        fixed = false,
        default = 0.0,
    )
    drive_q = QD.DeviceParameter(
        :drive_q;
        domain = QD.realdomain(),
        fixed = false,
        default = 0.0,
    )
    spec = QD.ModelSpec(
        :piccolo_carrier,
        [qubit];
        interactions = (
            QD.param(drive_i) * QD.op(:q, :x),
            QD.param(drive_q) * QD.op(:q, :x),
        ),
        parameters = (drive_i, drive_q),
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    i_path = QD.ParamPath(:drive_i)
    q_path = QD.ParamPath(:drive_q)
    gate = QD.GateSpec(
        :piccolo_carrier_x,
        spec;
        duration = 1.0,
        parameters = NamedTuple(),
        controls = Dict(
            i_path => QD.CarrierControl(
                time -> 0.4 * sinpi(time)^2;
                frequency = 0.2,
            ),
            q_path => QD.CarrierControl(
                time -> 0.0;
                frequency = 0.2,
                phase = -π / 2,
            ),
        ),
    )
    objective = QD.UnitaryObjective(
        ComplexF64[0 1; 1 0];
        states = [(0,), (1,)],
        frame_frequencies = (q = 0.2,),
    )
    result = QD.optimize(
        gate,
        objective,
        (i_path => (-1.0, 1.0), q_path => (-1.0, 1.0));
        backend = :piccolo,
        knots = 8,
        solve_kwargs = Dict(:max_iter => 1, :verbose => false, :print_level => 0),
        template_kwargs = Dict(:ddu_bound => 2.0),
    )
    @test Set(result.settings[:carrier_paths]) == Set((i_path, q_path))
    @test result.settings[:initial_controls][1, :] !=
          [gate.controls[sort([i_path, q_path]; by = path -> path.parts)[1]](
              time,
          ) for time in range(0.0, 1.0; length = 8)]
    @test result.gate.control_recipes === nothing
    @test result.gate.controls[i_path](0.0) == 0.0
    @test result.gate.controls[q_path](1.0) == 0.0
    verified = QD._evaluate_gate_objective(result.gate, objective)
    @test result.fidelity ≈ verified.fidelity atol = 1e-10

    @test_throws Exception QD.optimize(
        gate,
        objective,
        (i_path => (-1.0, 1.0), :duration => (0.8, 1.2));
        backend = :piccolo,
        knots = 5,
    )

    nonlinear_drive = QD.DeviceParameter(
        :nonlinear_drive;
        domain = QD.realdomain(),
        fixed = false,
        default = 0.0,
    )
    nonlinear_spec = QD.ModelSpec(
        :nonlinear_carrier,
        [qubit];
        interactions = (QD.param(nonlinear_drive)^2 * QD.op(:q, :x),),
        parameters = (nonlinear_drive,),
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    nonlinear_path = QD.ParamPath(:nonlinear_drive)
    nonlinear_gate = QD.GateSpec(
        :nonlinear_carrier,
        nonlinear_spec;
        duration = 1.0,
        parameters = NamedTuple(),
        controls = Dict(
            nonlinear_path => QD.CarrierControl(
                time -> 0.2 * sinpi(time)^2;
                frequency = 0.2,
            ),
        ),
    )
    @test_throws Exception QD.optimize(
        nonlinear_gate,
        ComplexF64[0 1; 1 0],
        nonlinear_path => (-1.0, 1.0);
        backend = :piccolo,
        states = [(0,), (1,)],
        knots = 5,
    )
end

@testset "Piccolo multi-control pulse conversion" begin
    qubit = QD.Component(QD.QubitSpec(0.0; fixed_frequency = true), :q)
    drive =
        QD.DeviceParameter(:drive; domain = QD.realdomain(), fixed = false, default = 0.0)
    detuning = QD.DeviceParameter(
        :detuning;
        domain = QD.realdomain(),
        fixed = false,
        default = 0.1,
    )
    spec = QD.ModelSpec(
        :piccolo_multicontrol,
        [qubit];
        interactions = (
            QD.param(drive) * QD.op(:q, :x),
            QD.param(detuning) * QD.op(:q, :z),
        ),
        parameters = (drive, detuning),
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    drive_path = QD.ParamPath(:drive)
    detuning_path = QD.ParamPath(:detuning)
    initial = zeros(2, 5)
    initial[1, :] .= 0.1
    initial[2, 2:4] .= 0.2
    problem = QD.PiccoloOptimizationProblem(
        :piccolo_multicontrol,
        spec,
        Dict(drive_path => (-1.0, 1.0), detuning_path => (-0.5, 0.5)),
        QD.UnitaryObjective(ComplexF64[0 1; 1 0]; states = [(0,), (1,)]);
        duration = 1.0,
        knots = 5,
        initial_controls = initial,
        template_kwargs = Dict(:ddu_bound => 2.0),
        solve_kwargs = Dict(:max_iter => 1, :verbose => false, :print_level => 0),
    )
    result = QD.optimize_gate(problem)
    @test result.gate isa QD.GateSpec
    @test result.gate.controls[drive_path](0.0) == 0.0
    @test result.gate.controls[drive_path](1.0) == 0.0
    @test result.gate.controls[detuning_path](0.0) == 0.1
    @test result.gate.controls[detuning_path](1.0) == 0.1
end
