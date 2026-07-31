using Test
import LinearAlgebra: I, norm
import QuantumDevices as QD
import Optim

function driven_qubit_fixture()
    qubit = QD.Component(QD.QubitSpec(0.0; fixed_frequency = true), :q)
    drive =
        QD.DeviceParameter(:drive; domain = QD.realdomain(), fixed = false, default = 0.0)
    spec = QD.ModelSpec(
        :driven_qubit,
        [qubit];
        interactions = (QD.param(drive) * QD.op(:q, :x),),
        parameters = (drive,),
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    (; spec, path = QD.ParamPath(:drive))
end

@testset "streamlined gate optimization" begin
    (; spec, path) = driven_qubit_fixture()
    gate = QD.GateSpec(
        :concise_x,
        spec;
        duration = 1.0,
        parameters = (amplitude = 0.2,),
        controls = Dict(path => ((p, time) -> p.amplitude * sin(π * time / p.duration)^2)),
    )
    target = ComplexF64[0 1; 1 0]
    inline = QD.optimize(gate, target, :amplitude => (0.0, 1.0); iterations = 40)
    @test inline.fidelity > 0.999
    @test inline.minimizer[:amplitude] ≈ 0.5 atol = 2e-3
    @test inline.gate.parameters.amplitude ≈ inline.minimizer[:amplitude]
    @test gate.parameters.amplitude == 0.2
    @test occursin("gate=concise_x", sprint(show, inline))
    @test occursin("amplitude=", sprint(show, inline))

    separate = QD.optimize(
        gate,
        target,
        :amplitude;
        bounds = (amplitude = (0.0, 1.0),),
        options = Optim.Options(iterations = 40),
    )
    @test separate.fidelity > 0.999
    @test separate.minimizer == inline.minimizer

    state_result = QD.optimize(
        gate,
        QD.StateTransferObjective((0,), (1,); states = [(0,), (1,)]),
        :amplitude => (0.0, 1.0);
        iterations = 40,
    )
    @test state_result.fidelity > 0.999

    duration_result = QD.optimize(
        gate,
        target,
        (:amplitude => (0.0, 1.0), :duration => (0.75, 1.25));
        states = [(0,), (1,)],
        iterations = 5,
    )
    @test 0.75 <= duration_result.gate.duration <= 1.25
    @test 0.0 <= duration_result.gate.parameters.amplitude <= 1.0

    @test_throws Exception QD.optimize(
        gate,
        target,
        :amplitude => (0.0, 1.0);
        options = Optim.Options(iterations = 2),
        iterations = 2,
    )
    @test_throws Exception QD.optimize(gate, target, :missing => (0.0, 1.0))
    @test_throws Exception QD.optimize(
        gate,
        target,
        path => (-1.0, 1.0);
        backend = :piccolo,
    )

    q1 = QD.Component(QD.QubitSpec(0.0; fixed_frequency = true), :q1)
    q2 = QD.Component(QD.QubitSpec(0.0; fixed_frequency = true), :q2)
    pair_spec = QD.ModelSpec(
        :ambiguous_pair,
        [q1, q2];
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    pair_gate = QD.GateSpec(:pair, pair_spec; duration = 1.0, parameters = (scale = 0.5,))
    @test_throws Exception QD.optimize(
        pair_gate,
        Matrix{ComplexF64}(I, 4, 4),
        :scale => (0.0, 1.0),
    )
end

@testset "virtual-Z frame correction" begin
    qubit = QD.Component(QD.QubitSpec(0.37; fixed_frequency = true), :q)
    spec = QD.ModelSpec(
        :frame_qubit,
        [qubit];
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    gate = QD.GateSpec(
        :frame_idle,
        spec;
        duration = 0.83,
        parameters = (unused = 0.0,),
    )
    states = [(0,), (1,)]
    target = Matrix{ComplexF64}(I, 2, 2)
    raw = QD._evaluate_gate_objective(
        gate,
        QD.UnitaryObjective(target; states),
    )
    corrected = QD._evaluate_gate_objective(
        gate,
        QD.UnitaryObjective(
            target;
            states,
            frame_frequencies = (q = 0.37,),
        ),
    )
    @test raw.fidelity < 0.5
    @test corrected.fidelity ≈ 1.0 atol = 1e-8
    @test corrected.leakage < 1e-8

    duration_result = QD.optimize(
        gate,
        target,
        :duration => (0.5, 1.2);
        states,
        frame_frequencies = (q = 0.37,),
        iterations = 5,
    )
    @test duration_result.fidelity ≈ 1.0 atol = 1e-8
    @test 0.5 <= duration_result.gate.duration <= 1.2

    bad_objective = QD.UnitaryObjective(
        target;
        states,
        frame_frequencies = (wrong = 0.37,),
    )
    @test_throws Exception QD._evaluate_gate_objective(gate, bad_objective)
    @test_throws Exception QD.UnitaryObjective(
        target;
        states,
        frame_frequencies = (q = Inf,),
    )

    q1 = QD.Component(QD.QubitSpec(0.31; fixed_frequency = true), :q1)
    q2 = QD.Component(QD.QubitSpec(0.47; fixed_frequency = true), :q2)
    pair_spec = QD.ModelSpec(
        :frame_pair,
        [q1, q2];
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    pair_gate = QD.GateSpec(:frame_pair_idle, pair_spec; duration = 0.73)
    pair_states = [(0, 0), (0, 1), (1, 0), (1, 1)]
    pair_metrics = QD._evaluate_gate_objective(
        pair_gate,
        QD.UnitaryObjective(
            Matrix{ComplexF64}(I, 4, 4);
            states = pair_states,
            frame_frequencies = (q1 = 0.31, q2 = 0.47),
        ),
    )
    @test pair_metrics.fidelity ≈ 1.0 atol = 1e-8
end

@testset "gate optimization inputs" begin
    @test_throws Exception QD.GateVariable(:x; initial = 0.0, lower = 1.0, upper = 0.0)
    @test_throws Exception QD.GateVariable(:x; initial = 2.0, lower = 0.0, upper = 1.0)

    variable = QD.GateVariable(:amplitude; initial = 0.2, lower = 0.0, upper = 1.0)
    objective = QD.UnitaryObjective(ComplexF64[0 1; 1 0]; states = [(0,), (1,)])
    @test_throws Exception QD.EnvelopeOptimizationProblem(
        identity,
        [variable, variable],
        objective,
    )
    @test_throws Exception QD.UnitaryObjective(
        Matrix{ComplexF64}(I, 3, 3);
        states = [(0,), (1,)],
    )
end

@testset "bounded envelope optimization" begin
    (; spec, path) = driven_qubit_fixture()
    variable = QD.GateVariable(:amplitude; initial = 0.2, lower = 0.0, upper = 1.0)
    builder =
        values -> QD.GateSpec(
            :x_gate,
            spec;
            duration = 1.0,
            controls = Dict(path => (time -> values.amplitude * sin(π * time)^2)),
        )
    unitary_problem = QD.EnvelopeOptimizationProblem(
        builder,
        [variable],
        QD.UnitaryObjective(ComplexF64[0 1; 1 0]; states = [(0,), (1,)]),
    )
    unitary_result = QD.optimize_gate(unitary_problem; iterations = 40)
    @test unitary_result.gate isa QD.GateSpec
    @test unitary_result.converged
    @test unitary_result.fidelity > 0.999
    @test unitary_result.leakage < 1e-6
    @test unitary_result.minimizer[:amplitude] ≈ 0.5 atol = 2e-3
    @test 0.0 <= unitary_result.minimizer[:amplitude] <= 1.0
    @test isempty(unitary_result.failures)

    state_problem = QD.EnvelopeOptimizationProblem(
        builder,
        [variable],
        QD.StateTransferObjective((0,), (1,); states = [(0,), (1,)]),
    )
    state_result = QD.optimize_gate(state_problem; iterations = 40)
    @test state_result.gate isa QD.GateSpec
    @test state_result.fidelity > 0.999
    @test state_result.minimizer[:amplitude] ≈ 0.5 atol = 2e-3

    repeat_result = QD.optimize_gate(unitary_problem; iterations = 40)
    @test repeat_result.minimizer == unitary_result.minimizer
    @test repeat_result.loss == unitary_result.loss
end

@testset "invalid candidates are recorded" begin
    (; spec) = driven_qubit_fixture()
    variable = QD.GateVariable(:duration; initial = 0.5, lower = 0.1, upper = 1.0)
    problem = QD.EnvelopeOptimizationProblem(
        _ -> error("invalid candidate"),
        [variable],
        QD.StateTransferObjective((0,), (1,)),
    )
    result = QD.optimize_gate(problem; iterations = 2)
    @test result.gate === nothing
    @test !result.converged
    @test isinf(result.loss)
    @test result.evaluations == 4
    @test length(result.failures) == 4
end

@testset "model control decomposition" begin
    (; spec, path) = driven_qubit_fixture()
    decomposition = QD._model_control_decomposition(QD.model(spec), [path])
    @test length(decomposition.terms) == 1
    @test decomposition.terms[1].linear_weights == [1.0]
    @test decomposition.drift ≈ zeros(ComplexF64, 2, 2)

    nonlinear_drive = QD.DeviceParameter(
        :nonlinear;
        domain = QD.realdomain(),
        fixed = false,
        default = 0.0,
    )
    nonlinear_spec = QD.ModelSpec(
        :nonlinear_qubit,
        [QD.Component(QD.QubitSpec(0.0; fixed_frequency = true), :q)];
        interactions = (QD.param(nonlinear_drive)^2 * QD.op(:q, :x),),
        parameters = (nonlinear_drive,),
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    nonlinear_path = QD.ParamPath(:nonlinear)
    nonlinear_model = QD.model(nonlinear_spec)
    nonlinear = QD._model_control_decomposition(nonlinear_model, [nonlinear_path])
    @test nonlinear.terms[1].linear_weights === nothing
    @test norm(
        QD._control_decomposition_matrix(nonlinear, [0.3]) -
        Matrix(QD.hamiltonian(nonlinear_model; param = Dict(nonlinear_path => 0.3)).data),
    ) < 1e-10

    piccolo_problem = QD.PiccoloOptimizationProblem(
        :x_gate,
        spec,
        Dict(path => (-1.0, 1.0)),
        QD.UnitaryObjective(ComplexF64[0 1; 1 0]; states = [(0,), (1,)]);
        duration = 1.0,
        knots = 5,
    )
    @test size(piccolo_problem.initial_controls) == (1, 5)
    @test piccolo_problem.initial_controls[:, 1] == [0.0]
    @test piccolo_problem.initial_controls[:, end] == [0.0]
end
