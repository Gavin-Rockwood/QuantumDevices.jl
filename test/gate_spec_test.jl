using Test
import LinearAlgebra: norm
import QuantumDevices as QD
import QuantumToolbox as qt

@testset "dynamic operator expressions" begin
    qubit = QD.Component(QD.QubitSpec(5.0), :q)
    amplitude = QD.DeviceParameter(:amplitude;
        domain = QD.realdomain(), fixed = false, default = t -> sin(t))
    expression = sin(QD.param(amplitude)) * (QD.op(:x) + QD.op(:z))^2

    static = QD.numerical(qubit, QD.op(:x))
    dynamic = QD.numerical(qubit, expression)
    expected_operator = (qt.sigmax() + qt.sigmaz())^2
    @test static isa qt.QuantumObject
    @test dynamic isa qt.QobjEvo
    @test norm(dynamic(nothing, 0.4) - sin(sin(0.4)) * expected_operator) < 1e-12

    phase = QD.DeviceParameter(:phase; fixed = false, default = t -> 1 + 1im * t)
    nonlinear = QD.numerical(qubit,
        conj(QD.param(phase)) / QD.param(phase)^2 * QD.op(:x))
    phase_value = 1 + 0.3im
    @test norm(nonlinear(nothing, 0.3) -
        conj(phase_value) / phase_value^2 * qt.sigmax()) < 1e-12

    dynamic_qubit = QD.Component(QD.QubitSpec(t -> 5.0 + sin(t)^2), :q)
    dynamic_spec = QD.ModelSpec(:dynamic_model, [dynamic_qubit];
        dressingspec = QD.DressingSpec(minimum_overlap = 0))
    dynamic_model = QD.model(dynamic_spec)
    model_expression = QD.param(dynamic_qubit.spec.frequency) * QD.op(:q, :z) / 2
    model_evolution = QD.numerical(dynamic_model, model_expression)
    @test model_evolution isa qt.QobjEvo
    @test norm(model_evolution(nothing, 0.4) -
        (5.0 + sin(0.4)^2) * qt.sigmaz() / 2) < 1e-12

    @test_throws Exception QD.DeviceParameter(:fixed_trajectory;
        fixed = true, default = t -> t)
    @test_throws Exception QD.QubitSpec(t -> 5.0; fixed_frequency = true)
    @test_throws Exception QD.DeviceParameter(:outside;
        domain = QD.positivedomain(), default = -1.0)

    bounded = QD.DeviceParameter(:bounded;
        domain = QD.interval(-1.0, 1.0), fixed = false, default = t -> 2t)
    bounded_evolution = QD.numerical(qubit, QD.param(bounded) * QD.op(:x))
    @test_throws Exception bounded_evolution(nothing, 1.0)
end

@testset "GateSpec Hamiltonians" begin
    duration = 2.0
    flux_default = 0.125
    flux_baseline = 0.2
    coupler = QD.Component(
        QD.FluxTunableTransmonSpec(0.25, 12.0, 8.0; flux = flux_default),
        :coupler,
    )
    flux_spec = QD.ModelSpec(
        :flux_model,
        [coupler],
        Dict(:coupler => 3);
        pretruncation_dims = Dict(:coupler => 5),
        defaults = Dict((:coupler, :flux) => flux_baseline),
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    flux_model = QD.model(flux_spec)
    flux = QD.DeviceParameter((:coupler, :flux);
        domain = QD.realdomain(),
        fixed = false,
        default = t -> flux_baseline + 0.05 * sin(pi * t / duration)^2,
    )
    flux_gate = QD.GateSpec(:flux_pulse, flux_spec;
        duration, parameters = (flux,))
    flux_hamiltonian = QD.numerical(flux_model, flux_gate)

    @test flux_hamiltonian isa qt.QobjEvo
    @test norm(flux_hamiltonian(nothing, 0.0) - flux_model.hamiltonian) < 1e-12
    @test norm(flux_hamiltonian(nothing, duration / 2) - flux_model.hamiltonian) > 1e-8
    @test norm(flux_hamiltonian(nothing, duration) - flux_model.hamiltonian) < 1e-12
    @test QD.numerical(flux_gate) isa qt.QobjEvo

    fixed_EC = QD.DeviceParameter((:coupler, :EC);
        domain = QD.positivedomain(), default = 0.3)
    @test_throws Exception QD.GateSpec(:fixed, flux_spec;
        duration, parameters = (fixed_EC,))

    bad_start = QD.DeviceParameter((:coupler, :flux);
        domain = QD.realdomain(), fixed = false,
        default = t -> flux_baseline + 0.01 + 0.05 * sin(pi * t / duration)^2)
    @test_throws Exception QD.numerical(flux_model,
        QD.GateSpec(:bad_start, flux_spec; duration, parameters = (bad_start,)))

    bad_end = QD.DeviceParameter((:coupler, :flux);
        domain = QD.realdomain(), fixed = false,
        default = t -> flux_baseline + 0.05 * t / duration)
    @test_throws Exception QD.numerical(flux_model,
        QD.GateSpec(:bad_end, flux_spec; duration, parameters = (bad_end,)))

    other_spec = QD.ModelSpec(:other, [QD.Component(QD.QubitSpec(5.0), :q)];
        dressingspec = QD.DressingSpec(minimum_overlap = 0))
    @test_throws Exception QD.numerical(QD.model(other_spec), flux_gate)
end

@testset "gate drives and interaction tuning" begin
    duration = 1.0
    qubit = QD.Component(QD.QubitSpec(5.0; fixed_frequency = true), :q)
    drive_spec = QD.ModelSpec(:drive_model, [qubit];
        dressingspec = QD.DressingSpec(minimum_overlap = 0))
    drive_model = QD.model(drive_spec)
    drive = QD.DeviceParameter(:drive;
        domain = QD.realdomain(), fixed = false,
        default = t -> 2.0 * sin(pi * t / duration))
    drive_gate = QD.GateSpec(:x_drive, drive_spec;
        duration,
        hamiltonian = QD.param(drive) * QD.op(:q, :x),
        parameters = (drive,),
    )
    @test fieldnames(typeof(drive_gate)) ==
        (:name, :modelspec, :duration, :hamiltonian, :parameters)
    @test all(value -> !(value isa Union{qt.QuantumObject,qt.QobjEvo}),
        (getfield(drive_gate, field) for field in fieldnames(typeof(drive_gate))))
    drive_hamiltonian = QD.numerical(drive_model, drive_gate)
    @test drive_hamiltonian isa qt.QobjEvo
    @test norm(drive_hamiltonian(nothing, 0.5) -
        (drive_model.hamiltonian + 2.0 * drive_model.operators[(:q, :x)])) < 1e-12

    idle = QD.GateSpec(:idle, drive_spec; duration)
    @test QD.numerical(drive_model, idle) isa qt.QuantumObject

    missing = QD.DeviceParameter(:missing; fixed = false, default = t -> 0.0)
    @test_throws Exception QD.GateSpec(:missing, drive_spec;
        duration, hamiltonian = QD.param(missing) * QD.op(:q, :x))

    q1 = QD.Component(QD.QubitSpec(5.0), :q1)
    q2 = QD.Component(QD.QubitSpec(5.1), :q2)
    g = QD.DeviceParameter(:g;
        domain = QD.realdomain(), fixed = false, default = 0.1)
    interaction = QD.param(g) * QD.op(:q1, :x) * QD.op(:q2, :x)
    pair_spec = QD.ModelSpec(:pair, [q1, q2];
        interactions = (interaction,), parameters = (g,),
        dressingspec = QD.DressingSpec(minimum_overlap = 0))
    pair_model = QD.model(pair_spec)
    gate_g = QD.DeviceParameter(:g;
        domain = QD.realdomain(), fixed = false,
        default = t -> 0.1 + 0.05 * sin(pi * t / duration)^2)
    interaction_gate = QD.GateSpec(:tunable_interaction, pair_spec;
        duration, parameters = (gate_g,))
    interaction_hamiltonian = QD.numerical(pair_model, interaction_gate)
    @test interaction_hamiltonian isa qt.QobjEvo
    @test norm(interaction_hamiltonian(nothing, 0.5) - pair_model.hamiltonian) > 1e-8
    @test norm(interaction_hamiltonian(nothing, duration) - pair_model.hamiltonian) < 1e-12
end
