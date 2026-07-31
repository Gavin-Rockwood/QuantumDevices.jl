using Test
import LinearAlgebra: norm
import QuantumDevices as QD
import QuantumToolbox as qt

@testset "dynamic operator expressions" begin
    qubit = QD.Component(QD.QubitSpec(5.0), :q)
    amplitude = QD.DeviceParameter(
        :amplitude;
        domain = QD.realdomain(),
        fixed = false,
        default = t -> sin(t),
    )
    expression = sin(QD.param(amplitude)) * (QD.op(:x) + QD.op(:z))^2

    static = QD.numerical(qubit, QD.op(:x))
    dynamic = QD.numerical(qubit, expression)
    expected_operator = (qt.sigmax() + qt.sigmaz())^2
    @test static isa qt.QuantumObject
    @test dynamic isa qt.QobjEvo
    @test norm(dynamic(nothing, 0.4) - sin(sin(0.4)) * expected_operator) < 1e-12

    phase = QD.DeviceParameter(:phase; fixed = false, default = t -> 1 + 1im * t)
    nonlinear = QD.numerical(qubit, conj(QD.param(phase)) / QD.param(phase)^2 * QD.op(:x))
    phase_value = 1 + 0.3im
    @test norm(nonlinear(nothing, 0.3) - conj(phase_value) / phase_value^2 * qt.sigmax()) <
          1e-12

    dynamic_qubit = QD.Component(QD.QubitSpec(t -> 5.0 + sin(t)^2), :q)
    dynamic_spec = QD.ModelSpec(
        :dynamic_model,
        [dynamic_qubit];
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    dynamic_model = QD.model(dynamic_spec)
    model_expression = QD.param(dynamic_qubit.spec.frequency) * QD.op(:q, :z) / 2
    model_evolution = QD.numerical(dynamic_model, model_expression)
    @test dynamic_model.hamiltonian isa qt.QobjEvo
    @test Set(keys(QD.parameters(dynamic_model))) == Set([QD.ParamPath(:q, :frequency)])
    @test norm(
        QD.hamiltonian(dynamic_model; t = 0.4) - (5.0 + sin(0.4)^2) * qt.sigmaz() / 2,
    ) < 1e-12
    @test model_evolution isa qt.QobjEvo
    @test norm(model_evolution(nothing, 0.4) - (5.0 + sin(0.4)^2) * qt.sigmaz() / 2) < 1e-12

    @test_throws Exception QD.DeviceParameter(
        :fixed_trajectory;
        fixed = true,
        default = t -> t,
    )
    @test_throws Exception QD.QubitSpec(t -> 5.0; fixed_frequency = true)
    @test_throws Exception QD.DeviceParameter(
        :outside;
        domain = QD.positivedomain(),
        default = -1.0,
    )

    bounded = QD.DeviceParameter(
        :bounded;
        domain = QD.interval(-1.0, 1.0),
        fixed = false,
        default = t -> 2t,
    )
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
    flux = QD.DeviceParameter(
        (:coupler, :flux);
        domain = QD.realdomain(),
        fixed = false,
        default = t -> flux_baseline + 0.05 * sin(pi * t / duration)^2,
    )
    flux_gate = QD.GateSpec(
        :flux_pulse,
        flux_spec;
        duration,
        controls = Dict(QD.ParamPath(:coupler, :flux) => flux.default),
    )
    flux_hamiltonian = QD.numerical(flux_model, flux_gate)

    @test flux_hamiltonian isa qt.QobjEvo
    @test norm(flux_hamiltonian(nothing, 0.0) - QD.hamiltonian(flux_model)) < 1e-12
    @test norm(flux_hamiltonian(nothing, duration / 2) - QD.hamiltonian(flux_model)) > 1e-8
    @test norm(flux_hamiltonian(nothing, duration) - QD.hamiltonian(flux_model)) < 1e-12
    @test QD.numerical(flux_gate) isa qt.QobjEvo

    QD.update!(coupler; params = (flux = 0.3,))
    rebuilt_flux_model = QD.model(flux_spec)
    rebuilt_flux_gate = QD.numerical(rebuilt_flux_model, flux_gate)
    @test flux_spec.defaults[(:coupler, :flux)] == flux_baseline
    @test norm(QD.hamiltonian(rebuilt_flux_model) - QD.hamiltonian(flux_model)) < 1e-12
    @test norm(rebuilt_flux_gate(nothing, 0.0) - QD.hamiltonian(rebuilt_flux_model)) < 1e-12

    fixed_EC =
        QD.DeviceParameter((:coupler, :EC); domain = QD.positivedomain(), default = 0.3)
    @test_throws Exception QD.GateSpec(
        :fixed,
        flux_spec;
        duration,
        controls = Dict(fixed_EC.path => fixed_EC.default),
    )

    @test_throws Exception QD.GateSpec(
        :bad_start,
        flux_spec;
        duration,
        controls = Dict(
            QD.ParamPath(:coupler, :flux) =>
                (t -> flux_baseline + 0.01 + 0.05 * sin(pi * t / duration)^2),
        ),
    )
    @test_throws Exception QD.GateSpec(
        :bad_end,
        flux_spec;
        duration,
        controls = Dict(
            QD.ParamPath(:coupler, :flux) => (t -> flux_baseline + 0.05 * t / duration),
        ),
    )

    other_spec = QD.ModelSpec(
        :other,
        [QD.Component(QD.QubitSpec(5.0), :q)];
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    @test_throws Exception QD.numerical(QD.model(other_spec), flux_gate)
end

@testset "gate drives and interaction tuning" begin
    duration = 1.0
    qubit = QD.Component(QD.QubitSpec(5.0; fixed_frequency = true), :q)
    drive =
        QD.DeviceParameter(:drive; domain = QD.realdomain(), fixed = false, default = 0.0)
    drive_expression = QD.param(drive) * QD.op(:q, :x)
    drive_spec = QD.ModelSpec(
        :drive_model,
        [qubit];
        interactions = (drive_expression,),
        parameters = (drive,),
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    drive_model = QD.model(drive_spec)
    drive_control = t -> 2.0 * sin(pi * t / duration)
    drive_gate = QD.GateSpec(
        :x_drive,
        drive_spec;
        duration,
        controls = Dict(QD.ParamPath(:drive) => drive_control),
    )
    @test fieldnames(typeof(drive_gate)) ==
          (:name, :modelspec, :duration, :controls, :parameters, :control_recipes)
    @test all(
        value -> !(value isa Union{qt.QuantumObject,qt.QobjEvo}),
        (getfield(drive_gate, field) for field in fieldnames(typeof(drive_gate))),
    )
    drive_hamiltonian = QD.numerical(drive_model, drive_gate)
    @test drive_hamiltonian isa qt.QobjEvo
    @test norm(
        drive_hamiltonian(nothing, 0.5) - QD.hamiltonian(
            drive_model;
            param = Dict(QD.ParamPath(:drive) => drive_control),
            t = 0.5,
        ),
    ) < 1e-12

    idle = QD.GateSpec(:idle, drive_spec; duration)
    @test norm(QD.numerical(drive_model, idle) - QD.hamiltonian(drive_model)) < 1e-12

    @test_throws Exception QD.GateSpec(
        :missing,
        drive_spec;
        duration,
        controls = Dict(QD.ParamPath(:missing) => (t -> 0.0)),
    )
    @test_throws Exception QD.GateSpec(
        :bad_domain,
        drive_spec;
        duration,
        controls = Dict(QD.ParamPath(:drive) => (t -> t == 0 ? 0.0 : "bad")),
    )

    q1 = QD.Component(QD.QubitSpec(5.0), :q1)
    q2 = QD.Component(QD.QubitSpec(5.1), :q2)
    g = QD.DeviceParameter(:g; domain = QD.realdomain(), fixed = false, default = 0.1)
    interaction = QD.param(g) * QD.op(:q1, :x) * QD.op(:q2, :x)
    pair_spec = QD.ModelSpec(
        :pair,
        [q1, q2];
        interactions = (interaction,),
        parameters = (g,),
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    pair_model = QD.model(pair_spec)
    gate_g = QD.DeviceParameter(
        :g;
        domain = QD.realdomain(),
        fixed = false,
        default = t -> 0.1 + 0.05 * sin(pi * t / duration)^2,
    )
    interaction_gate = QD.GateSpec(
        :tunable_interaction,
        pair_spec;
        duration,
        controls = Dict(QD.ParamPath(:g) => gate_g.default),
    )
    interaction_hamiltonian = QD.numerical(pair_model, interaction_gate)
    @test interaction_hamiltonian isa qt.QobjEvo
    @test norm(interaction_hamiltonian(nothing, 0.5) - QD.hamiltonian(pair_model)) > 1e-8
    @test norm(interaction_hamiltonian(nothing, duration) - QD.hamiltonian(pair_model)) <
          1e-12
end

@testset "parameterized GateSpec recipes" begin
    (; spec, path) = let
        qubit = QD.Component(QD.QubitSpec(0.0; fixed_frequency = true), :q)
        drive = QD.DeviceParameter(
            :drive;
            domain = QD.realdomain(),
            fixed = false,
            default = 0.0,
        )
        spec = QD.ModelSpec(
            :recipe_qubit,
            [qubit];
            interactions = (QD.param(drive) * QD.op(:q, :x),),
            parameters = (drive,),
            dressingspec = QD.DressingSpec(minimum_overlap = 0),
        )
        (; spec, path = QD.ParamPath(:drive))
    end
    recipe = (p, time) -> p.amplitude * sin(π * time / p.duration)^2
    gate = QD.GateSpec(
        :parameterized_x,
        spec;
        duration = 1.0,
        parameters = (amplitude = 0.2, sigma = 0.1),
        controls = Dict(path => recipe),
    )
    @test gate.parameters == (amplitude = 0.2, sigma = 0.1)
    @test gate.controls[path](0.5) ≈ 0.2

    candidate = QD.with_parameters(gate; amplitude = 0.4, duration = 2.0)
    @test candidate !== gate
    @test candidate.duration == 2.0
    @test candidate.parameters.amplitude == 0.4
    @test candidate.controls[path](1.0) ≈ 0.4
    @test gate.duration == 1.0
    @test gate.parameters.amplitude == 0.2
    @test gate.controls[path](0.5) ≈ 0.2

    time_only = QD.GateSpec(
        :time_only,
        spec;
        duration = 1.0,
        parameters = (amplitude = 0.2,),
        controls = Dict(path => (time -> 0.2 * sin(π * time)^2)),
    )
    @test time_only.controls[path](0.5) ≈ 0.2

    ambiguous(args...) = 0.0
    @test_throws Exception QD.GateSpec(
        :ambiguous,
        spec;
        duration = 1.0,
        parameters = (amplitude = 0.2,),
        controls = Dict(path => ambiguous),
    )
    @test_throws Exception QD.GateSpec(
        :incompatible,
        spec;
        duration = 1.0,
        parameters = (amplitude = 0.2,),
        controls = Dict(path => ((_, _, _) -> 0.0)),
    )
    @test_throws Exception QD.GateSpec(
        :reserved,
        spec;
        duration = 1.0,
        parameters = (duration = 1.0,),
    )
    @test_throws Exception QD.with_parameters(gate; missing = 1.0)
end

@testset "laboratory-frame carrier controls" begin
    qubit = QD.Component(QD.QubitSpec(1.0; fixed_frequency = true), :q)
    drive = QD.DeviceParameter(
        :drive;
        domain = QD.realdomain(),
        fixed = false,
        default = 0.0,
    )
    spec = QD.ModelSpec(
        :carrier_qubit,
        [qubit];
        interactions = (QD.param(drive) * QD.op(:q, :x),),
        parameters = (drive,),
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    path = QD.ParamPath(:drive)
    carrier = QD.CarrierControl(
        (p, time) -> p.amplitude * sinpi(time / p.duration)^2;
        frequency = 1.0,
        phase = 0.0,
    )
    gate = QD.GateSpec(
        :carrier_x,
        spec;
        duration = 1.0,
        parameters = (amplitude = 0.2,),
        controls = Dict(path => carrier),
    )
    @test gate.controls[path](0.0) == 0.0
    @test gate.controls[path](0.5) ≈ -0.2
    @test gate.controls[path](1.0) ≈ 0.0 atol = 1e-14
    @test QD._carrier_envelope(gate, path).envelope(0.5) ≈ 0.2
    @test QD._carrier_envelope(gate, path).modulation(0.5) ≈ -1.0

    candidate = QD.with_parameters(gate; amplitude = 0.4)
    @test candidate.controls[path](0.5) ≈ -0.4
    built = QD.model(spec)
    direct = QD.hamiltonian(
        built;
        param = Dict(path => gate.controls[path]),
        t = 0.37,
    )
    @test norm(QD.numerical(built, gate)(nothing, 0.37) - direct) < 1e-12

    @test_throws Exception QD.CarrierControl(time -> time; frequency = 0.0)
    @test_throws Exception QD.CarrierControl(time -> time; frequency = Inf)
    @test_throws Exception QD.CarrierControl(time -> time; frequency = 1.0, phase = Inf)
    @test_throws Exception QD.GateSpec(
        :scalar_envelope,
        spec;
        duration = 1.0,
        parameters = NamedTuple(),
        controls = Dict(path => QD.CarrierControl(0.0; frequency = 1.0)),
    )
end
