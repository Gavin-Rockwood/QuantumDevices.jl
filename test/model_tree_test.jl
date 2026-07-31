using Test
import LinearAlgebra: norm
import QuantumDevices as QD
import QuantumToolbox as qt

struct LadderSpec <: QD.AbstractComponentSpec
    dimension::QD.Dimension
    operators::Dict{Any,Any}
    hamiltonian::QD.OperatorExpr
end

function LadderSpec(dimension::Int)
    operators = Dict{Any,Any}(
        :n => (; dimension, kwargs...) -> qt.num(only(dimension)),
        :x => (; dimension, kwargs...) -> qt.destroy(only(dimension)) + qt.create(only(dimension)),
    )
    LadderSpec(QD.Dimension(dimension), operators, QD.op(:n))
end

parameter(name, default) = QD.DeviceParameter(name;
    domain = QD.realdomain(), default)

@testset "flat model behavior" begin
    q1 = QD.Component(LadderSpec(3), :q1)
    q2 = QD.Component(LadderSpec(3), :q2)
    g = parameter(:g, 0.1)
    spare = parameter(:spare, 2.0)
    interaction = QD.param(g) * QD.op(:q1, :x) * QD.op(:q2, :x)

    @test_throws Exception QD.ModelSpec(:missing, [q1, q2]; interactions = (interaction,))
    @test_throws Exception QD.ModelSpec(:duplicate, [q1, QD.Component(LadderSpec(2), :q1)])
    @test_throws Exception QD.ModelSpec(:bad_default, [q1];
        defaults = Dict((:q1, :missing) => 1.0))
    @test_throws Exception QD.ModelSpec(:bad_dimension, [q1], Dict(:missing => 2))

    qubit = QD.Component(QD.QubitSpec(5.0), :qubit)
    default_spec = QD.ModelSpec(:defaulted, [qubit];
        defaults = Dict((:qubit, :frequency) => 7.0),
        dressingspec = QD.DressingSpec(minimum_overlap = 0))
    @test default_spec.parameters[QD.ParamPath(:qubit, :frequency)] === qubit.spec.frequency
    @test default_spec.defaults == Dict{Tuple{Vararg{Symbol}},Any}(
        (:qubit, :frequency) => 7.0,
    )
    @test sort(collect(values(QD.model(default_spec).energies))) ≈ [-3.5, 3.5]

    resonator = QD.Component(QD.ResonatorSpec(5.0; dimension = 3), :resonator)
    @test_throws Exception QD.ModelSpec(:bad_domain, [resonator];
        defaults = Dict((:resonator, :frequency) => -1.0))

    spec = QD.ModelSpec(
        :pair,
        [q2, q1],
        Dict(:q1 => 2, :q2 => 2);
        interactions = (interaction,),
        parameters = (g, spare),
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    @test keys(spec.subsystems) == (:q2, :q1)
    @test spec.dimension == QD.Dimension((2, 2))
    @test spec.interactions == (interaction,)
    @test spec.parameters[QD.ParamPath(:g)] === g
    @test spec.parameters[QD.ParamPath(:spare)] === spare
    @test haskey(spec.parameters, QD.ParamPath(:q1, :n)) == false
    @test spec.defaults == Dict{Tuple{Vararg{Symbol}},Any}(
        (:g,) => 0.1,
        (:spare,) => 2.0,
    )

    overridden_model_defaults = QD.ModelSpec(
        :overridden_model_defaults,
        [q2, q1],
        Dict(:q1 => 2, :q2 => 2);
        interactions = (interaction,),
        parameters = (g, spare),
        defaults = Dict((:g,) => 0.2, QD.ParamPath(:spare) => 3.0),
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    @test overridden_model_defaults.defaults[(:g,)] == 0.2
    @test overridden_model_defaults.defaults[(:spare,)] == 3.0

    built = QD.model(spec)
    @test size(built.hamiltonian) == (4, 4)
    @test built.spec === spec
    @test built.dressing_res isa QD.StateTrackingResult
    @test (:q1, :x) in keys(built.operators)
    @test norm(QD.hamiltonian(built).data - QD.hamiltonian(QD.model(spec)).data) < 1e-12

    selected = QD.ModelSpec(
        :selected,
        [q1, q2],
        Dict(:q1 => 2, :q2 => 2);
        states_to_keep = [(0, 0), (1, 0)],
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    selected_model = QD.model(selected)
    @test size(selected_model.hamiltonian) == (2, 2)
    @test Set(keys(selected_model.states)) == Set([(0, 0), (1, 0)])

    transmon = QD.Component(QD.TransmonSpec(0.2, 20.0), :transmon)
    transmon_spec = QD.ModelSpec(
        :transmon_model,
        [transmon],
        Dict(:transmon => 3);
        pretruncation_dims = Dict(:transmon => 5),
        defaults = Dict((:transmon, :ng) => 0.1),
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    @test transmon_spec.subsystems.transmon.spec.dimension == QD.Dimension(5)
    @test transmon_spec.subsystems.transmon.spec isa QD.GenericSpec
    @test transmon_spec.subsystems.transmon.spec.source === transmon.spec
    @test transmon_spec.defaults[(:transmon, :ng)] == 0.1
    @test transmon_spec.defaults[(:transmon, :EC)] == 0.2
    transmon_before_update = QD.hamiltonian(QD.model(transmon_spec))
    QD.update!(transmon; params = (EC = 0.3,))
    @test transmon_spec.defaults[(:transmon, :EC)] == 0.2
    @test norm(
        QD.hamiltonian(QD.model(transmon_spec)) - transmon_before_update,
    ) < 1e-12
    @test_throws Exception QD.ModelSpec(:unbounded, [transmon], Dict(:transmon => 3))

    effective = QD.component(selected_model; name = :effective)
    nested_spec = QD.ModelSpec(:nested, [effective];
        dressingspec = QD.DressingSpec(minimum_overlap = 0))
    @test effective.name == :effective
    @test effective.spec isa QD.GenericSpec
    @test effective.spec.source === selected_model.spec
    @test effective.spec.source.subsystems.q1.spec === q1.spec
    @test size(QD.model(nested_spec).hamiltonian) == (2, 2)
end

@testset "ModelSpec default snapshots" begin
    qubit = QD.Component(QD.QubitSpec(5.0), :q)
    old_spec = QD.ModelSpec(
        :old,
        [qubit];
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    old_hamiltonian = QD.hamiltonian(QD.model(old_spec))

    QD.update!(qubit; params = (frequency = 6.0,))
    new_spec = QD.ModelSpec(
        :new,
        [qubit];
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )

    @test old_spec.parameters[QD.ParamPath(:q, :frequency)] ===
        qubit.spec.frequency
    @test old_spec.defaults[(:q, :frequency)] == 5.0
    @test new_spec.defaults[(:q, :frequency)] == 6.0
    @test norm(QD.hamiltonian(QD.model(old_spec)) - old_hamiltonian) < 1e-12
    @test norm(
        QD.hamiltonian(QD.model(new_spec)) - 3.0 * qt.sigmaz(),
    ) < 1e-12
end

@testset "parameterized model Hamiltonians" begin
    q1 = QD.Component(QD.QubitSpec(5.0), :q1)
    q2 = QD.Component(QD.QubitSpec(5.1; fixed_frequency = true), :q2)
    g = QD.DeviceParameter(:g;
        domain = QD.interval(-1.0, 1.0),
        fixed = false,
        default = 0.1)
    spare = QD.DeviceParameter(:spare; fixed = false, default = 2.0)
    interaction = QD.param(g) * QD.op(:q1, :x) * QD.op(:q2, :x)
    spec = QD.ModelSpec(
        :parameterized,
        [q1, q2];
        interactions = (interaction,),
        parameters = (g, spare),
        defaults = Dict((:q1, :frequency) => 6.0),
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    built = QD.model(spec)
    frequency_path = QD.ParamPath("q1/frequency")
    fixed_path = QD.ParamPath(:q2, :frequency)
    g_path = QD.ParamPath(:g)
    expected = 3.0 * built.operators[(:q1, :z)] +
        2.55 * built.operators[(:q2, :z)] +
        0.1 * built.operators[(:q1, :x)] * built.operators[(:q2, :x)]

    @test built.hamiltonian isa qt.QobjEvo
    @test Set(keys(QD.parameters(built))) == Set((frequency_path, g_path))
    @test QD.parameters(built)[frequency_path] === q1.spec.frequency
    @test QD.parameters(built)[g_path] === g
    @test spec.defaults == Dict{Tuple{Vararg{Symbol}},Any}(
        (:q1, :frequency) => 6.0,
        (:q2, :frequency) => 5.1,
        (:g,) => 0.1,
        (:spare,) => 2.0,
    )
    @test norm(QD.hamiltonian(built) - expected) < 1e-12
    @test norm(built.hamiltonian(Dict{QD.ParamPath,Any}(), 0.0) - expected) < 1e-12

    g.default = 0.3
    @test spec.defaults[(:g,)] == 0.1
    @test norm(QD.hamiltonian(QD.model(spec)) - expected) < 1e-12

    overridden = QD.hamiltonian(built;
        param = Dict(frequency_path => 7.0, g_path => 0.2))
    expected_override = 3.5 * built.operators[(:q1, :z)] +
        2.55 * built.operators[(:q2, :z)] +
        0.2 * built.operators[(:q1, :x)] * built.operators[(:q2, :x)]
    @test norm(overridden - expected_override) < 1e-12
    @test norm(built.hamiltonian(Dict(g_path => 0.2), 0.0) -
        QD.hamiltonian(built; param = Dict(g_path => 0.2))) < 1e-12
    @test norm(QD.hamiltonian(built;
        param = Dict(g_path => t -> 0.1 + t), t = 0.5) -
        (expected + 0.5 * built.operators[(:q1, :x)] *
            built.operators[(:q2, :x)])) < 1e-12

    @test_throws Exception QD.hamiltonian(built; param = (g = 0.2,))
    @test_throws Exception QD.hamiltonian(built; param = Dict(:g => 0.2))
    @test_throws Exception QD.hamiltonian(built;
        param = Dict(QD.ParamPath(:missing) => 0.2))
    @test_throws Exception QD.hamiltonian(built; param = Dict(fixed_path => 5.2))
    @test_throws Exception QD.hamiltonian(built; param = Dict(g_path => 2.0))
    @test_throws Exception built.hamiltonian(Dict(QD.ParamPath(:missing) => 0.2), 0.0)

    selected = QD.model(QD.ModelSpec(
        :parameterized_selected,
        [q1, q2];
        interactions = (interaction,),
        parameters = (g,),
        states_to_keep = [(0, 0), (1, 0)],
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    ))
    @test selected.hamiltonian isa qt.QobjEvo
    @test size(QD.hamiltonian(selected)) == (2, 2)
    @test norm(QD.hamiltonian(selected; param = Dict(g_path => 0.2)) -
        QD.hamiltonian(selected)) > 1e-8

    fixed = QD.model(QD.ModelSpec(
        :fixed,
        [QD.Component(QD.QubitSpec(5.0; fixed_frequency = true), :q)];
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    ))
    @test fixed.hamiltonian isa qt.QuantumObject
    @test isempty(QD.parameters(fixed))
    @test QD.hamiltonian(fixed) === fixed.hamiltonian
    @test_throws Exception QD.hamiltonian(fixed;
        param = Dict(QD.ParamPath(:q, :frequency) => 5.2))
end
