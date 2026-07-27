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

    built = QD.model(spec)
    @test size(built.hamiltonian) == (4, 4)
    @test built.spec === spec
    @test built.dressing_res isa QD.StateTrackingResult
    @test (:q1, :x) in keys(built.operators)
    @test norm(built.hamiltonian.data - QD.model(spec).hamiltonian.data) < 1e-12

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
    @test size(QD.model(transmon_spec).hamiltonian) == (3, 3)
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
