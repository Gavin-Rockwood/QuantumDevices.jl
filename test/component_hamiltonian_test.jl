using Test
import LinearAlgebra: norm
import QuantumDevices as QD
import QuantumToolbox as qt

struct DimensionProbeSpec <: QD.AbstractComponentSpec
    dimension::QD.Dimension
    operators
    hamiltonian::QD.OperatorExpr
end

struct ParameterCollisionSpec <: QD.AbstractComponentSpec
    dimension::QD.Dimension
    operators
    hamiltonian::QD.OperatorExpr
end

struct OrphanParameterSpec <: QD.AbstractComponentSpec
    orphan::QD.DeviceParameter
    dimension::QD.Dimension
    operators
    hamiltonian::QD.OperatorExpr
end

struct CopyableParameterSpec <: QD.AbstractComponentSpec
    gain::QD.DeviceParameter
    dimension::QD.Dimension
    operators
    hamiltonian::QD.OperatorExpr
end

struct NonReconstructibleSpec <: QD.AbstractComponentSpec
    gain::QD.DeviceParameter
    dimension::QD.Dimension
    operators
    hamiltonian::QD.OperatorExpr

    function NonReconstructibleSpec()
        gain = QD.DeviceParameter(:gain; default = 1.0)
        new(
            gain,
            QD.Dimension(2),
            (x = (; kwargs...) -> qt.sigmax(),),
            QD.param(gain) * QD.op(:x),
        )
    end
end

DimensionProbeSpec() = DimensionProbeSpec(
    QD.Dimension(3),
    Dict(:a => (; dimension, kwargs...) -> fill(prod(dimension), 1, 1),
         :b => (; dimension, kwargs...) -> fill(2 * prod(dimension), 1, 1)),
    QD.op(:a) + QD.op(:b),
)

@testset "component Hamiltonians" begin
    built_component = QD.Component(QD.QubitSpec(5.0), :built_qubit)
    @test built_component.spec.frequency.default == 5.0
    @test QD.hamiltonian(built_component) === built_component.spec.hamiltonian
    @test norm(QD.numerical(built_component) - 2.5 * qt.sigmaz()) ≈ 0

    default_spec = QD.QubitSpec(5.0)
    qubit_spec = QD.QubitSpec(5.0;
        hamiltonian = QD.param(default_spec.frequency) * QD.op(:x) / 2)
    component = QD.Component(qubit_spec, :q1)

    @test component.name === :q1
    @test !hasproperty(component, :id)
    @test !hasproperty(component, :uuid)
    @test fieldnames(typeof(component)) == (:spec, :name)
    @test_throws MethodError QD.Component(qubit_spec, "q1")
    @test haskey(QD.parameters(component), QD.ParamPath(:frequency))
    @test !haskey(QD.parameters(component), QD.ParamPath(:dimension))
    @test component.spec.dimension == QD.Dimension(2)
    @test :x in QD.operators(component)
    @test QD.hamiltonian(component) isa QD.OperatorExpr
    @test !hasproperty(component, :hamiltonian_expr)

    @test norm(QD.numerical(component, :x; dimension = (2,)) - QD.numerical_operator(component.spec, :x; dimension = (2,))) ≈ 0
    @test norm(QD.numerical(component, :x) - QD.numerical_operator(component.spec, :x; dimension = (2,))) ≈ 0
    @test norm(QD.numerical(component; dimension = (2,)) - 2.5 * qt.sigmax()) ≈ 0
    @test norm(QD.numerical(component) - 2.5 * qt.sigmax()) ≈ 0
    @test norm(QD.numerical(component, QD.hamiltonian(component); dimension = (2,)) - 2.5 * qt.sigmax()) ≈ 0
    @test norm(QD.numerical(component, QD.hamiltonian(component)) - 2.5 * qt.sigmax()) ≈ 0
    @test norm(QD.numerical(component, QD.hamiltonian(component); dimension = (2,), frequency = 10.0) - 5.0 * qt.sigmax()) ≈ 0
    @test norm(QD.numerical(component, QD.hamiltonian(component); dimension = (2,), params = Dict(:frequency => 10.0)) - 5.0 * qt.sigmax()) ≈ 0
    @test norm(QD.numerical(component, QD.hamiltonian(component); dimension = (2,), frequency = 8.0, params = Dict(:frequency => 10.0)) - 4.0 * qt.sigmax()) ≈ 0

    probe = QD.Component(DimensionProbeSpec(), :probe)

    @test QD.numerical(probe; dimension = (7,)) == fill(21, 1, 1)
    @test QD.numerical(probe) == fill(9, 1, 1)
    @test_throws ArgumentError QD.numerical(component; dimension = 2)
    @test_throws Exception QD.numerical(component, :bad; dimension = (2,))
    @test_throws Exception QD.numerical(component, QD.op(:q2, :x); dimension = (2,))

    first_g = QD.DeviceParameter(:g; default = 1.0)
    second_g = QD.DeviceParameter(:g; default = 2.0)
    collision_spec = ParameterCollisionSpec(
        QD.Dimension(2),
        (x = (; kwargs...) -> qt.sigmax(),),
        QD.param(first_g) * QD.op(:x) + QD.param(second_g) * QD.op(:x),
    )
    @test QD.parameters(collision_spec.hamiltonian)[QD.ParamPath(:g)] === second_g
    collision_parameters = QD.parameters(QD.Component(collision_spec, :collision))
    @test collision_parameters[QD.ParamPath(:g)] === second_g

    orphan = QD.DeviceParameter(:orphan; default = 1.0)
    orphan_component = QD.Component(
        OrphanParameterSpec(
            orphan,
            QD.Dimension(2),
            (x = (; kwargs...) -> qt.sigmax(),),
            QD.op(:x),
        ),
        :orphan,
    )
    @test isempty(QD.parameters(orphan_component))
end

@testset "component copies" begin
    original = QD.Component(
        QD.QubitSpec(5.0; metadata = Dict(:origin => "original")),
        :q1,
    )
    frequency = QD.ParamPath(:frequency)

    same_name = copy(original)
    @test same_name.name === :q1
    @test same_name !== original
    @test same_name.spec !== original.spec
    @test QD.parameters(same_name)[frequency].default == 5.0
    @test QD.parameters(same_name)[frequency] !== QD.parameters(original)[frequency]
    @test same_name.spec.frequency === QD.parameters(same_name)[frequency]
    @test same_name.spec.frequency.metadata == original.spec.frequency.metadata
    @test same_name.spec.frequency.metadata !== original.spec.frequency.metadata

    renamed = copy(original; name = :q2, params = (frequency = 7.0,))
    @test renamed.name === :q2
    @test original.spec.frequency.default == 5.0
    @test renamed.spec.frequency.default == 7.0
    @test only(QD.get_params(QD.hamiltonian(renamed))) === renamed.spec.frequency
    @test norm(QD.numerical(original) - 2.5 * qt.sigmaz()) ≈ 0
    @test norm(QD.numerical(renamed) - 3.5 * qt.sigmaz()) ≈ 0

    @test copy(original; params = Dict(:frequency => 6.0)).spec.frequency.default == 6.0
    @test copy(original; params = Dict(frequency => 6.5)).spec.frequency.default == 6.5
    @test copy(original; params = Dict((:frequency,) => 6.75)).spec.frequency.default == 6.75
    @test copy(original; params = Dict("frequency" => 7.5)).spec.frequency.default == 7.5

    gain = QD.DeviceParameter(
        :gain;
        domain = QD.interval(0.0, 4.0),
        fixed = false,
        default = 1.0,
    )
    custom = QD.Component(
        CopyableParameterSpec(
            gain,
            QD.Dimension(2),
            (x = (; kwargs...) -> qt.sigmax(),),
            QD.param(gain) * QD.op(:x),
        ),
        :custom,
    )
    custom_copy = copy(custom; name = :custom_copy, params = Dict((:gain,) => 2.0))
    @test custom_copy.spec isa CopyableParameterSpec
    @test custom_copy.spec.gain.default == 2.0
    @test custom_copy.spec.gain !== gain
    @test only(QD.get_params(custom_copy.spec.hamiltonian)) === custom_copy.spec.gain
    @test norm(QD.numerical(custom_copy) - 2.0 * qt.sigmax()) ≈ 0

    shared = QD.DeviceParameter(:g; default = 1.0)
    repeated = QD.Component(
        ParameterCollisionSpec(
            QD.Dimension(2),
            (x = (; kwargs...) -> qt.sigmax(),),
            QD.param(shared) * QD.op(:x) + QD.param(shared) * QD.op(:x),
        ),
        :repeated,
    )
    repeated_copy = copy(repeated; params = Dict(:g => 3.0))
    @test length(QD.get_params(repeated_copy.spec.hamiltonian)) == 1
    @test only(QD.get_params(repeated_copy.spec.hamiltonian)).default == 3.0

    first_g = QD.DeviceParameter(:g; default = 1.0)
    second_g = QD.DeviceParameter(:g; default = 2.0)
    collision = QD.Component(
        ParameterCollisionSpec(
            QD.Dimension(2),
            (x = (; kwargs...) -> qt.sigmax(),),
            QD.param(first_g) * QD.op(:x) + QD.param(second_g) * QD.op(:x),
        ),
        :collision,
    )
    collision_copy = copy(collision; params = Dict(:g => 3.0))
    copied_collisions = QD.get_params(collision_copy.spec.hamiltonian)
    @test length(copied_collisions) == 2
    @test all(parameter -> parameter.default == 3.0, copied_collisions)
    @test copied_collisions[1] !== first_g
    @test copied_collisions[2] !== second_g
    @test QD.parameters(collision_copy)[QD.ParamPath(:g)] === copied_collisions[2]
    @test norm(QD.numerical(collision_copy) - 6.0 * qt.sigmax()) ≈ 0

    orphan = QD.DeviceParameter(:orphan; default = 1.0)
    orphan_component = QD.Component(
        OrphanParameterSpec(
            orphan,
            QD.Dimension(2),
            (x = (; kwargs...) -> qt.sigmax(),),
            QD.op(:x),
        ),
        :orphan,
    )
    @test copy(orphan_component).spec.orphan !== orphan

    duplicate_paths = Dict{Any,Any}(:frequency => 6.0, frequency => 7.0)
    @test_throws Exception copy(original; params = Dict(:missing => 1.0))
    @test_throws Exception copy(original; params = duplicate_paths)
    @test_throws Exception copy(custom; params = Dict(:gain => 5.0))
    @test_throws Exception copy(original; params = 1)
    @test_throws ArgumentError copy(original; name = "q2")

    reconstruction_error = try
        copy(QD.Component(NonReconstructibleSpec(), :locked))
        nothing
    catch exception
        exception
    end
    @test reconstruction_error isa ErrorException
    @test occursin(
        "must support construction from its fields",
        sprint(showerror, reconstruction_error),
    )
end

@testset "component default updates" begin
    qubit = QD.Component(QD.QubitSpec(5.0), :q)
    parameter = qubit.spec.frequency

    @test QD.update!(qubit; params = (frequency = 5.5,)) === qubit
    @test qubit.spec.frequency === parameter
    @test parameter.default == 5.5
    @test norm(QD.numerical(qubit) - 2.75 * qt.sigmaz()) ≈ 0

    QD.update!(qubit; params = Dict("frequency" => 6.0))
    @test parameter.default == 6.0
    @test QD.update!(qubit) === qubit
    duplicate_paths =
        Dict{Any,Any}(:frequency => 6.5, QD.ParamPath(:frequency) => 7.0)
    @test_throws Exception QD.update!(qubit; params = duplicate_paths)
    @test parameter.default == 6.0

    fixed = QD.Component(QD.QubitSpec(5.0; fixed_frequency = true), :fixed)
    QD.update!(fixed; params = Dict(QD.ParamPath(:frequency) => 5.25))
    @test fixed.spec.frequency.default == 5.25
    @test_throws Exception QD.update!(
        fixed;
        params = Dict(:frequency => t -> 5.0 + t),
    )
    @test fixed.spec.frequency.default == 5.25

    transmon = QD.Component(QD.TransmonSpec(0.2, 20.0), :transmon)
    @test_throws Exception QD.update!(
        transmon;
        params = (EC = 0.3, EJ = -1.0),
    )
    @test transmon.spec.EC.default == 0.2
    @test transmon.spec.EJ.default == 20.0
    @test_throws Exception QD.update!(transmon; params = Dict(:missing => 1.0))
    @test_throws Exception QD.update!(transmon; params = 1)

    first_g = QD.DeviceParameter(
        :g;
        domain = QD.interval(0.0, 4.0),
        default = 1.0,
    )
    second_g = QD.DeviceParameter(
        :g;
        domain = QD.interval(0.0, 2.0),
        default = 1.5,
    )
    collision = QD.Component(
        ParameterCollisionSpec(
            QD.Dimension(2),
            (x = (; kwargs...) -> qt.sigmax(),),
            QD.param(first_g) * QD.op(:x) + QD.param(second_g) * QD.op(:x),
        ),
        :collision,
    )
    QD.update!(collision; params = Dict(:g => 2.0))
    @test first_g.default == 2.0
    @test second_g.default == 2.0
    @test_throws Exception QD.update!(collision; params = Dict(:g => 3.0))
    @test first_g.default == 2.0
    @test second_g.default == 2.0
end
