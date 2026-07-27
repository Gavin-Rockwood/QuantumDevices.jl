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
