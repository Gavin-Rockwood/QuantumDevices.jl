using Test
import LinearAlgebra: norm
import QuantumDevices as QD
import QuantumToolbox as qt

struct DimensionProbeSpec <: QD.AbstractComponentSpec
    dimension::QD.DeviceParameter
end

DimensionProbeSpec() = DimensionProbeSpec(
    QD.DeviceParameter(
        QD.ParamPath(:dimension),
        QD.integerrange(3, 3),
        true,
        true,
        3,
        "Probe Dimension",
        Dict{Symbol, Any}(),
    ),
)

QD.available_operators(::DimensionProbeSpec) = Set([:a, :b])

function QD.numerical_operator(::DimensionProbeSpec, operator; dimension = nothing, kwargs...)
    if operator == :a
        return fill(dimension, 1, 1)
    elseif operator == :b
        return fill(2 * dimension, 1, 1)
    else
        error("Unknown probe operator $operator.")
    end
end

@testset "component Hamiltonians" begin
    built_component = QD.Component(:built_qubit, QD.QubitSpec, 5.0)
    @test built_component.spec.frequency.default == 5.0
    @test built_component.hamiltonian === built_component.spec.hamiltonian
    @test norm(QD.numerical(built_component) - 2.5 * qt.sigmaz()) ≈ 0

    component = QD.Component(
        QD.QubitSpec(5.0),
        :q1;
        hamiltonian = QD.param(:frequency) * QD.op(:x) / 2,
    )

    @test haskey(component.parameters, :frequency)
    @test :x in component.operators
    @test component.hamiltonian isa QD.OperatorExpr
    @test !hasproperty(component, :hamiltonian_expr)

    expr = QD.absolute_path(component.hamiltonian, component)
    @test expr isa QD.OperatorExpr
    @test expr.args[1].args == ((:components, :q1, :params, :frequency),)
    @test expr.args[2].args[1].path == (:components, :q1, :operators, :x)

    @test norm(QD.numerical(component, :x; dimension = 2) - QD.numerical_operator(component.spec, :x; dimension = 2)) ≈ 0
    @test norm(QD.numerical(component, :x) - QD.numerical_operator(component.spec, :x; dimension = 2)) ≈ 0
    @test norm(QD.numerical(component; dimension = 2) - 2.5 * qt.sigmax()) ≈ 0
    @test norm(QD.numerical(component) - 2.5 * qt.sigmax()) ≈ 0
    @test norm(QD.numerical(component, component.hamiltonian; dimension = 2) - 2.5 * qt.sigmax()) ≈ 0
    @test norm(QD.numerical(component, component.hamiltonian) - 2.5 * qt.sigmax()) ≈ 0
    @test norm(QD.numerical(component, component.hamiltonian; dimension = 2, frequency = 10.0) - 5.0 * qt.sigmax()) ≈ 0
    @test norm(QD.numerical(component, component.hamiltonian; dimension = 2, params = Dict(:frequency => 10.0)) - 5.0 * qt.sigmax()) ≈ 0

    probe = QD.Component(
        DimensionProbeSpec(),
        :probe;
        hamiltonian = QD.op(:a) + QD.op(:b),
    )

    @test QD.numerical(probe; dimension = 7) == fill(21, 1, 1)
    @test QD.numerical(probe) == fill(9, 1, 1)
    @test_throws Exception QD.numerical(component, :bad; dimension = 2)
    @test_throws Exception QD.numerical(component, QD.op(:q2, :x); dimension = 2)
end
