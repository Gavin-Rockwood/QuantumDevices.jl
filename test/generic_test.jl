using Test
import LinearAlgebra: diag, norm
import QuantumDevices as QD

@testset "generic component" begin
    spec = QD.GenericSpec([0.0, 1.0, 2.0])
    component = QD.Component(spec, :generic)

    @test spec.operators isa NamedTuple
    @test QD.available_operators(spec) == Set([:identity, :energy])
    @test component.operators == Set{Any}([:identity, :energy])
    @test !(:hamiltonian in component.operators)
    @test component.hamiltonian == QD.op(:energy)

    hamiltonian = QD.numerical(component)
    identity = QD.numerical(component, :identity)
    @test size(hamiltonian) == (3, 3)
    @test size(identity) == (3, 3)
    @test real.(diag(hamiltonian.data)) == [0.0, 1.0, 2.0]

    shifted = QD.Component(
        QD.GenericSpec([0.0, 1.0, 2.0]; extra_terms = QD.op(:identity)),
        :shifted,
    )
    @test shifted.spec.hamiltonian == QD.op(:energy) + QD.op(:identity)
    @test norm(QD.numerical(shifted) - (hamiltonian + identity)) ≈ 0

    custom = QD.Component(
        QD.GenericSpec([0.0, 1.0, 2.0]; hamiltonian = QD.op(:identity)),
        :custom,
    )
    @test custom.hamiltonian == QD.op(:identity)
    @test norm(QD.numerical(custom) - identity) ≈ 0
end
