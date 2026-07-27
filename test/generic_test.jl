using Test
import LinearAlgebra: diag, norm
import QuantumDevices as QD

@testset "generic component" begin
    spec = QD.GenericSpec([0.0, 1.0, 2.0])
    component = QD.Component(spec, :generic)

    @test spec.operators isa Dict
    @test QD.operators(spec) == Set([:identity, :energy])
    @test QD.operators(component) == Set{Any}([:identity, :energy])
    @test !(:hamiltonian in QD.operators(component))
    @test QD.hamiltonian(component) == QD.op(:energy)
    @test isempty(QD.parameters(component))
    @test spec.dimension == QD.Dimension(3)
    @test !haskey(QD.parameters(component), QD.ParamPath(:dimension))
    @test spec.spectrum == [0.0, 1.0, 2.0]

    hamiltonian = QD.numerical(component)
    identity = QD.numerical(component, :identity)
    @test size(hamiltonian) == (3, 3)
    @test size(identity) == (3, 3)
    @test real.(diag(hamiltonian.data)) == [0.0, 1.0, 2.0]

    projected = ComplexF64[0 1; 1 0]
    projected_spec = QD.GenericSpec([0.0, 1.0];
        operators = Dict(:x => (; dimension, kwargs...) -> begin
            keep = only(dimension)
            projected[1:keep, 1:keep]
        end),
        source = spec)
    @test projected_spec.source === spec
    @test QD.numerical(QD.Component(projected_spec, :projected), :x; dimension = (2,)) == projected

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
    @test QD.hamiltonian(custom) == QD.op(:identity)
    @test norm(QD.numerical(custom) - identity) ≈ 0
end
