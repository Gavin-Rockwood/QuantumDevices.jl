using Test
import LinearAlgebra: diag, norm
import QuantumDevices as QD
import QuantumToolbox as qt

@testset "resonator component" begin
    spec = QD.ResonatorSpec(6.25; dimension = 8)
    resonator = QD.Component(spec, :cavity)

    @test resonator.parameters[:frequency].default == 6.25
    @test resonator.parameters[:dimension].default == 8
    @test resonator.hamiltonian == QD.param(:frequency) * QD.op(:n)
    @test spec.operators isa NamedTuple
    @test resonator.operators == Set{Any}([:a, :adag, :n, :number, :q, :p])

    hamiltonian = QD.numerical(resonator)
    annihilation = QD.numerical(resonator, :a)
    creation = QD.numerical(resonator, :adag)
    number = QD.numerical(resonator, :n)
    @test size(hamiltonian) == (8, 8)
    @test real.(diag(hamiltonian.data)) ≈ 6.25 .* collect(0:7)
    @test norm(creation - adjoint(annihilation)) ≈ 0
    @test norm(number - QD.numerical(resonator, :number)) ≈ 0
    @test norm(QD.numerical(resonator, :q) - adjoint(QD.numerical(resonator, :q))) ≈ 0
    @test norm(QD.numerical(resonator, :p) - adjoint(QD.numerical(resonator, :p))) ≈ 0

    @test size(QD.numerical(resonator; dimension = 4)) == (4, 4)
    @test real.(diag(QD.numerical(resonator; frequency = 7.0).data)) ≈
          7.0 .* collect(0:7)
    @test real.(diag(QD.numerical(
        resonator;
        params = Dict(:frequency => 7.0),
    ).data)) ≈ 7.0 .* collect(0:7)

    @test_throws Exception QD.numerical(resonator, :bad)
    @test_throws Exception QD.ResonatorSpec(0.0)
    @test_throws Exception QD.ResonatorSpec(6.0; dimension = 0)

    built = QD.Component(QD.ResonatorSpec(6.0; dimension = 5), :mode1)
    @test built.id == :mode1
    @test size(QD.numerical(built)) == (5, 5)

    extra = QD.param(:frequency) * QD.op(:number) / 100
    with_extra = QD.Component(QD.ResonatorSpec(
        6.0;
        dimension = 4,
        extra_terms = extra,
    ), :mode)
    @test with_extra.spec.hamiltonian == QD.param(:frequency) * QD.op(:n) + extra
    @test real.(diag(QD.numerical(with_extra).data)) ≈
          6.06 .* collect(0:3)

    compact = sprint(show, spec)
    plain = sprint(show, MIME"text/plain"(), spec)
    @test compact == "ResonatorSpec(2 parameters, 6 operators)"
    @test occursin("frequency = 6.25", plain)
    @test occursin("Hamiltonian: frequency × n", plain)
    @test !occursin("var\"#", plain)
    @test !occursin("Matrix", plain)
end
