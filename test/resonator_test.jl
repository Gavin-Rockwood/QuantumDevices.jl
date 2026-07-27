using Test
import LinearAlgebra: diag, norm
import QuantumDevices as QD
import QuantumToolbox as qt

@testset "resonator component" begin
    spec = QD.ResonatorSpec(6.25; dimension = 8)
    resonator = QD.Component(spec, :cavity)

    @test QD.parameters(resonator)[QD.ParamPath(:frequency)].default == 6.25
    @test spec.dimension == QD.Dimension(8)
    @test !haskey(QD.parameters(resonator), QD.ParamPath(:dimension))
    @test QD.hamiltonian(resonator) == QD.param(spec.frequency) * QD.op(:n)
    @test spec.operators isa NamedTuple
    @test QD.operators(resonator) == Set{Any}([:a, :adag, :n, :number, :q, :p])

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

    @test size(QD.numerical(resonator; dimension = (4,))) == (4, 4)
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
    @test built.name == :mode1
    @test size(QD.numerical(built)) == (5, 5)

    shift = QD.DeviceParameter(:shift; default = 0.06)
    extra = QD.param(shift) * QD.op(:number)
    with_extra = QD.Component(QD.ResonatorSpec(
        6.0;
        dimension = 4,
        extra_terms = extra,
    ), :mode)
    @test with_extra.spec.hamiltonian == QD.param(with_extra.spec.frequency) * QD.op(:n) + extra
    @test QD.parameters(with_extra)[QD.ParamPath(:shift)] === shift
    @test real.(diag(QD.numerical(with_extra).data)) ≈
          6.06 .* collect(0:3)

    compact = sprint(show, spec)
    plain = sprint(show, MIME"text/plain"(), spec)
    @test compact == "ResonatorSpec(1 parameter, 6 operators)"
    @test occursin("Dimension: (8,)", plain)
    @test occursin("frequency = 6.25", plain)
    @test occursin("Hamiltonian: frequency × n", plain)
    @test !occursin("var\"#", plain)
    @test !occursin("Matrix", plain)
end
