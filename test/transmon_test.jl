using Test
import LinearAlgebra: diag, norm
import QuantumDevices as QD

@testset "transmon component" begin
    spec = QD.TransmonSpec(0.2, 20.0; ng = 0.1)
    component = QD.Component(spec, :transmon)
    default_hamiltonian =
        4 * QD.param(spec.EC) * (QD.op(:n) - QD.param(spec.ng) * QD.op(:identity))^2 -
        QD.param(spec.EJ) * QD.op(:tunneling) / 2

    @test QD.parameters(component)[QD.ParamPath(:EC)].default == 0.2
    @test QD.parameters(component)[QD.ParamPath(:EJ)].default == 20.0
    @test QD.parameters(component)[QD.ParamPath(:ng)].default == 0.1
    @test spec.dimension == QD.Dimension(Inf)
    @test !haskey(QD.parameters(component), QD.ParamPath(:dimension))
    @test !haskey(QD.parameters(component), QD.ParamPath(:n_cutoff))
    @test spec.operators isa NamedTuple
    @test QD.operators(spec) == Set([:identity, :n, :tunneling])
    @test QD.operators(component) == Set{Any}([:identity, :n, :tunneling])
    @test !(:charge in QD.operators(component))
    @test !(:hamiltonian in QD.operators(component))
    @test QD.hamiltonian(component) == default_hamiltonian

    hamiltonian = QD.numerical(component; dimension = (9,))
    charge = QD.numerical(component, :n; dimension = (9,))
    @test size(hamiltonian) == (9, 9)
    @test size(charge) == (9, 9)
    @test norm(hamiltonian - adjoint(hamiltonian)) < 1e-12
    @test norm(charge - adjoint(charge)) < 1e-12
    @test real.(diag(charge.data)) == collect(-4:4)
    @test real.(diag(hamiltonian.data)) ≈ 4 * 0.2 .* (collect(-4:4) .- 0.1).^2
    @test hamiltonian.data[1, 2] ≈ -20.0 / 2
    @test hamiltonian.data[2, 1] ≈ -20.0 / 2

    shifted = QD.numerical(component; dimension = (9,), ng = 0.2)
    shifted_from_params = QD.numerical(
        component;
        dimension = (9,),
        params = Dict(:ng => 0.2),
    )
    @test norm(shifted - shifted_from_params) ≈ 0
    @test norm(shifted - hamiltonian) > 1e-8
    @test norm(QD.numerical(component; dimension = (9,), EC = 0.25) - hamiltonian) > 1e-8

    custom = QD.Component(
        QD.TransmonSpec(0.2, 20.0; hamiltonian = QD.op(:n)),
        :custom,
    )
    @test QD.hamiltonian(custom) == QD.op(:n)
    @test norm(QD.numerical(custom; dimension = (9,)) - charge) ≈ 0

    extra = QD.Component(
        QD.TransmonSpec(0.2, 20.0; ng = 0.1, extra_terms = QD.op(:n)),
        :with_extra,
    )
    extra_default =
        4 * QD.param(extra.spec.EC) *
            (QD.op(:n) - QD.param(extra.spec.ng) * QD.op(:identity))^2 -
        QD.param(extra.spec.EJ) * QD.op(:tunneling) / 2
    @test extra.spec.hamiltonian == extra_default + QD.op(:n)
    @test norm(QD.numerical(extra; dimension = (9,)) - (hamiltonian + charge)) ≈ 0

    @test_throws Exception QD.numerical(component)
    @test_throws Exception QD.numerical(component, :bad; dimension = (9,))
    @test_throws Exception QD.numerical(component, :charge; dimension = (9,))
    @test_throws Exception QD.numerical(component; dimension = (8,))
    @test_throws Exception QD.TransmonSpec(0.0, 20.0)
    @test_throws Exception QD.TransmonSpec(0.2, 0.0)

    child = QD.ModelSpec(
        :t1,
        [component],
        Dict(:transmon => 4);
        pretruncation_dims = Dict(:transmon => 9),
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    child_model = QD.model(child)
    @test size(child_model.hamiltonian) == (4, 4)
    @test length(child_model.states) == 4
    @test length(child_model.energies) == 4
    @test child_model.spec === child
    effective = QD.component(child_model)
    @test (:transmon, :n) in QD.operators(effective)
    @test size(QD.numerical(effective, (:transmon, :n))) == (4, 4)

    compact = sprint(show, spec)
    plain = sprint(show, MIME"text/plain"(), spec)
    @test compact == "TransmonSpec(3 parameters, 3 operators)"
    @test occursin("EC = 0.2", plain)
    @test occursin("Dimension: (Inf,)", plain)
    @test !occursin("n_cutoff", plain)
    @test occursin("Hamiltonian: 4 × EC × (n - ng × identity)² - EJ × tunneling / 2", plain)
    @test !occursin("Matrix", plain)
end
