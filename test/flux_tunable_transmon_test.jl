using Test
import LinearAlgebra: diag, norm
import QuantumDevices as QD

@testset "flux-tunable transmon component" begin
    spec = QD.FluxTunableTransmonSpec(0.2, 12.0, 8.0; flux = 0.1, ng = 0.05)
    component = QD.Component(spec, :ft)
    ej_eff =
        (
            ((QD.param(spec.EJ1) + QD.param(spec.EJ2)) * cos(pi * QD.param(spec.flux)))^2 +
            ((QD.param(spec.EJ1) - QD.param(spec.EJ2)) * sin(pi * QD.param(spec.flux)))^2
        )^(1 / 2)
    default_hamiltonian =
        4 * QD.param(spec.EC) * (QD.op(:n) - QD.param(spec.ng) * QD.op(:identity))^2 -
        ej_eff * QD.op(:tunneling) / 2

    @test QD.parameters(component)[QD.ParamPath(:EC)].default == 0.2
    @test QD.parameters(component)[QD.ParamPath(:EJ1)].default == 12.0
    @test QD.parameters(component)[QD.ParamPath(:EJ2)].default == 8.0
    @test QD.parameters(component)[QD.ParamPath(:flux)].default == 0.1
    @test QD.parameters(component)[QD.ParamPath(:ng)].default == 0.05
    @test spec.dimension == QD.Dimension(Inf)
    @test !haskey(QD.parameters(component), QD.ParamPath(:dimension))
    @test QD.operators(spec) == Set([:identity, :n, :tunneling])
    @test QD.operators(component) == Set{Any}([:identity, :n, :tunneling])
    @test !(:charge in QD.operators(component))
    @test !(:hamiltonian in QD.operators(component))
    @test QD.hamiltonian(component) == default_hamiltonian

    hamiltonian = QD.numerical(component; dimension = (9,))
    charge = QD.numerical(component, :n; dimension = (9,))
    @test size(hamiltonian) == (9, 9)
    @test norm(hamiltonian - adjoint(hamiltonian)) < 1e-12
    @test real.(diag(charge.data)) == collect(-4:4)
    @test real.(diag(hamiltonian.data)) ≈ 4 * 0.2 .* (collect(-4:4) .- 0.05).^2

    zero_flux = QD.Component(
        QD.FluxTunableTransmonSpec(0.2, 12.0, 8.0; flux = 0.0, ng = 0.05),
        :zero_flux,
    )
    fixed = QD.Component(QD.TransmonSpec(0.2, 20.0; ng = 0.05), :fixed)
    @test norm(
        QD.numerical(zero_flux; dimension = (9,)) -
        QD.numerical(fixed; dimension = (9,)),
    ) ≈ 0

    symmetric_half_flux = QD.Component(
        QD.FluxTunableTransmonSpec(0.2, 10.0, 10.0; flux = 0.5),
        :symmetric,
    )
    symmetric_hamiltonian = QD.numerical(symmetric_half_flux; dimension = (9,))
    @test abs(symmetric_hamiltonian.data[1, 2]) < 1e-12
    @test abs(symmetric_hamiltonian.data[2, 1]) < 1e-12

    shifted = QD.numerical(component; dimension = (9,), flux = 0.25)
    shifted_from_params = QD.numerical(
        component;
        dimension = (9,),
        params = Dict(:flux => 0.25),
    )
    @test norm(shifted - shifted_from_params) ≈ 0
    @test norm(shifted - hamiltonian) > 1e-8
    @test norm(QD.numerical(component; dimension = (9,), EC = 0.25) - hamiltonian) > 1e-8
    @test norm(QD.numerical(component; dimension = (9,), EJ1 = 13.0) - hamiltonian) > 1e-8
    @test norm(QD.numerical(component; dimension = (9,), EJ2 = 9.0) - hamiltonian) > 1e-8
    @test norm(QD.numerical(component; dimension = (9,), ng = 0.1) - hamiltonian) > 1e-8

    custom = QD.Component(
        QD.FluxTunableTransmonSpec(0.2, 12.0, 8.0; hamiltonian = QD.op(:n)),
        :custom,
    )
    @test QD.hamiltonian(custom) == QD.op(:n)
    @test norm(QD.numerical(custom; dimension = (9,)) - charge) ≈ 0

    extra = QD.Component(
        QD.FluxTunableTransmonSpec(0.2, 12.0, 8.0; flux = 0.1, ng = 0.05, extra_terms = QD.op(:n)),
        :with_extra,
    )
    extra_ej_eff =
        (
            ((QD.param(extra.spec.EJ1) + QD.param(extra.spec.EJ2)) *
                cos(pi * QD.param(extra.spec.flux)))^2 +
            ((QD.param(extra.spec.EJ1) - QD.param(extra.spec.EJ2)) *
                sin(pi * QD.param(extra.spec.flux)))^2
        )^(1 / 2)
    extra_default =
        4 * QD.param(extra.spec.EC) *
            (QD.op(:n) - QD.param(extra.spec.ng) * QD.op(:identity))^2 -
        extra_ej_eff * QD.op(:tunneling) / 2
    @test extra.spec.hamiltonian == extra_default + QD.op(:n)
    @test norm(QD.numerical(extra; dimension = (9,)) - (hamiltonian + charge)) ≈ 0

    child = QD.ModelSpec(
        :ft_child,
        [component],
        Dict(:ft => 4);
        pretruncation_dims = Dict(:ft => 9),
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    child_model = QD.model(child)
    @test size(child_model.hamiltonian) == (4, 4)
    @test child_model.spec === child

    compact = sprint(show, spec)
    plain = sprint(show, MIME"text/plain"(), spec)
    @test compact == "FluxTunableTransmonSpec(5 parameters, 3 operators)"
    @test occursin("flux = 0.1", plain)
    @test occursin("Hamiltonian:", plain)

    @test_throws Exception QD.numerical(component)
    @test_throws Exception QD.numerical(component, :bad; dimension = (9,))
    @test_throws Exception QD.numerical(component; dimension = (8,))
    @test_throws Exception QD.FluxTunableTransmonSpec(0.0, 12.0, 8.0)
    @test_throws Exception QD.FluxTunableTransmonSpec(0.2, 0.0, 8.0)
    @test_throws Exception QD.FluxTunableTransmonSpec(0.2, 12.0, 0.0)
    @test_throws Exception QD.ModelSpec(
        :even_flux_tunable,
        [component],
        Dict(:ft => 3);
        pretruncation_dims = Dict(:ft => 8),
    )
end
