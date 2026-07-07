using Test
import LinearAlgebra: diag, norm
import QuantumDevices as QD

@testset "flux-tunable transmon component" begin
    spec = QD.FluxTunableTransmonSpec(0.2, 12.0, 8.0; flux = 0.1, ng = 0.05)
    component = QD.Component(spec, :ft)
    ej_eff =
        (
            ((QD.param(:EJ1) + QD.param(:EJ2)) * cos(pi * QD.param(:flux)))^2 +
            ((QD.param(:EJ1) - QD.param(:EJ2)) * sin(pi * QD.param(:flux)))^2
        )^(1 / 2)
    default_hamiltonian =
        4 * QD.param(:EC) * (QD.op(:n) - QD.param(:ng) * QD.op(:identity))^2 -
        ej_eff * QD.op(:tunneling) / 2

    @test component.parameters[:EC].default == 0.2
    @test component.parameters[:EJ1].default == 12.0
    @test component.parameters[:EJ2].default == 8.0
    @test component.parameters[:flux].default == 0.1
    @test component.parameters[:ng].default == 0.05
    @test component.parameters[:dimension].default == Inf
    @test QD.available_operators(spec) == Set([:identity, :n, :tunneling])
    @test component.operators == Set{Any}([:identity, :n, :tunneling])
    @test !(:charge in component.operators)
    @test !(:hamiltonian in component.operators)
    @test component.hamiltonian == default_hamiltonian

    hamiltonian = QD.numerical(component; dimension = 9)
    charge = QD.numerical(component, :n; dimension = 9)
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
        QD.numerical(zero_flux; dimension = 9) -
        QD.numerical(fixed; dimension = 9),
    ) ≈ 0

    symmetric_half_flux = QD.Component(
        QD.FluxTunableTransmonSpec(0.2, 10.0, 10.0; flux = 0.5),
        :symmetric,
    )
    symmetric_hamiltonian = QD.numerical(symmetric_half_flux; dimension = 9)
    @test abs(symmetric_hamiltonian.data[1, 2]) < 1e-12
    @test abs(symmetric_hamiltonian.data[2, 1]) < 1e-12

    shifted = QD.numerical(component; dimension = 9, flux = 0.25)
    shifted_from_params = QD.numerical(
        component;
        dimension = 9,
        params = Dict(:flux => 0.25),
    )
    @test norm(shifted - shifted_from_params) ≈ 0
    @test norm(shifted - hamiltonian) > 1e-8
    @test norm(QD.numerical(component; dimension = 9, EC = 0.25) - hamiltonian) > 1e-8
    @test norm(QD.numerical(component; dimension = 9, EJ1 = 13.0) - hamiltonian) > 1e-8
    @test norm(QD.numerical(component; dimension = 9, EJ2 = 9.0) - hamiltonian) > 1e-8
    @test norm(QD.numerical(component; dimension = 9, ng = 0.1) - hamiltonian) > 1e-8

    custom = QD.Component(
        QD.FluxTunableTransmonSpec(0.2, 12.0, 8.0; hamiltonian = QD.op(:n)),
        :custom,
    )
    @test custom.hamiltonian == QD.op(:n)
    @test norm(QD.numerical(custom; dimension = 9) - charge) ≈ 0

    extra = QD.Component(
        QD.FluxTunableTransmonSpec(0.2, 12.0, 8.0; flux = 0.1, ng = 0.05, extra_terms = QD.op(:n)),
        :with_extra,
    )
    @test extra.spec.hamiltonian == default_hamiltonian + QD.op(:n)
    @test norm(QD.numerical(extra; dimension = 9) - (hamiltonian + charge)) ≈ 0

    child = QD.ModelSpec(
        :ft_child;
        components = Dict(:ft => component),
        dims = Dict(:ft => 4),
        initialization_dims = Dict(:ft => 9),
    )
    child_model = QD.model(child)
    @test size(child_model.hamiltonian) == (4, 4)
    @test haskey(child_model.spec.children, :ft)

    compact = sprint(show, spec)
    plain = sprint(show, MIME"text/plain"(), spec)
    @test compact == "FluxTunableTransmonSpec(6 parameters, 3 operators)"
    @test occursin("flux = 0.1", plain)
    @test occursin("Hamiltonian:", plain)

    @test_throws Exception QD.numerical(component)
    @test_throws Exception QD.numerical(component, :bad; dimension = 9)
    @test_throws Exception QD.numerical(component; dimension = 8)
    @test_throws Exception QD.FluxTunableTransmonSpec(0.0, 12.0, 8.0)
    @test_throws Exception QD.FluxTunableTransmonSpec(0.2, 0.0, 8.0)
    @test_throws Exception QD.FluxTunableTransmonSpec(0.2, 12.0, 0.0)
    @test_throws Exception QD.model(QD.ModelSpec(
        :even_flux_tunable;
        components = Dict(:ft => component),
        dims = Dict(:ft => 3),
        initialization_dims = Dict(:ft => 8),
    ))
end
