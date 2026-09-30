@testset "Model construction" begin
    a, b = make_qubit("a", 1.0), make_qubit("b", 1.7)
    interaction = param(:g) * op(:a_x) * op(:b_x)
    model = make_model([a, b], interaction, (g = 0.03,))
    direct = 0.5 * tensor(sigmaz(), qeye(2)) + 0.85 * tensor(qeye(2), sigmaz()) + 0.03 * tensor(sigmax(), sigmax())
    @test model.H ≈ direct
    @test model.coupling_parameters == (g = 0.03,)
    @test model.max_dimension == 10^4
    @test model.states === model.eigensystem.states
    @test :confidence in propertynames(model)
    @test_throws ArgumentError model.nonexistent

    model2 = setpath(model, "components/1/parameters/ν", 1.2)
    fresh_model2 = make_model([make_qubit("a", 1.2), b], interaction, (g = 0.03,))
    @test model2.H ≈ fresh_model2.H
    @test model2.parameters == fresh_model2.parameters
    for key in keys(model2.states)
        @test abs2(dot(model2.states[key], fresh_model2.states[key])) ≈ 1 atol = 1e-10
    end
    @test model.components[1].parameters.ν == 1.0

    stronger = setpath(model, "coupling_parameters/g", 0.05)
    @test stronger.coupling_parameters.g == 0.05
    @test stronger.H ≈ make_model([a, b], interaction, (g = 0.05,)).H
    replacement_interaction = 2 * param(:g) * op(:a_x) * op(:b_x)
    replaced = setpath(model, "interactions", replacement_interaction)
    @test replaced.H ≈ make_model([a, b], replacement_interaction, (g = 0.03,)).H

    gate = DeviceGate((drive = 0.1,), param(:drive) * op(:a_x), 1.0)
    with_gate = setpath(model, "gates/test", gate)
    @test !haskey(model.gates, :test)
    @test with_gate.gates[:test] === gate
    @test with_gate.H === model.H
    @test with_gate.eigensystem === model.eigensystem
    rebuilt_with_gate = setpath(with_gate, "components/1/parameters/ν", 1.2)
    @test rebuilt_with_gate.gates[:test] === gate
    @test rebuilt_with_gate.gates !== with_gate.gates

    @test_throws ArgumentError setpath(model, "parameters", model.parameters)
    @test_throws ArgumentError setpath(model, "H", model.H)
    @test_throws ArgumentError setpath(model, "states", model.states)
    @test_throws ArgumentError setpath(model, "components/1/name", "renamed")
    @test_throws ArgumentError setpath(model, "components", [a])
    @test_throws ArgumentError setpath(model, "components", [b, a])
    @test_throws ArgumentError setpath(model, "max_dimension", 3)
    @test_throws ArgumentError make_model([a, a], interaction, (g = 0.03,))
    @test_throws ArgumentError make_model([a, b], interaction, (g = 0.03,); max_dimension = 3)
    @test_throws ArgumentError make_model([a, b], interaction, (g = 0.03,); truncation_dimensions = Dict(a => 3))
    @test_throws ArgumentError make_model([a, b], interaction, (g = 0.03,); truncation_dimensions = Dict(a => 0))
    @test_throws ArgumentError make_model([a, b], interaction, (g = 0.03, a_ν = 1.0))
    @test_throws ArgumentError make_model([a, b], interaction, (;))
    @test_throws ArgumentError make_model([a], op(:missing), (;))
    c = make_transmon("c", 0.2, 5.0, 9)
    Hc = numerical(c.hamiltonian, c.operators, c.parameters)
    P = QD.projection(9, 1:3)
    reduced = QD.truncate(Hc, P)
    mix = make_model([c, b], 0 * op(:b_x), (;); truncation_dimensions = Dict(c => 3))
    @test mix.H ≈ tensor(reduced, qeye(2)) + tensor(qeye(3), 0.85 * sigmaz())
    remapped = setpath(mix, "components/1/parameters/EJ", 6.0)
    @test get(remapped.truncation_dimensions, remapped.components[1], nothing) == 3
    @test !haskey(remapped.truncation_dimensions, c)
    fresh_c = make_transmon("c", 0.2, 6.0, 9)
    fresh_remapped = make_model(
        [fresh_c, b],
        0 * op(:b_x),
        (;);
        truncation_dimensions = Dict(fresh_c => 3),
    )
    @test remapped.H ≈ fresh_remapped.H
    full = make_model([c], 0 * op(:c_charge), (;); truncation_dimensions = Dict(c => 3))
    @test full.H ≈ reduced
    naive = numerical(full.hamiltonian, full.operators, full.parameters)
    @test norm(full.H - naive) > 1e-4
    word_model = make_model([c, b], 0.07 * op(:c_charge) * op(:b_x) * op(:c_charge), (;);
        truncation_dimensions = Dict(c => 3))
    @test word_model.H ≈ mix.H + 0.07 * tensor(QD.truncate(c.operators.charge^2, P), sigmax())
end
