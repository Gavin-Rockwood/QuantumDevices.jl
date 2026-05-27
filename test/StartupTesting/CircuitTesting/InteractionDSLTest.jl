@testset "Typed Interaction DSL" begin
    transmon = QD.Circuits.init_transmon(0.2, 1.0, 3; name = "transmon", n_full = 12)
    resonator = QD.Circuits.init_resonator(5.0, 3; name = "resonator")
    components = [transmon, resonator]

    @test norm(QD.Circuits.get_operator(transmon, :n) - transmon.n_op) < 1e-12
    @test norm(QD.Circuits.get_operator(resonator, :a) - resonator.a_op) < 1e-12

    interactions = [QD.Circuits.coupling(0.1, QD.Circuits.op(:transmon, :n) * QD.Circuits.op(:resonator, :a))]
    circuit = QD.Circuits.init_circuit(
        components,
        interactions;
        operators_to_add = Dict{String, Any}("na" => QD.Circuits.op(:transmon, :n) * QD.Circuits.op(:resonator, :a)),
        dressed_kwargs = Dict(:step_number => 3),
    )
    expected = qt.tensor(transmon.n_op, resonator.a_op)
    @test norm(circuit.ops["na"] - expected) < 1e-12

    expr = QD.Circuits.op(:transmon, :n) * QD.Circuits.op(:resonator, :a)
    circuit_hc = QD.Circuits.init_circuit(
        components,
        [QD.Circuits.coupling(0.1, expr; hc = true)];
        dressed_kwargs = Dict(:step_number => 3),
    )
    circuit_explicit = QD.Circuits.init_circuit(
        components,
        [QD.Circuits.interaction(0.1 * expr + (0.1 * expr)')];
        dressed_kwargs = Dict(:step_number => 3),
    )
    @test norm(circuit_hc.H_op - circuit_explicit.H_op) < 1e-12

    @test_throws ErrorException QD.Circuits.get_operator(transmon, :a)
    @test_throws ErrorException QD.Circuits.init_circuit(
        components,
        [QD.Circuits.coupling(0.1, QD.Circuits.op(:missing, :n))],
    )
end
