@testset "charge-basis transmon components" begin
    dimension = 15
    EC, EJ, ng = 0.27, 18.4, 0.23
    transmon = make_transmon(:t, EC, EJ; ng, dimension, metadata=(source=:test,))
    @test transmon.hamiltonian isa Sym
    @test all(g -> g isa OperatorGenerator, values(transmon.operator_generators))
    @test transmon.label == "t"
    @test transmon.metadata == (source=:test,)
    @test make_transmon(:t; EC, EJ, label=nothing).label === nothing
    @test make_transmon(; name=:t, EC, EJ, label=:short).label == "short"
    @test transmon.operator_generators["charge"]() isa QT.QuantumObject
    @test size(transmon.operator_generators["jump"](7; ignored=true)) == (7, 7)
    @test make_transmon(:default, EC, EJ).operator_generators["charge"].default_dimension == 101
    for bad in (0, -1, 2, 10)
        @test_throws ArgumentError make_transmon(:bad, EC, EJ; dimension=bad)
        @test_throws ArgumentError transmon.operator_generators["charge"](bad)
        @test_throws ArgumentError transmon.operator_generators["jump"](bad)
    end

    charge = QT.num(dimension) - (dimension - 1) ÷ 2
    jump = QT.tunneling(dimension, 1)
    reference = 4 * EC * (charge - ng)^2 - EJ * jump / 2
    device = QuantumDevice(:transmons)
    add_component!(device, transmon)
    idle = build(device, Leaf(component=:t))
    @test idle.spectrum ≈ real.(QT.eigenenergies(reference))
    @test input_basis_map(idle)' * reference * input_basis_map(idle) ≈ QuantumDevices.spectrum_operator(idle.spectrum)
    for (name, operator) in ((:charge, charge), (:jump, jump))
        @test idle.operators[OperatorRef(:t, name)] ≈ operator
    end
    changed = build(device, Leaf(component=:t); parameters=Dict(transmon.variables["EJ"]=>14.0))
    @test changed.spectrum ≈ real.(QT.eigenenergies(4 * EC * (charge - ng)^2 - 14 * jump / 2))
    @test transmon.parameters["EJ"] == EJ
    @test build(device, Leaf(component=:t)).spectrum ≈ idle.spectrum

    tunable = make_flux_tunable_transmon(:f, EC, EJ; ng, d=-0.17, dimension)
    add_component!(device, tunable)
    @test Set(keys(tunable.variables)) == Set(["EC", "EJmax", "d", "ng", "flux"])
    @test tunable.label == "f"
    @test make_flux_tunable_transmon(; name=:g, EC, EJmax=EJ, label=nothing).label === nothing
    @test make_flux_tunable_transmon(:g; EC, EJmax=EJ, label=:short).label == "short"
    @test build(device, Leaf(component=:f)).spectrum ≈ idle.spectrum
    for flux in (0.5, 0.28, 1.28), asymmetry in (-0.17, 0.0)
        effective_EJ = EJ * sqrt(cospi(flux)^2 + asymmetry^2 * sinpi(flux)^2)
        model = build(device, Leaf(component=:f); parameters=Dict(
            tunable.variables["flux"]=>flux, tunable.variables["d"]=>asymmetry))
        expected = 4 * EC * (charge - ng)^2 - effective_EJ * jump / 2
        @test model.spectrum ≈ real.(QT.eigenenergies(expected))
        @test model.operators[OperatorRef(:f, :charge)] ≈ charge
    end
    @test tunable.parameters["flux"] == 0.0
    @test tunable.parameters["d"] == -0.17
    @test_throws ArgumentError make_flux_tunable_transmon(:bad, EC, EJ; dimension=8)

    # Scalar function coefficients travel through the interaction term compiler.
    p = Symbolics.variable(:angle)
    coefficient = sqrt(cospi((p))^2 + sinpi((p))^2)
    @test evaluate(coefficient, Dict(p=>0.31), Dict()) ≈ 1
    @test occursin("sqrt", sprint(show, coefficient))
    @test_throws ArgumentError cospi(op(:t, :charge))
    @test_throws ArgumentError sinpi(op(:t, :charge))
    @test_throws ArgumentError sqrt(op(:t, :charge))
    interaction = Interaction(:capacitive, 0.03 * coefficient * op(:t, :charge) * op(:f, :charge); parameters=Dict(p=>0.31))
    add_interaction!(device, interaction)
    add_recipe!(device, :t_leaf, Leaf(component=:t))
    add_recipe!(device, :f_leaf, Leaf(component=:f))
    coupled = build(device, Combine(inputs=[Input(:t, :t_leaf, [0, 1, 2]),
        Input(:f, :f_leaf, [0, 1, 2])], interactions=[:capacitive]))
    embedding = QuantumDevices.selection(dimension, [1, 2, 3])
    projected_charge = embedding' * input_basis_map(idle)' * idle.operators[OperatorRef(:t, :charge)] * input_basis_map(idle) * embedding
    energies = QuantumDevices.spectrum_operator(idle.spectrum[1:3])
    expected = QT.tensor(energies, QT.qeye(3)) + QT.tensor(QT.qeye(3), energies) +
        0.03 * QT.tensor(projected_charge, projected_charge)
    @test coupled.spectrum ≈ real.(QT.eigenenergies(expected))
    roundoff = QT.sigmaz() + 1e-14 * QT.destroy(2)
    energies_roundoff, _ = QuantumDevices.diagonalize(roundoff)
    @test energies_roundoff ≈ real.(QT.eigenenergies((roundoff + roundoff') / 2))
    @test_throws ArgumentError QuantumDevices.diagonalize(QT.sigmaz() + 1e-3 * QT.destroy(2))
end
