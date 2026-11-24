@testset "gate recipes" begin
    device = QuantumDevice(:gate_recipes)
    transmon = add_component!(device, make_flux_tunable_transmon(:f, 0.25, 18; d=0.15, dimension=11))
    flux = transmon.variables["flux"]
    recipe = Recipe(device, [transmon]; states=Dict(transmon=>0:3), parameters=Dict(flux=>0.2))
    pulse = Pulse(GaussianZero(1;duration=4); duration=4, amplitude=0.1, offset=0.2, delay=1)
    gate = GateRecipe(recipe; controls=[Control(:flux,flux,pulse)])
    prepared = build(device,gate)
    @test gate.duration == 5
    @test prepared.model.spectrum ≈ build(device,recipe).spectrum
    expected = make_hamiltonian(prepared.model,[flux])
    actual = make_hamiltonian(prepared)
    for t in (0.,1.,2.,4.,5.,7.)
        @test actual(nothing,t) ≈ expected(Dict(flux=>pulse),t)
    end
    @test check_endpoints(prepared;atol=1e-14)
    @test pulse_tstops(prepared) == [0,1,5]
    @test_throws ArgumentError GateRecipe(recipe)
    @test_throws ArgumentError GateRecipe(recipe;controls=[Control(:flux,flux,pulse)],duration=4)
    @test_throws ArgumentError GateRecipe(recipe;controls=[Control(:flux,flux,pulse),Control(:other,flux,pulse)])
    @test_throws ArgumentError GateRecipe(recipe;controls=[Control(:same,flux,pulse),Control(:same,Symbolics.variable(:unknown),pulse)])
    @test_throws ArgumentError build(device,GateRecipe(recipe;controls=[Control(:unknown,Symbolics.variable(:unknown),pulse)]))
    @test make_hamiltonian(build(device,GateRecipe(recipe;duration=3))) ≈ make_hamiltonian(build(device,recipe))

    amplitude = Symbolics.variable(:drive)
    n = recipe.inputs.f.operators.charge
    gate_recipe = draft(recipe)
    gate_recipe.overrides[amplitude] = 0.
    drive = Interaction(:charge, amplitude*gate_recipe.inputs.f.operators.charge)
    drivepulse = Pulse(Constant(;duration=2);duration=2,amplitude=0.03)
    microwave = GateRecipe(gate_recipe;controls=[Control(:drive,amplitude,drivepulse)],terms=[drive])
    driven = build(device,microwave)
    @test driven.model.spectrum ≈ prepared.model.spectrum
    @test driven.model.eigensystem.tracking.keys == prepared.model.eigensystem.tracking.keys
    operator = driven.model.operators[QuantumDevices.physical_ref(only(QuantumDevices.refs(n,:op)))]
    @test make_hamiltonian(driven)(nothing,1) ≈ make_hamiltonian(driven.model)+0.03operator
    @test_throws ArgumentError check_endpoints(driven)

    # A complete product must be projected before truncation, including at nested nodes.
    square = Interaction(:square,(amplitude)*n^2)
    nested = Recipe(device,[transmon];inputs=(group=recipe,),states=Dict(:group=>[0,2]))
    nested.overrides[amplitude] = 0.
    @test_throws ArgumentError build(device,GateRecipe(nested;terms=[square],duration=2))
    nested_n = nested.inputs.group.operators.charge
    square = Interaction(:square,(amplitude)*nested_n^2)
    squaregate = build(device,GateRecipe(nested;terms=[square],duration=2))
    ref = QuantumDevices.physical_ref(only(QuantumDevices.refs(n,:op)))
    product = squaregate.model.operator_products[((ref,false),(ref,false))]
    @test norm(product-squaregate.model.operators[ref]^2) > 1e-4
    basis_map(m) = reduce(+, (ket*QT.basis(length(m.spectrum),j-1)' for (j,ket) in enumerate(m.states)))
    leaf = build(device,Leaf(component=:f);parameters=Dict(flux=>0.2))
    intermediate = build(device,recipe)
    embedding = basis_map(leaf)*QuantumDevices.selection(11,1:4)*basis_map(intermediate)*QuantumDevices.selection(4,[1,3])*basis_map(squaregate.model)
    rawcharge = transmon.operator_generators[string(ref.name)](11)
    @test product ≈ embedding'*rawcharge^2*embedding
    leafdrive = Interaction(:leafdrive,(amplitude)*op(:f,:charge))
    @test make_hamiltonian(build(device,GateRecipe(Leaf(component=:f,parameters=Dict(amplitude=>0.));terms=[leafdrive],duration=1))) isa QT.QuantumObject
    @test make_hamiltonian(build(device,GateRecipe(nested;terms=[square],duration=2);parameters=Dict(amplitude=>0.2))) ≈ make_hamiltonian(squaregate.model)+0.2product

    # Captured working data and custom-envelope arrays are independent.
    input_values=[0.,1.,0.]
    samples = InterpolatedEnvelope([0.,1.,2.],input_values;duration=2.)
    mutablepulse = Pulse(samples;duration=2,offset=0.2)
    working = GateRecipe(recipe;controls=[Control(:flux,flux,mutablepulse)],metadata=(notes=["original"],))
    saved = build(device,working)
    before = make_hamiltonian(saved)(nothing,1)
    input_values[2] = 3
    @test samples.values[2]==1.
    @test setpath(samples,"values/2",3.).values[2]==3.
    recipe.overrides[flux] = 0.1
    working.metadata.notes[1] = "edited"
    @test make_hamiltonian(saved)(nothing,1) ≈ before
    @test saved.recipe.metadata.notes[1] == "original"
    editable = draft(saved.recipe)
    @test editable.model_recipe.overrides[flux] == 0.2

    foreign = Interaction(:foreign,(amplitude)*op(:absent,:n))
    @test_throws ArgumentError build(device,GateRecipe(gate_recipe;terms=[foreign],duration=1))
    unknown = Interaction(:unknown,(amplitude)*op(:f,:absent))
    @test_throws ArgumentError build(device,GateRecipe(gate_recipe;terms=[unknown],duration=1))
    conflict = draft(recipe)
    conflict.overrides[flux] = 0.7
    @test_throws ArgumentError make_hamiltonian(prepared.model,GateRecipe(conflict;duration=1))

    complexpulse = Pulse(Constant(;duration=2);duration=2,amplitude=0.03im)
    complexterm = Interaction(:complex,amplitude*gate_recipe.inputs.f.operators.charge+adjoint(amplitude*gate_recipe.inputs.f.operators.charge))
    complexgate = build(device,GateRecipe(gate_recipe;terms=[complexterm],controls=[Control(:complex,amplitude,complexpulse)]))
    @test make_hamiltonian(complexgate)(nothing,1) ≈ make_hamiltonian(complexgate.model)
    sequence = PulseSequence(Pulse(Constant(;duration=1);duration=1,offset=0.2),Pulse(Constant(;duration=2);duration=2,delay=3,offset=0.2);delay=2,offset=0.2)
    both = build(device,GateRecipe(gate_recipe;terms=[drive],controls=[Control(:flux,flux,sequence),Control(:drive,amplitude,drivepulse)],duration=10))
    @test both.recipe.duration == 10
    @test pulse_tstops(both) == [0,2,3,6,8,10]
    @test both.recipe.controls[1].pulse(4) == 0.2
    @test make_hamiltonian(both)(nothing,4) isa QT.QuantumObject
    @test_throws ArgumentError make_hamiltonian(both;herm_check=-1)
end
