@testset "construction basis and explicit eigenstate projections" begin
    d = QuantumDevice(:construction_basis)
    q = add_component!(d,make_transmon(:q,.25,15.;dimension=5))
    r = add_component!(d,make_transmon(:r,.28,17.;dimension=7))
    leaf = build(d,Leaf(component=:q))
    generators = Dict(ref=>generator(5) for (ref,generator) in QuantumDevices.component_operator_generators(q))
    rawH = evaluate(q.hamiltonian,leaf.resolved_parameters,generators)
    @test leaf.H ≈ rawH
    @test make_hamiltonian(leaf) ≈ rawH
    @test leaf.operators[OperatorRef(:q,:charge)] ≈ generators[OperatorRef(:q,:charge)]
    @test !isdiag(Matrix(leaf.H.data))
    for (energy,ket) in zip(leaf.energies,leaf.states)
        @test leaf.H*ket ≈ energy*ket
    end
    @test states(leaf) == collect(leaf.states)
    @test state(leaf,1) === leaf.states(1)
    left = Recipe(d,[q];states=Dict(q=>0:2))
    right = Recipe(d,[r];states=Dict(r=>0:2))
    coupled = Recipe((left=left,right=right);interactions=[Interaction(:coupling,
        .03*q.operators["charge"]*r.operators["charge"])])
    amplitude = Symbolics.variable(:basis_drive)
    gate_model = draft(coupled)
    gate_model.overrides[amplitude] = 0.
    charge = gate_model.inputs.left.operators.charge
    term = Interaction(:square,amplitude*charge^2)
    prepared = build(d,GateRecipe(gate_model;terms=[term],duration=.05))
    model = prepared.model
    lm,rm = build(d,left),build(d,right)
    P1,P2 = input_basis_map(lm),input_basis_map(rm)
    n1,n2 = lm.operators[OperatorRef(:q,:charge)],rm.operators[OperatorRef(:r,:charge)]
    expectedcharge = QT.tensor(P1'*n1*P1,QT.qeye(3))
    expectedH = QT.tensor(QuantumDevices.spectrum_operator(lm.spectrum),QT.qeye(3))+
        QT.tensor(QT.qeye(3),QuantumDevices.spectrum_operator(rm.spectrum))+
        .03*QT.tensor(P1'*n1*P1,P2'*n2*P2)
    @test model.H ≈ expectedH
    @test !isdiag(Matrix(model.H.data))
    @test model.operators[OperatorRef(:q,:charge)] ≈ expectedcharge
    # Complete products cross projections as products in upstream coordinates.
    ref = OperatorRef(:q,:charge)
    leafprojector = QuantumDevices.projector_from_states(collect(leaf.states[1:3]))
    leftsquare = leafprojector' * generators[ref]^2 * leafprojector
    expected_square = QT.tensor(P1'*leftsquare*P1,QT.qeye(3))
    @test model.operator_products[((ref,false),(ref,false))] ≈ expected_square
    @test norm(expected_square-expectedcharge^2)>1e-5
    for key in model.eigensystem.keys
        ket = state(model,key)
        @test ket === model.states(key)
        @test model.H*ket ≈ stateenergy(model,key)*ket
        @test model.states[stateindex(model,key)+1] === ket
    end
    labels = [(1,0),(0,0),(0,1)]
    reduction = Recipe(coupled;states=labels,
        dressing_options=(trajectory=x->error("pure reduction cannot run a tracking ramp"),))
    final = build(d,reduction)
    indices = sort([findfirst(==(key),model.eigensystem.keys) for key in labels])
    projector = QuantumDevices.projector_from_states([model.states[i] for i in indices])
    @test final.H ≈ projector'*model.H*projector
    @test final.H ≈ QuantumDevices.spectrum_operator(final.spectrum)
    @test final.operators[ref] ≈ projector'*model.operators[ref]*projector
    @test all(isapprox.(states(final),[QT.basis(3,i) for i in 0:2]))
    @test final.labeling.total_steps == model.labeling.total_steps
    @test final.labeling.diagnostics[:] == model.labeling.diagnostics[indices]
end

@testset "construction-coordinate gate states, idle phases and calibration bases" begin
    d = QuantumDevice(:basis_evolution)
    a = add_component!(d,make_qubit(:a;frequency=.7))
    b = add_component!(d,make_qubit(:b;frequency=1.1))
    recipe = Recipe(d,[a,b];interactions=[Interaction(:coupled,.2*a.operators["x"]*b.operators["x"])])
    prepared = build(d,GateRecipe(recipe;duration=.07))
    model = prepared.model
    @test make_hamiltonian(model,prepared.recipe) ≈ make_hamiltonian(prepared)
    H = Matrix(model.H.data)
    exact = exp(-2pi*im*.07*H)
    result = evolve(prepared;abstol=1e-11,reltol=1e-11)
    @test unitary(result).data ≈ exact atol=1e-8
    @test unitary(result).dims == model.H.dims
    U = input_basis_map(model).data
    @test gate_matrix(result).data ≈ U'*exact*U atol=1e-8
    @test norm(gate_matrix(result).data - Matrix{ComplexF64}(I,4,4)) > 0.1
    @test result.reference.tracking.states[1] !== model.states[1]
    captured = copy(gate_matrix(result).data)
    original = copy(model.states[1].data)
    model.states[1].data[1] += .2im
    @test gate_matrix(result).data == captured
    model.states[1].data .= original
    labels = [(1,0),(0,0)]
    selected = evolve(prepared;states=labels,energy_shift=:auto,abstol=1e-11,reltol=1e-11)
    kets = [state(model,key) for key in labels]
    columns = hcat((ket.data for ket in kets)...)
    @test selected.propagators[end].data ≈ exact*columns atol=1e-8
    @test selected.solver.energy_shift ≈ sum(stateenergy(model,key) for key in labels)/2
    @test gate_matrix(selected).data ≈ columns'*exact*columns atol=1e-8
    superposition = [(kets[1]+im*kets[2])/sqrt(2)]
    super = evolve(prepared;states=superposition,energy_shift=:auto,abstol=1e-11,reltol=1e-11)
    @test super.solver.energy_shift ≈ real(dot(superposition[1].data,H*superposition[1].data))
    @test gate_matrix(super).data ≈ reshape([dot(superposition[1].data,exact*superposition[1].data)],1,1) atol=1e-8
    # A tiny objective evaluation, without an optimizer or calibration scan.
    target = Matrix{ComplexF64}(I,4,4)
    drive = Symbolics.variable(:construction_drive)
    gate_model = draft(recipe)
    gate_model.overrides[drive] = 0.
    term = Interaction(:drive,drive*gate_model.inputs.a.operators.x)
    driven = build(d,GateRecipe(gate_model;terms=[term],controls=[Control(:drive,drive,Pulse(Constant(;duration=.07);duration=.07,amplitude=0.))]))
    problem = CalibrationProblem(driven,["drive/amplitude"],target;
        evolution_settings=(;abstol=1e-10,reltol=1e-10))
    @test all(isapprox.(problem.states,states(driven.model)))
    expected_loss = gate_infidelity(target, U'*exact*U)
    @test problem.problem.f([0.],nothing) ≈ expected_loss atol=1e-8
    @test expected_loss > 1e-3
    @test_throws ArgumentError CalibrationProblem(driven,["drive/amplitude"],target;remove_idle=true)
end
