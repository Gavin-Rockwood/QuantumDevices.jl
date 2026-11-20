@testset "gate evolution" begin
    device = QuantumDevice(:evolution)
    component = add_component!(device,make_transmon(:q,0.25,15;dimension=5))
    recipe = Recipe(device,[component];states=Dict(component=>0:2))
    prepared = build(device,GateRecipe(recipe;duration=0.13))
    evolution = evolve(prepared;abstol=1e-11,reltol=1e-11)
    energies = collect(prepared.model.spectrum)
    exact = exp(-2pi*im*0.13*Diagonal(energies))
    @test evolution.times == [0.13]
    @test Matrix(unitary(evolution).data) ≈ exact atol=1e-8
    @test Matrix(gate_matrix(evolution).data) ≈ exact atol=1e-8
    @test norm(exact - Matrix{ComplexF64}(I, 3, 3)) > 0.1
    @test_throws MethodError gate_matrix(evolution; remove_idle=true)
    @test collect(idle_phases(evolution)) ≈ 2pi*0.13*energies
    labels = [(1,),(0,)]
    selected = evolve(prepared;states=labels,saveat=[0.,0.04],abstol=1e-11,reltol=1e-11,dtmax=0.01)
    @test selected.times == [0.,0.04,0.13]
    @test Matrix(gate_matrix(selected).data) ≈ exact[[2,1],[2,1]] atol=1e-8
    @test Matrix(gate_matrix(evolution;states=labels).data) ≈ Matrix(gate_matrix(selected).data) atol=1e-8
    centered = evolve(prepared;states=labels,saveat=[0.,0.04],energy_shift=:auto,
        abstol=1e-11,reltol=1e-11)
    @test centered.solver.energy_shift ≈ sum(energies[[2,1]])/2
    @test centered.times == selected.times
    for step in eachindex(centered.times)
        @test centered.propagators[step].data ≈ selected.propagators[step].data atol=1e-8
    end
    shifted = evolve(prepared;energy_shift=energies[1],abstol=1e-11,reltol=1e-11)
    @test unitary(shifted).data ≈ exact atol=1e-8
    @test gate_matrix(centered).data ≈ exact[[2,1],[2,1]] atol=1e-8
    @test_throws ArgumentError evolve(prepared;energy_shift=Inf)
    @test_throws ArgumentError evolve(prepared;energy_shift=:invalid)
    @test idle_phases(selected)((1,)) ≈ 2pi*0.13*energies[2]
    @test Matrix(unitary(evolve(prepared;saveat=[0.04,0.13]);step=1).data) ≈ exp(-2pi*im*0.04*Diagonal(energies)) atol=1e-5
    kets = [(QT.basis(3,0)+im*QT.basis(3,1))/sqrt(2),QT.basis(3,2)]
    complexresult = evolve(prepared;states=kets,abstol=1e-11,reltol=1e-11)
    @test complexresult.propagators[end].data ≈ exact*hcat((s.data for s in kets)...) atol=1e-8
    outputs = [QT.basis(3,0)]
    projected = gate_matrix(complexresult;output_states=outputs)
    @test size(projected) == (1,2)
    @test sum(abs2,projected.data) ≈ 0.5 atol=1e-8
    @test gate_matrix(complexresult;output_states=outputs).data ≈ outputs[1].data'*exact*hcat((s.data for s in kets)...) atol=1e-8
    @test_throws ArgumentError unitary(selected)
    @test_throws ArgumentError gate_matrix(selected;states=labels)
    @test_throws ArgumentError evolve(prepared;states=[])
    @test_throws ArgumentError evolve(prepared;states=[QT.basis(3,0),QT.basis(3,0)])
    @test_throws DimensionMismatch evolve(prepared;states=[QT.basis(2,0)])
    @test_throws ArgumentError evolve(prepared;saveat=[0.1,0.02])
    @test_throws ArgumentError evolve(prepared;saveat=[NaN])
    @test_throws ArgumentError evolve(prepared;tstops=[1.])
    @test_throws ArgumentError evolve(prepared;save_end=false)
    @test_throws ArgumentError evolve(prepared;save_idxs=[1])
    @test_throws ArgumentError evolve(build(device,GateRecipe(recipe;duration=0)))
    @test_throws ErrorException evolve(prepared;maxiters=1,dt=1e-8)
    before = collect(idle_phases(evolution))
    getfield(prepared.model.spectrum,:_storage)[1] = 100
    @test collect(idle_phases(evolution)) == before
    @test evolution.reference.tracking !== prepared.model.eigensystem.tracking
end

@testset "driven gate evolution and pulse stops" begin
    device = QuantumDevice(:driven_evolution)
    q = add_component!(device,make_qubit(:q;frequency=0.7))
    recipe = Recipe(device,[q])
    amplitude = Symbolics.variable(:amplitude)
    term = Interaction(:drive,(amplitude)*recipe.inputs.q.operators.x;parameters=Dict(amplitude=>0.))
    pulse = Pulse(Constant(;duration=0.2);duration=0.2,amplitude=0.3)
    prepared = build(device,GateRecipe(recipe;terms=[term],controls=[Control(:drive,amplitude,pulse)]))
    H = Matrix(make_hamiltonian(prepared)(nothing,0.1).data)
    result = evolve(prepared;abstol=1e-11,reltol=1e-11)
    @test unitary(result).data ≈ exp(-2pi*im*0.2*H) atol=1e-8
    @test result.solver.options.abstol == 1e-11
    @test result.solver.options.tstops == pulse_tstops(prepared)
    # Noncommuting idle and driven intervals exercise delayed onset and gaps.
    sequence = PulseSequence(Pulse(Constant(;duration=0.03);duration=0.03,amplitude=0.3,delay=0.8),
        Pulse(Constant(;duration=0.04);duration=0.04,amplitude=-0.2,delay=0.5),
        Pulse(Constant(;duration=0.02);duration=0.02,amplitude=0.1);delay=0.2)
    prepared = build(device,GateRecipe(recipe;terms=[term],controls=[Control(:drive,amplitude,sequence)],duration=2.))
    result = evolve(prepared;tstops=[0.5],abstol=1e-10,reltol=1e-10)
    stops = pulse_tstops(prepared)
    exact = Matrix{ComplexF64}(I,2,2)
    for (a,b) in zip(stops[1:end-1],stops[2:end])
        exact = exp(-2pi*im*(b-a)*Matrix(make_hamiltonian(prepared)(nothing,(a+b)/2).data))*exact
    end
    @test unitary(result).data ≈ exact atol=2e-7
    centered = evolve(prepared;tstops=[0.5],energy_shift=:auto,abstol=1e-10,reltol=1e-10)
    @test unitary(centered).data ≈ exact atol=2e-7
    @test centered.solver.options.tstops == result.solver.options.tstops
    @test 0.5 in result.solver.options.tstops
    @test all(t -> t in result.solver.options.tstops,stops)
    @test evolve(prepared;saveat=[0,1]).times == [0.,1.,2.]
    @test_throws BoundsError unitary(result;step=2)
    @test_throws ArgumentError evolve(prepared;states=[2QT.basis(2,0)])
end

@testset "interacting reference retains physical conditional phase" begin
    device = QuantumDevice(:physical_reference_phase)
    q1 = add_component!(device, make_qubit(:q1; frequency=0.7))
    q2 = add_component!(device, make_qubit(:q2; frequency=1.1))
    recipe = Recipe(device, [q1,q2]; interactions=[
        Interaction(:conditional, 0.17*q1.operators["z"]*q2.operators["z"])])
    prepared = build(device, GateRecipe(recipe; duration=0.23))
    labels = [(0,0),(0,1),(1,0),(1,1)]
    result = evolve(prepared; states=labels, energy_shift=:auto, abstol=1e-11, reltol=1e-11)
    energies = [stateenergy(prepared.model,label) for label in labels]
    expected = exp.(-2pi*im*0.23.*energies)
    physical = gate_matrix(result).data
    @test physical ≈ Diagonal(expected) atol=1e-8
    phase = angle(physical[4,4]*physical[1,1]*conj(physical[3,3]*physical[2,2]))
    @test phase ≈ angle(expected[4]*expected[1]*conj(expected[3]*expected[2])) atol=1e-8
    @test abs(phase) > 0.1
    @test gate_infidelity(Matrix{ComplexF64}(I,4,4),physical) > 0.01
end
