@testset "Batched state evolution and gate matrices" begin
    states = [basis(3, 0), basis(3, 1)]
    H = QuantumObject(ComplexF64[0.1 0.05im 0.02; -0.05im 0.2 0.03; 0.02 0.03 -0.1])
    times = [0.0, 0.1, 0.2]
    batch = sesolve(2pi * H, states, times; progress_bar=false, abstol=1e-10, reltol=1e-10)
    @test batch isa QuantumToolbox.TimeEvolutionSol
    @test size(last(batch.states)) == (3, 2)
    @test batch.times_states == times
    for (j, state) in enumerate(states)
        single = sesolve(2pi * H, state, times; progress_bar=false, abstol=1e-10, reltol=1e-10)
        @test last(batch.states).data[:, j] ≈ last(single.states).data atol=1e-8
    end
    full = get_unitary(H, last(times); abstol=1e-10, reltol=1e-10)
    @test last(batch.states).data ≈ full.data[:, 1:2] atol=1e-8
    matrix = get_gate_matrix(batch, states)
    @test matrix ≈ full.data[1:2, 1:2] atol=1e-8
    @test get_gate_matrix(batch, reverse(states)) ≈ reverse(matrix; dims=1)
    @test get_gate_matrix(batch, states; frame=:interaction, idle_hamiltonian=H) ≈
        Matrix{ComplexF64}(I, 2, 2) atol=1e-8

    problem = sesolveProblem(2pi * H, states, times; progress_bar=false)
    @test size(problem.prob.u0) == (3, 2)
    @test sesolve(problem).states[end].data ≈ last(batch.states).data atol=1e-6
    # The elapsed time begins at tlist[1], even when only the final operator is saved.
    shifted = sesolve(2pi * H, states, [1.0, 1.2]; saveat=[1.2], save_start=false,
        progress_bar=false, abstol=1e-10, reltol=1e-10)
    @test length(shifted.states) == 1
    @test get_gate_matrix(shifted, states; frame=:interaction, idle_hamiltonian=H) ≈
        Matrix{ComplexF64}(I, 2, 2) atol=1e-8
    @test_throws ArgumentError sesolve(H, QuantumObject[], times)
    @test_throws ArgumentError sesolve(H, [qeye(3)], times)
    @test_throws DimensionMismatch sesolve(H, [basis(3, 0), basis(2, 0)], times)
    @test_throws ArgumentError get_gate_matrix(batch, [states[1], states[1]])
    @test_throws ArgumentError get_gate_matrix(batch, states; frame=:unknown)
    @test_throws ArgumentError get_gate_matrix(batch, states; frame=:interaction)
    ket_solution = sesolve(H, first(states), times; progress_bar=false)
    @test_throws ArgumentError get_gate_matrix(ket_solution, states)

    model = make_model([make_qubit("q", 0.0)], val(0), (;))
    gate = DeviceGate((drive=Pulse(Constant(); duration=1.0, amplitude=0.25),),
        param(:drive)*op(:q_x))
    basis_states = [basis(2, 0), basis(2, 1)]
    options = (; abstol=1e-10, reltol=1e-10)
    actual = numerical(model, gate, basis_states; options...)
    @test actual isa QuantumObject
    @test isoper(actual)
    @test actual ≈ -im * sigmax() atol=1e-8
    @test get_gate_matrix(model, gate, basis_states; options...) ≈ actual
    @test numerical(model, gate, basis_states; output_states=reverse(basis_states), options...) ≈
        -im * qeye(2) atol=1e-8
    @test gate_infidelity(model, gate, sigmax(); states=basis_states, evolution_kwargs=options) < 1e-8
    @test gate_infidelity(model, gate, sigmay(); states=basis_states,
        include_phases=false, evolution_kwargs=options) < 1e-8
    @test_throws DimensionMismatch numerical(model, gate, basis_states; output_states=basis_states[1:1])
    @test_throws ArgumentError numerical(model, gate, [2basis_states[1], basis_states[2]])
    @test_throws ArgumentError numerical(model, gate, basis_states;
        saveat=Float64[], save_start=false, save_end=false)

    drift_model = make_model([make_qubit("drift", 0.7)], val(0), (;))
    idle_gate = DeviceGate((;), 0, 0.3)
    idle_lab = numerical(drift_model, idle_gate, basis_states; options...)
    @test idle_lab ≈ get_unitary(drift_model, idle_gate; options...)
    @test numerical(drift_model, idle_gate, basis_states; frame=:interaction, options...) ≈
        qeye(2) atol=1e-8
    @test gate_infidelity(drift_model, idle_gate, qeye(2); frame=:interaction,
        evolution_kwargs=options) < 1e-8

    # A selected subspace has its own dimensions, rather than the full model's.
    larger_model = make_model([make_resonator("r", 0.7, 3)], val(0), (;))
    reduced_states = [basis(3, 0), basis(3, 1)]
    reduced_lab = numerical(larger_model, idle_gate, reduced_states; options...)
    reduced_interaction = numerical(larger_model, idle_gate, reduced_states;
        frame=:interaction, options...)
    @test reduced_lab isa QuantumObject
    @test reduced_lab.dimensions == sigmax().dimensions
    @test reduced_lab ≈ QuantumObject(Diagonal([1, cis(-2pi * 0.7 * idle_gate.duration)])) atol=1e-8
    @test reduced_interaction isa QuantumObject
    @test reduced_interaction ≈ qeye(2) atol=1e-8
    @test unitary_fidelity(qeye(2), reduced_interaction) ≈ 1 atol=1e-8
end

@testset "Coherent and probability-only scores" begin
    X, Y, Z = sigmax(), sigmay(), sigmaz()
    @test unitary_fidelity(X, Y) ≈ 0 atol=1e-14
    @test unitary_fidelity(X, Y; include_phases=false) ≈ 1
    @test unitary_infidelity(X, Y; include_phases=false) ≈ 0
    @test unitary_fidelity(qeye(2), Z; include_phases=false) ≈ 1
    @test unitary_fidelity(qeye(2), Z) ≈ 0 atol=1e-14
    hadamard = (X + Z) / sqrt(2)
    @test unitary_fidelity(hadamard, qeye(2); include_phases=false) ≈ 0.5
    for include_phases in (true, false)
        @test unitary_fidelity(hadamard, sqrt(0.5)*hadamard; include_phases) ≈ 0.5
    end
    partial = ComplexF64[1 0; 0 0]
    @test unitary_fidelity(qeye(2), partial) ≈ 0.25
    @test unitary_fidelity(qeye(2), partial; include_phases=false) ≈ 0.5
    @test unitary_fidelity(qeye(2), zeros(2, 2); include_phases=false) == 0
end

@testset "Prepared calibration and dynamic overrides" begin
    model = make_model([make_qubit("q", 0.2)], val(0), (;))
    pulse = Pulse(SineSquared(); duration=1.0, delay=0.1, amplitude=0.1)
    gate = DeviceGate((drive=pulse,), param(:drive)*op(:q_x))
    states = [basis(2, 0), basis(2, 1)]
    options = (; abstol=1e-9, reltol=1e-9)
    problem = CalibrationProblem(model, gate, ["drive/amplitude", "drive/duration"], sigmax();
        states, frame=:interaction, include_phases=false, evolution_kwargs=options,
        lb=[0.0, 0.5], ub=[0.5, 2.0])
    @test problem.problem.lb == [0.0, 0.5]
    for values in ([0.1, 1.0], [0.25, 1.2])
        candidate = calibrated_gate(problem, values)
        actual = gate_infidelity(model, candidate, sigmax(); states, frame=:interaction,
            include_phases=false, evolution_kwargs=options)
        @test problem.problem.f(values, nothing) ≈ actual atol=1e-12
    end
    @test_throws ArgumentError CalibrationProblem(model, gate, ["drive/amplitude"], sigmax();
        states, objective=(_, _, _) -> 0.0)
    @test_throws ArgumentError CalibrationProblem(model, gate, ["drive/amplitude"], sigmax();
        include_phases=false, objective=(_, _, _) -> 0.0)
    @test_throws DimensionMismatch CalibrationProblem(model, gate, ["drive/amplitude"], qeye(3); states)

    # An idle-parameter override must remain dynamic when constants are grouped.
    flux = Pulse(SineSquared(); duration=1.0, delay=0.1, amplitude=0.3, offset=0.2)
    flux_gate = DeviceGate((q_ν=flux,), 0)
    flux_problem = CalibrationProblem(model, flux_gate, ["q_ν/amplitude"], qeye(2);
        states, evolution_kwargs=options)
    cached = QD._prepared_gate_numerical(model, flux_gate)
    for amplitude in (0.1, 0.3)
        candidate = calibrated_gate(flux_problem, [amplitude])
        H = numerical(model, candidate)
        for t in (0.0, 0.1, 0.4, 1.1)
            @test H(t) ≈ (flux.offset + amplitude *
                (0.1 <= t <= 1.1 ? envelope_value(flux.envelope, t-0.1, 1.0) : 0)) * sigmaz()/2
            @test cached(candidate)(t) ≈ H(t)
        end
        @test flux_problem.problem.f([amplitude], nothing) ≈
            gate_infidelity(model, candidate, qeye(2); states, evolution_kwargs=options)
    end
    # Grouped and original symbolic terms agree without dropping unresolved parameters.
    compiled = QD._projected_terms(gate.hamiltonian, model.components, [2])
    combined = QD._gate_parameters(model, gate)
    old = QobjEvo(tuple(((O, QD.scalar_function(c, combined))
        for (O, c) in vcat(model.compiled_terms, compiled))...))
    new = numerical(model, gate)
    for t in (0.0, 0.15, 0.6, gate.duration)
        @test new(t) ≈ old(t)
    end
    unresolved = numerical(op(:q_x)*param(:missing), model.operators, (;))
    @test unresolved((missing=0.3,), 0.0) ≈ 0.3sigmax()
end

# Exercise the cached, in-place ODE derivative after compilation.
function _gate_rhs_allocations(problem)
    f = problem.prob.f
    u = problem.prob.u0
    du = similar(u)
    p = problem.prob.p
    t = 0.3
    f(du, u, p, t)
    return @allocated for _ in 1:1000
        f(du, u, p, t)
    end
end

@testset "Scaled gate Hamiltonians avoid derivative allocations" begin
    model = make_model([make_qubit("q", 0.7)], val(0), (;))
    pulse = Pulse(Gaussian(0.25); amplitude=0.1, duration=1.0,
        carrier=SineCarrier(0.7; phase=0.2))
    gate = DeviceGate((; drive=pulse), 0.5param(:drive)*op(:q_x))
    states = [basis(2, 0), basis(2, 1)]
    H = numerical(model, gate; scalar=2pi)
    for t in (0.0, 0.3, 1.0)
        @test H(t) ≈ 2pi*numerical(model, gate)(t)
    end
    cached = QD._prepared_gate_numerical(model, gate; scalar=2pi)
    @test cached(gate)(0.3) ≈ H(0.3)
    problem = sesolveProblem(H, states, [0.0, gate.duration]; progress_bar=false)
    _gate_rhs_allocations(problem)
    @test _gate_rhs_allocations(problem) == 0
    options = (; abstol=1e-10, reltol=1e-10)
    direct = sesolve(H, states, [0.0, gate.duration]; progress_bar=false, options...)
    expected = get_gate_matrix(direct, states)
    @test numerical(model, gate, states; options...).data ≈ expected atol=1e-8
    @test get_unitary(model, gate; options...).data ≈ expected atol=1e-8
end

@testset "Automatic energy centering preserves lab-frame solutions" begin
    options = (; progress_bar=false, abstol=1e-10, reltol=1e-10)
    H = QuantumObject(ComplexF64[8.0 0.08im 0.0; -0.08im 8.3 0.02; 0.0 0.02 9.0])
    states = [basis(3, 0), basis(3, 1)]
    times = [1.0, 1.1, 1.2]
    centered = sesolve(H, states, times; options...)
    original = sesolve(H, states, times; energy_shift=0, options...)
    for (a, b) in zip(centered.states, original.states)
        @test a ≈ b atol=1e-8
    end
    @test centered.times == original.times
    @test centered.times_states == original.times_states
    @test centered.retcode == original.retcode
    @test centered.alg == original.alg
    final_only = sesolve(H, states, times; saveat=[last(times)], save_start=false, options...)
    @test only(final_only.states) ≈ last(original.states) atol=1e-8
    custom = sesolve(H, states, times; energy_shift=7.0, options...)
    @test last(custom.states) ≈ last(original.states) atol=1e-8
    @test_throws ArgumentError sesolve(H, states, times; energy_shift=:bad, options...)
    @test_throws ArgumentError sesolve(H, states, times; energy_shift=Inf, options...)
    # Dynamic coefficients receive native params and absolute solver time.
    time_H = (H, (QuantumObject(ComplexF64[0 1 0; 1 0 0; 0 0 0]), (p,t)->p.a*cos(t)))
    dynamic = sesolve(time_H, states, times; params=(; a=0.1), options...)
    baseline = sesolve(time_H, states, times; energy_shift=0, params=(; a=0.1), options...)
    @test last(dynamic.states) ≈ last(baseline.states) atol=1e-8

    model = make_model([make_resonator("r", 0.7, 3)], val(0), (;))
    gate = DeviceGate((; drive=Pulse(Constant(); duration=0.3, amplitude=0.0)),
        param(:drive)*op(:r_a))
    basis_states = [basis(3, 0), basis(3, 1)]
    grid = [0.0, 0.1, 0.3]
    native = sesolve(numerical(model, gate; scalar=2pi), basis_states, grid;
        energy_shift=0, options...)
    device = sesolve(model, gate, basis_states, grid; options...)
    for (a, b) in zip(device.states, native.states)
        @test a ≈ b atol=1e-8
    end
    single = sesolve(model, gate, basis_states[2], grid; options...)
    @test last(single.states).data ≈ last(native.states).data[:,2] atol=1e-8
    @test last(sesolve(model, gate, basis_states[2], grid; energy_shift=0.2, options...).states) ≈
        last(single.states) atol=1e-8
    @test numerical(model, gate, basis_states; energy_shift=0, abstol=1e-10, reltol=1e-10) ≈
        numerical(model, gate, basis_states; abstol=1e-10, reltol=1e-10) atol=1e-8
    @test numerical(model, gate, basis_states; frame=:interaction, abstol=1e-10, reltol=1e-10) ≈
        qeye(2) atol=1e-8
    @test_throws ArgumentError numerical(model, gate, basis_states; energy_shift=NaN)

    # Centering is prepared once per calibration and never changes native solve dispatch.
    for frame in (:lab, :interaction), include_phases in (true, false)
        cp = CalibrationProblem(model, gate, ["drive/amplitude"], qeye(2);
            states=basis_states, frame, include_phases,
            evolution_kwargs=(; abstol=1e-10, reltol=1e-10))
        unshifted = CalibrationProblem(model, gate, ["drive/amplitude"], qeye(2);
            states=basis_states, frame, include_phases, energy_shift=0,
            evolution_kwargs=(; abstol=1e-10, reltol=1e-10))
        @test cp.problem.f([0.0], nothing) ≈ unshifted.problem.f([0.0], nothing) atol=1e-8
    end
    @test_throws ArgumentError CalibrationProblem(model, gate, ["drive/amplitude"], qeye(2);
        states=basis_states, energy_shift=:bad)
end

@testset "Dense gate construction and evolution" begin
    model = make_model([make_resonator("r", 0.7, 4)], val(0), (;))
    pulse = Pulse(Gaussian(0.25); amplitude=0.1, duration=1.0,
        carrier=SineCarrier(0.7; phase=0.2))
    gate = DeviceGate((; drive=pulse), param(:drive)*(op(:r_a)+op(:r_adag)))
    states = [basis(4, 0), basis(4, 1)]
    H = numerical(model, gate; scalar=2pi, dense=true)
    # Check the stored matrices, rather than only QobjEvo's evaluated snapshot.
    for term in H.data.ops
        matrix = hasproperty(term, :L) ? term.L.A : term.A
        @test matrix isa Matrix
    end
    original = numerical(model, gate; scalar=2pi)
    for t in (0.0, 0.3, 1.0)
        @test H(t) ≈ original(t)
    end
    problem = sesolveProblem(H, states, [0.0, gate.duration]; progress_bar=false)
    _gate_rhs_allocations(problem)
    @test _gate_rhs_allocations(problem) == 0
    cached = QD._prepared_gate_numerical(model, gate; scalar=2pi, energy_shift=0.35, dense=true)
    @test cached(gate)(0.3) ≈ H(0.3)-2pi*0.35qeye_like(model.H)
    centered_problem = sesolveProblem(cached(gate), states, [0.0, gate.duration]; progress_bar=false)
    _gate_rhs_allocations(centered_problem)
    @test _gate_rhs_allocations(centered_problem) == 0
    options = (; abstol=1e-10, reltol=1e-10)
    @test numerical(model, gate, states; dense=true, options...) ≈
        numerical(model, gate, states; options...) atol=1e-8
    # A purely static gate returns a dense QuantumObject rather than QobjEvo.
    idle = DeviceGate((;), 0, 1.0)
    @test numerical(model, idle; dense=true).data isa Matrix

    grid = [0.0, 0.3, gate.duration]
    for initial in (states, states[2])
        dense_solution = sesolve(model, gate, initial, grid; dense=true,
            saveat=[0.3, gate.duration], save_start=false, progress_bar=false, options...)
        original_solution = sesolve(model, gate, initial, grid;
            saveat=[0.3, gate.duration], save_start=false, progress_bar=false, options...)
        @test dense_solution.times_states == original_solution.times_states
        for (a, b) in zip(dense_solution.states, original_solution.states)
            @test a ≈ b atol=1e-8
        end
    end
    for frame in (:lab, :interaction), include_phases in (true, false)
        dense_problem = CalibrationProblem(model, gate, ["drive/amplitude"], qeye(2);
            states, frame, include_phases, dense=true, evolution_kwargs=options)
        original_problem = CalibrationProblem(model, gate, ["drive/amplitude"], qeye(2);
            states, frame, include_phases, evolution_kwargs=options)
        for amplitude in (0.05, 0.12)
            @test dense_problem.problem.f([amplitude], nothing) ≈
                original_problem.problem.f([amplitude], nothing) atol=1e-8
        end
    end
    @test gate_infidelity(model, gate, qeye(2); states, dense=true, evolution_kwargs=options) ≈
        gate_infidelity(model, gate, qeye(2); states, evolution_kwargs=options) atol=1e-8
end
